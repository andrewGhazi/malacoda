#' Evaluate prior ratios of MLE estimates
#'
#' @description This function evaluates prior ratios of MLE estimates under
#'   user-supplied conditional and marginal prior distributions. If this ratio
#'   is above one (or equivalently if the log is above 0), this implies that the conditional prior improves upon the
#'   marginal prior for that particular parameter.
#'
#' @param mpra_data a data frame of MPRA data
#' @param marg_prior a marginal prior
#' @param cond_prior a conditional prior
#' @param n_cores number of cores to use in parallel
#' @details the inputs \code{marg_prior} and \code{cond_prior} can be taken
#'   directly as the outputs of fit_marg_prior and fit_cond_prior.
#'
#'   This output is returned as the difference of log densities.
#' @return a data frame of MLE mean and dispersion parameters by variant_id,
#'   sample_id, and barcode giving prior ratio for each
#' @note The priors for DNA counts are always the same (since we have no prior
#'   knowledge of DNA counts)
#'
#'   The current maximum likelihood estimates are sub-par and will be improved. Interpret with caution.
#' @export
get_prior_ratios = function(mpra_data,
                            marg_prior,
                            cond_prior,
                            n_cores = 1) {

  sample_depths = get_sample_depths(mpra_data)

  print('Determining well-represented variants, see plot...')
  well_represented = get_well_represented(mpra_data,
                                          sample_depths,
                                          rep_cutoff = .15,
                                          plot_rep_cutoff = FALSE)


  #### First get the ratios for the mean parameters ----
  mean_dna_abundance = mpra_data %>%
    select(.data$variant_id, .data$allele, .data$barcode, matches('DNA')) %>%
    gather('sample_id', 'counts', matches('DNA|RNA')) %>%
    left_join(sample_depths,
              by = 'sample_id') %>%
    mutate(depth_adj_count = .data$counts / .data$depth_factor) %>%
    group_by(.data$barcode) %>%
    summarise(mean_depth_adj_count = mean(.data$depth_adj_count))

  rna_cond_mu_prior = cond_prior$rna_priors %>%
    select(-.data$annotation_weights, -.data$variant_p_prior) %>%
    unnest %>% # this leaves a bunch of messy column names
    select(-.data$prior, -matches('acid_type')) %>%
    select(-.data$prior_type)

  count_remnants = mpra_data %>%
    filter(.data$barcode %in% well_represented$barcode) %>%
    select(-matches('DNA')) %>%
    gather('sample_id', 'counts', matches('RNA')) %>%
    left_join(sample_depths, by = 'sample_id') %>%
    left_join(mean_dna_abundance, by = 'barcode') %>%
    mutate(count_remnant = .1 + .data$counts / .data$depth_factor / .data$mean_depth_adj_count)

  marg_mu = marg_prior %>%
    select(-.data$prior, -.data$acid_type) %>%
    filter(grepl('mu', .data$prior_type))

  mu_ratios = count_remnants %>% # the variability of the count after accounting for depth and DNA input
    left_join(marg_mu,
              by = c('allele')) %>%
    rename(marg_alpha = .data$alpha_est,
           marg_beta = .data$beta_est) %>%
    select(-.data$prior_type) %>%
    left_join(rna_cond_mu_prior,
              by = c('variant_id', 'allele')) %>%
    mutate(cond_dens = parallel::mcmapply(dgamma,
                                          .data$count_remnant, .data$alpha_est, .data$beta_est,
                                          MoreArgs = list(log = TRUE),
                                          mc.cores = n_cores,
                                          SIMPLIFY = TRUE),
           marg_dens = parallel::mcmapply(dgamma,
                                          .data$count_remnant, .data$marg_alpha, .data$marg_beta,
                                          MoreArgs = list(log = TRUE),
                                          mc.cores = n_cores,
                                          SIMPLIFY = TRUE),
           log_prior_ratio = .data$cond_dens - .data$marg_dens,
           param_type = 'mean') %>%
    rename(mle_estimate = .data$count_remnant) %>%
    select(.data$variant_id, .data$allele, .data$barcode, .data$sample_id,
           .data$param_type, .data$mle_estimate, .data$log_prior_ratio, .data$cond_dens, .data$marg_dens, .data$marg_alpha:.data$beta_est)

  #### Repeat for dispersion parameters ----

  size_guesses = mpra_data %>%
    select(.data$variant_id, .data$allele, .data$barcode, matches('RNA')) %>%
    filter(.data$barcode %in% well_represented$barcode) %>%
    gather('sample_id', 'counts', matches('DNA|RNA')) %>%
    group_by(.data$allele, .data$barcode) %>%
    summarise(mean_est = mean(.data$counts),
              var_est = var(.data$counts),
              size_guess = .data$mean_est^2 / (.data$var_est - .data$mean_est)) %>% # this works fine for fitting priors but ideally we'd use an actual MLE function
    filter(.data$size_guess > 0 & is.finite(.data$size_guess)) %>% # negative size guess = var < mean --> underdispersed
    filter(.data$size_guess < quantile(.data$size_guess, probs = .99)) %>%
    ungroup

  rna_cond_phi_prior =  cond_prior$rna_priors %>%
    select(-.data$annotation_weights, -.data$variant_m_prior) %>%
    unnest %>% # this leaves a bunch of messy column names
    select(-.data$prior, -matches('acid_type')) %>%
    select(-.data$prior_type)

  variants_barcodes = mpra_data %>%
    select(.data$variant_id, .data$barcode) %>%
    unique

  marg_phi = marg_prior %>%
    select(-.data$prior, -.data$acid_type) %>%
    filter(grepl('phi', .data$prior_type))

  across_samples_size_ratios = size_guesses %>%
    left_join(variants_barcodes, by = 'barcode') %>%
    left_join(marg_phi,
              by = c('allele')) %>%
    rename(marg_alpha = .data$alpha_est,
           marg_beta = .data$beta_est) %>%
    select(-.data$prior_type) %>%
    left_join(rna_cond_phi_prior,
              by = c('variant_id', 'allele')) %>%
    mutate(cond_dens = parallel::mcmapply(dgamma,
                                          .data$size_guess, .data$alpha_est, .data$beta_est,
                                          MoreArgs = list(log = TRUE),
                                          mc.cores = n_cores,
                                          SIMPLIFY = TRUE),
           marg_dens = parallel::mcmapply(dgamma,
                                          .data$size_guess, .data$marg_alpha, .data$marg_beta,
                                          MoreArgs = list(log = TRUE),
                                          mc.cores = n_cores,
                                          SIMPLIFY = TRUE),
           log_prior_ratio = .data$cond_dens - .data$marg_dens,
           param_type = 'dispersion') %>%
    rename(mle_estimate = .data$size_guess) %>%
    select(.data$variant_id, .data$allele, .data$barcode, .data$param_type, .data$mle_estimate,
           .data$log_prior_ratio, .data$cond_dens,
           .data$marg_dens, .data$marg_alpha:.data$beta_est)

  size_ratios = data_frame(sample_id = unique(mu_ratios$sample_id),
                           size_guesses = list(across_samples_size_ratios)) %>%
    unnest



  prior_ratios = bind_rows(mu_ratios, size_ratios)
  return(prior_ratios)
}

#' Score an annotation source
#'
#' @description This function takes MPRA data and an annotation source (or
#'   alternatively pre-fit priors) and checks how well the conditional prior
#'   improves upon the marginal prior by a prior ratio at the maximum likelihood
#'   estimates. A conditional:marginal ratio > 1 indicates that the conditional
#'   prior does better. A higher fraction of alleles for which this is true
#'   indicates an improved conditional prior.
#' @param mpra_data a data frame of MPRA data
#' @param annotations a data frame containing annotations of the same alleles in mpra_data
#' @param nb_params an optional data frame of pre-fit maximum likelihood estimates for allele parameters
#' @param conditional_prior an optional data frame containing pre-fit conditional priors
#' @param marginal_prior an optional data frame containing a pre-fit marginal prior
#' @return A fraction between 0 and 1 indicating the number of parameter
#'   estimates that have higher density under the conditional prior than under
#'   the marginal prior.
#' @note Note that even with near-perfect annotations or conditional priors, the
#'   fraction of parameter estimate ratios above 1 will not approach 1. It will
#'   instead approach some difficult-to-determine limit defined by the noise in
#'   the assay. Thus this functionality should only be used to COMPARE
#'   annotation sources, rather than make discreate claims / measurements about
#'   individual annotation sources.
score_annotation = function(mpra_data,
                            annotations,
                            nb_params = NULL,
                            conditional_prior = NULL,
                            marginal_prior = NULL){

}
