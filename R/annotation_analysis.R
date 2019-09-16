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
#' @param verbose logical indicating whether to print messages
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
#' @examples
#' get_prior_ratios(umpra_example, marg_prior_example, cond_prior_example)
#' @export
get_prior_ratios = function(mpra_data,
                            marg_prior,
                            cond_prior,
                            n_cores = 1,
                            verbose = TRUE) {

  sample_depths = get_sample_depths(mpra_data)

  if (verbose) {
    message('Determining well-represented variants, see plot...')
  }

  well_represented = get_well_represented(mpra_data,
                                          sample_depths,
                                          rep_cutoff = .15,
                                          plot_rep_cutoff = FALSE,
                                          verbose = verbose)


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
    unnest(c(.data$variant_m_prior)) %>% # this leaves a bunch of messy column names
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
    unnest(c(.data$variant_p_prior)) %>% # this leaves a bunch of messy column names
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

  size_ratios = tibble(sample_id = unique(mu_ratios$sample_id),
                       size_guesses = list(across_samples_size_ratios)) %>%
    unnest(c(.data$size_guesses))



  prior_ratios = bind_rows(mu_ratios, size_ratios)
  return(prior_ratios)
}

#' Get Kullback-Leibler Divergences
#'
#' @description Assess the prior improvement using Kullback-Leibler divergence
#' @param mpra_data a data frame of mpra data
#' @param cond_prior a conditional prior (see fit_cond_prior())
#' @param marg_prior a marginal prior (see fit_marg_prior())
#' @param n_cores number of cores for parallelization
#' @param verbose logical indicating whether to print messages
#' @return a data frame giving the KL between the observed normalized counts and
#'   the two priors for each allele of each variant_id
#' @note The KL divergences are approximations because gamma kernel density
#'   estimation is used to obtain the empirical distribution of the normalized
#'   counts.
get_kl_divergences = function(mpra_data,
                              cond_prior,
                              marg_prior,
                              n_cores,
                              verbose = TRUE) {

  sample_depths = get_sample_depths(mpra_data)

  if (verbose) {
    message('Determining well-represented variants, see plot...')
  }

  well_represented = get_well_represented(mpra_data,
                                          sample_depths,
                                          rep_cutoff = .15,
                                          plot_rep_cutoff = FALSE,
                                          verbose = verbose)


  #### First get the ratios for the mean parameters ----
  mean_dna_abundance = mpra_data %>%
    select(.data$variant_id, .data$allele, .data$barcode, matches('DNA')) %>%
    gather('sample_id', 'counts', matches('DNA|RNA')) %>%
    left_join(sample_depths,
              by = 'sample_id') %>%
    mutate(depth_adj_count = .data$counts / .data$depth_factor) %>%
    group_by(.data$barcode) %>%
    summarise(mean_depth_adj_count = mean(.data$depth_adj_count))

  count_remnants = mpra_data %>%
    filter(.data$barcode %in% well_represented$barcode) %>%
    select(-matches('DNA')) %>%
    gather('sample_id', 'counts', matches('RNA')) %>%
    left_join(sample_depths, by = 'sample_id') %>%
    left_join(mean_dna_abundance, by = 'barcode') %>%
    mutate(count_remnant = .1 + .data$counts / .data$depth_factor / .data$mean_depth_adj_count)

  if (verbose){
    message('Computing conditional prior KL divergences...')
  }

  cond_kl = cond_prior$rna_priors %>%
    select(.data$variant_id, .data$variant_m_prior) %>% # Only implementing for mu_prior for now # TODO reapply to dispersions
    unnest(c(.data$variant_m_prior)) %>%
    right_join(count_remnants,
               by = c('variant_id', 'allele')) %>%
    group_by(.data$variant_id,
             .data$allele) %>%
    nest() %>%
    dplyr::rename('remnants' = 'data') %>%
    ungroup() %>%
    mutate(cond_kl = unlist(parallel::mclapply(.data$remnants,
                                     compute_kl_est,
                                     mc.cores = n_cores))) %>%
    select(-.data$remnants)

  if (verbose) {
    message('Computing marginal prior KL divergences...')
  }

  marg_kl = marg_prior %>%
    filter(.data$prior_type == 'mu_prior', .data$acid_type == 'RNA') %>%
    right_join(count_remnants,
               by = c('allele')) %>%
    group_by(.data$variant_id,
             .data$allele) %>%
    nest() %>%
    dplyr::rename('remnants' = 'data') %>%
    ungroup %>%
    mutate(marg_kl = unlist(parallel::mclapply(.data$remnants,
                                     compute_kl_est,
                                     mc.cores = n_cores))) %>%
    select(-.data$remnants)

  both_kl = cond_kl %>%
    full_join(marg_kl, by = c('variant_id', 'allele'),
              suffix = c('_cond', '_marg'))

}

#' Compute KL divergence estimate
#'
#' @description Approximate the KL divergence between a set of normalized counts
#'   and a prior. An empirical density
#' @param remnant_set a data frame of count remnants (i.e. post sample/DNA
#'   abundance normalization) with columns specifying the prior
#' @return the Kullback-Leibler divergence between the empirical estimate and
#'   the prior
compute_kl_est = function(remnant_set){
  # KL = int(p*log(p/q), x)
  # p is the empirical distribution of count remnants
  # q is the prior
  # Since the count remnants are individual observations, the empirical distribution p is a sum of delta functions.
  # These integrate to one, so the continuous integral reduces to a finite sum.
  # I think.
  # This function will be called with dplyr::do()

  max_support = 2*max(remnant_set$count_remnant)

  count_dens_estimate_fun = kdensity::kdensity(remnant_set$count_remnant,
                                               kernel = 'gamma',
                                               support = c(0, max_support),
                                               normalized = TRUE,
                                               adjust = 2)

  estimate_support = seq(0, max_support, length.out = 1000)[-1] # exclude 0
  estimate_dens = count_dens_estimate_fun(estimate_support)

  prior_dens = dgamma(estimate_support,
                      shape = remnant_set$alpha_est[1],
                      rate = remnant_set$beta_est[1])

  data_frame(kl = sum(estimate_dens * log(estimate_dens / prior_dens)))

  # remnant_set %>%
  #   mutate(point_contribution = log(1 / dgamma(.data$count_remnant,
  #                                              shape = .data$alpha_est[1],
  #                                              rate = .data$beta_est[1]))) %>%
  #   summarise(kl = sum(.data$point_contribution))
  #
  # # "Am I unbelievably sick?" # edit - No

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
#'   annotation sources, rather than make discrete claims / measurements about
#'   individual annotation sources.
score_annotation = function(mpra_data,
                            annotations,
                            nb_params = NULL,
                            conditional_prior = NULL,
                            marginal_prior = NULL){

}
