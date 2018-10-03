#' Generate a distance matrix from a matrix of annotations
#'
#' Given an nxd matrix of variant annotations, produce an nxn distance matrix
#' describing the inter-variant distances in annotation space
#'
#' @param annotation_dat an n x d data frame of annotations
#' @param log_distance a logical indicating to use the log1p of the distances (TRUE) or the raw euclidean distances (FALSE)
#'
#'
#' @importFrom magrittr %>%
generate_distance_matrix = function(annotations,
                                    log_distance = TRUE){

  annotations =
  if (log_distance) {
    annotations %>%
      as.data.frame() %>%
      as.matrix %>%
      scale %>%
      dist %>%
      as.matrix %>%
      log1p
  } else {
    annotations %>%
      as.data.frame() %>%
      as.matrix %>%
      scale %>%
      dist %>%
      as.matrix
  }
}

find_prior_weights = function(annotations){

}

get_representation_cutoff = function(mpra_data,
                                     sample_depths,
                                     rep_cutoff,
                                     plot_rep_cutoff = FALSE){

  all_dna = mpra_data %>%
    select(variant_id, allele, barcode, matches('DNA')) %>%
    gather(sample_id, counts, matches('DNA|RNA')) %>%
    left_join(sample_depths,
              by = 'sample_id') %>%
    mutate(depth_adj_count = counts / depth_factor) %>%
    group_by(barcode) %>%
    summarise(mean_depth_adj_count = mean(depth_adj_count))


  if (plot_rep_cutoff){
     rep_cutoff_plot = all_dna %>%
      ggplot(aes(mean_depth_adj_count)) +
      geom_histogram(bins = 40,
                     color = 'black',
                     fill = 'grey50') +
      geom_vline(xintercept = quantile(all_dna$mean_depth_adj_count,
                                       probs = rep_cutoff),
                 lty = 2) +
      scale_x_log10() +
      labs(x = 'Mean Depth Adjusted DNA barcode count',
           title = 'DNA barcode abundance and cutoff',
           subtitle = paste0('Using depth-adjusted DNA barcode count cutoff of ',
                             round(quantile(all_dna$mean_depth_adj_count,
                                            probs = rep_cutoff),
                                   digits = 3))) +
       geom_rug(alpha = .01)

     print(rep_cutoff_plot)
  }

  well_represented = all_dna %>%
    filter(mean_depth_adj_count > quantile(all_dna$mean_depth_adj_count,
                                      probs = rep_cutoff)) %>%
    select(barcode) %>%
    unique

  return(well_represented)

}

#' Fit a marginal prior
#'
#' @param mpra_data a data frame of mpra data
#' @param n_cores number of cores to parallelize across
#' @param plot_rep_cutoff logical indicating whether to plot the representation cutoff used
#' @param rep_cutoff fraction indicating the depth-adjusted DNA count quantile to use as the cutoff
fit_marg_prior = function(mpra_data,
                          n_cores = 1,
                          plot_rep_cutoff = TRUE,
                          rep_cutoff = .15){

  sample_depths = mpra_data %>%
    gather(sample_id, counts, matches('DNA|RNA')) %>%
    group_by(sample_id) %>%
    summarise(depth_factor = sum(counts) / 1e6)

  print('Determining well-represented variants, see plot...')
  well_represented = get_representation_cutoff(mpra_data,
                                               sample_depths,
                                               rep_cutoff = rep_cutoff,
                                               plot_rep_cutoff = plot_rep_cutoff)


  print('Fitting marginal DNA prior...')

  dna_nb_fits = mpra_data %>%
    filter(barcode %in% well_represented$barcode) %>%
    select(variant_id, allele, matches('DNA')) %>%
    gather(sample_id, counts, matches('DNA|RNA')) %>%
    group_by(variant_id, allele, sample_id) %>%
    nest(.key = count_dat) %>%
    filter(map_lgl(count_dat, ~!all(.x$counts == 0))) %>% # some borderline barcodes are 0 in some samples
    mutate(nb_fit = mclapply(count_dat, fit_nb, mc.cores = n_cores),
           converged = map_lgl(nb_fit, ~.x$convergence == 0))

  if (!all(dna_nb_fits$converged)) {
    warning(paste0(sum(!dna_nb_fits$converged),
                   ' out of ',
                   nrow(dna_nb_fits),
                   ' (', round(sum(!dna_nb_fits$converged) / nrow(dna_nb_fits) * 100, digits = 3), '%)',
                   ' DNA-allele-samples failed to converge when fitting negative binomial parameters. A small fraction (<5%) failing is acceptable.'))
  }

  dna_nb_fits %<>%
    filter(converged) %>%
    left_join(sample_depths, by = 'sample_id') %>%
    mutate(depth_adj_mu_est = map2_dbl(nb_fit, depth_factor, ~.x$par[1] / .y),
           phi_est = map_dbl(nb_fit, ~.x$par[2]),
           acid_type = factor(stringr::str_extract(sample_id, 'DNA|RNA'))) %>%
    filter(phi_est < quantile(phi_est, probs = .995)) # cut out severely underdispersed alleles

  dna_gamma_prior = dna_nb_fits %>%
    summarise(mu_prior = list(fit_gamma(depth_adj_mu_est)),
              phi_prior = list(fit_gamma(phi_est))) %>%
    ungroup %>%
    gather(prior_type, prior, matches('prior')) %>%
    mutate(alpha_est = map_dbl(prior, ~.x$par[1]),
           beta_est = map_dbl(prior, ~.x$par[2]),
           acid_type = 'DNA') # doesn't line up :(

  #### Fit RNA prior ----
  mean_dna_abundance = mpra_data %>%
    select(variant_id, allele, barcode, matches('DNA')) %>%
    gather(sample_id, counts, matches('DNA|RNA')) %>%
    left_join(sample_depths,
              by = 'sample_id') %>%
    mutate(depth_adj_count = counts / depth_factor) %>%
    group_by(barcode) %>%
    summarise(mean_depth_adj_count = mean(depth_adj_count))

  print('Fitting marginal RNA mean priors...')
  rna_m_prior = mpra_data %>%
    select(variant_id, allele, barcode, matches('RNA')) %>%
    filter(barcode %in% well_represented$barcode) %>%
    gather(sample_id, counts, matches('DNA|RNA')) %>%
    left_join(sample_depths , by = 'sample_id') %>%
    left_join(mean_dna_abundance, by = 'barcode') %>%
    mutate(count_remnant = .1 + counts / depth_factor / mean_depth_adj_count) %>%
    group_by(allele) %>%
    summarise(mu_prior = list(fit_gamma(count_remnant))) %>%
    gather(prior_type, prior, matches('prior')) %>%
    mutate(alpha_est = map_dbl(prior, ~.x$par[1]),
           beta_est = map_dbl(prior, ~.x$par[2])) %>%
    mutate(acid_type = 'RNA')

  # the +.1 is to give some non-infinite log-density to 0's. See plot below. It seems to work well.

  print('Fitting marginal RNA dispersion priors...')
  rna_p_prior = mpra_data %>%
    select(variant_id, allele, barcode, matches('RNA')) %>%
    filter(barcode %in% well_represented$barcode) %>%
    gather(sample_id, counts, matches('DNA|RNA')) %>%
    group_by(allele, barcode) %>%
    summarise(mean_est = mean(counts),
              var_est = var(counts),
              size_guess = mean_est^2 / (var_est - mean_est)) %>%
    filter(size_guess > 0 & is.finite(size_guess)) %>% # negative size guess = var < mean = underdispersed
    filter(size_guess < quantile(size_guess, probs = .99)) %>% # HUGE size guess = underdispersed, cut out barcodes that are TOO consistent i.e. underdispersed
    summarise(phi_prior = list(fit_gamma(size_guess))) %>%
    gather(prior_type, prior, matches('prior')) %>%
    mutate(alpha_est = map_dbl(prior, ~.x$par[1]),
           beta_est = map_dbl(prior, ~.x$par[2])) %>%
    mutate(acid_type = 'RNA')

  # There is room for improvement with this phi prior estimation. It disregards
  # sequencing sample depth V and cuts out observed larged phi values.
  # It still covers the vast majority of alleles (98%)

  tot_prior = bind_rows(dna_gamma_prior,
                        rna_m_prior,
                        rna_p_prior)

  return(tot_prior)

  # mpra_data %>%
  #   select(variant_id, allele, barcode, matches('RNA')) %>%
  #   gather(sample_id, counts, matches('DNA|RNA')) %>%
  #   left_join(sample_depths , by = 'sample_id') %>%
  #   left_join(mean_dna_abundance, by = 'barcode') %>%
  #   mutate(count_remnant = .1 + counts / depth_factor / mean_depth_adj_count) %>%
  #   ggplot(aes(count_remnant)) +
  #   geom_density(aes(color = sample_id)) +
  #   scale_x_log10() +
  #   stat_function(fun = dgamma, args = list(shape = 1.04, rate = 1.128))

  # ^ per-sample priors might be worthwhile

}

fit_cond_prior = function(mpra_data, annotations){

  print('Fitting marginal DNA prior...')
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
