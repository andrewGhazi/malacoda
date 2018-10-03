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
    select(variant_id, allele, matches('DNA')) %>%
    gather(sample_id, counts, matches('DNA|RNA')) %>%
    left_join(sample_depths,
              by = 'sample_id') %>%
    mutate(depth_adj_count = counts / depth_factor)


  if (plot_rep_cutoff){
     rep_cutoff_plot = all_dna %>%
      ggplot(aes(depth_adj_count)) +
      geom_histogram(bins = 40,
                     color = 'black',
                     fill = 'grey50') +
      geom_vline(xintercept = quantile(all_dna$depth_adj_count,
                                       probs = rep_cutoff),
                 lty = 2) +
      scale_x_log10() +
      facet_grid(sample_id ~ .) +
      labs(x = 'Depth Adjusted DNA barcode count',
           title = 'DNA barcode abundance and cutoff',
           subtitle = paste0('Using depth-adjusted DNA barcode count cutoff of ',
                             round(quantile(all_dna$depth_adj_count,
                                            probs = rep_cutoff),
                                   digits = 3)))

     print(rep_cutoff_plot)
  }

  well_represented = all_dna %>%
    filter(depth_adj_count > quantile(all_dna$depth_adj_count,
                                      probs = rep_cutoff)) %>%
    select(variant_id) %>%
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
                          rep_cutoff = .2){

  sample_depths = mpra_data %>%
    gather(sample_id, counts, matches('DNA|RNA')) %>%
    group_by(sample_id) %>%
    summarise(depth_factor = sum(counts) / 1e6)

  well_represented = get_representation_cutoff(mpra_data,
                                               sample_depths,
                                               rep_cutoff = rep_cutoff,
                                               plot_rep_cutoff = plot_rep_cutoff)

  print('Fitting marginal DNA prior...')

  dna_nb_fits = mpra_data %>%
    select(variant_id, allele, matches('DNA')) %>%
    gather(sample_id, counts, matches('DNA|RNA')) %>%
    group_by(variant_id, allele, sample_id) %>%
    nest(.key = count_dat) %>%
    mutate(nb_fit = mclapply(count_dat, fit_nb, mc.cores = n_cores),
           converged = map_lgl(nb_fit, ~.x$convergence == 0))

  if (!all(dna_nb_fits$converged)) {
    warning(paste0(sum(!dna_nb_fits$converged),
                   ' out of ',
                   nrow(dna_nb_fits),
                   ' DNA-allele-samples failed to converge when fitting negative binomial parameters. A small fraction failing is acceptable.'))
  }

  dna_nb_fits %<>%
    filter(converged) %>%
    left_join(sample_depths, by = 'sample_id') %>%
    mutate(depth_adj_mu_est = map2_dbl(nb_fit, depth_factor, ~.x$par[1] / .y),
           phi_est = map_dbl(nb_fit, ~.x$par[2]),
           acid_type = factor(stringr::str_extract(sample_id, 'DNA|RNA')))

  dna_gamma_prior = dna_nb_fits %>%
    group_by(allele, sample_id) %>%
    summarise(mu_prior = list(fit_gamma(depth_adj_mu_est)),
              phi_prior = list(fit_gamma(phi_est))) %>%
    ungroup %>%
    gather(prior_type, prior, matches('prior')) %>%
    mutate(alpha_est = map_dbl(prior, ~.x$par[1]),
           beta_est = map_dbl(prior, ~.x$par[2])) # doesn't line up :(

  rna_nb_fits = mpra_data %>%
    select(variant_id, allele, matches('RNA')) %>%
    gather(sample_id, counts, matches('DNA|RNA')) %>%
    group_by(variant_id, allele, sample_id) %>%
    nest(.key = count_dat) %>%
    mutate(nb_fit = mclapply(count_dat, fit_nb, mc.cores = n_cores),
           converged = map_lgl(nb_fit, ~.x$convergence == 0))

  if (!all(rna_nb_fits$converged)) {
    warning(paste0(sum(!dna_nb_fits$converged),
                   ' out of ',
                   nrow(dna_nb_fits),
                   ' DNA-allele-samples failed to converge when fitting negative binomial parameters. A small fraction failing is acceptable.'))
  }

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
                            marginal_prior = NULL)