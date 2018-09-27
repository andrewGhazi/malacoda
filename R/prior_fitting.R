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

fit_marg_prior = function(mpra_data, n_cores = 1){

  sample_depths = mpra_data %>%
    gather(sample_id, counts, matches('DNA|RNA')) %>%
    group_by(sample_id) %>%
    summarise(depth_factor = sum(counts) / 1e6)

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



}

fit_cond_prior = function(mpra_data, annotations){

  print('Fitting marginal DNA prior...')
}
