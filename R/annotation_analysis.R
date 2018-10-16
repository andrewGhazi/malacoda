#' Evaluate prior ratios of MLE estimates
#'
#' @description This function evaluates prior ratios of MLE estimates under
#'   user-supplied conditional and marginal prior distributions. If this ratio
#'   is above one, this implies that the conditional prior improves upon the
#'   marginal prior for that particular parameter.
#'
#' @param mpra_data a data frame of MPRA data
#' @param marg_prior a marginal prior
#' @param cond_prior a conditional prior
#' @details the inputs \code{marg_prior} and \code{cond_prior} can be taken
#'   directly as the outputs of fit_marg_prior and fit_cond_prior
#' @return a data frame of MLE mean and dispersion parameters by variant_id,
#'   sample_id, and barcode giving prior ratio for each
#' @note The priors for DNA counts are always the same (since we have no prior
#'   knowledge of DNA counts)
#' @export
get_prior_ratios = function(mpra_data,
                            marg_prior,
                            cond_prior) {

  sample_depths = get_sample_depths(mpra_data)

  print('Determining well-represented variants, see plot...')
  well_represented = get_well_represented(mpra_data,
                                          sample_depths,
                                          rep_cutoff = rep_cutoff,
                                          plot_rep_cutoff = plot_rep_cutoff)


  mean_dna_abundance = mpra_data %>%
    select(variant_id, allele, barcode, matches('DNA')) %>%
    gather(sample_id, counts, matches('DNA|RNA')) %>%
    left_join(sample_depths,
              by = 'sample_id') %>%
    mutate(depth_adj_count = counts / depth_factor) %>%
    group_by(barcode) %>%
    summarise(mean_depth_adj_count = mean(depth_adj_count))


  mpra_data %>%
    filter(barcode %in% well_represented$barcode) %>%
    select(-matches('DNA')) %>%
    left_join(annotation_weights, by = 'variant_id') %>%
    gather(sample_id, counts, matches('RNA')) %>%
    left_join(sample_depths, by = 'sample_id') %>%
    left_join(mean_dna_abundance, by = 'barcode') %>%
    mutate(count_remnant = .1 + counts / depth_factor / mean_depth_adj_count)

}
