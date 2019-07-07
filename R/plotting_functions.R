#' Create a posterior beeswarm plot
#'
#' @description Create a combined beeswarm & activity posterior plot
#'
#' @param sampler_result a stanfit object
#' @param variant_activities a dataframe of activities (only for the same variant as sampler_result)
#' @param color_by_sample optional logical indicating whether to color dots by
#'   sample
#' @return a posterior beeswarm ggplot object
#' @note sampler_result objects are written by fit_mpra_model to the out_dir argument for each variant_id
#' @examples
#' variant_activities = activities_example[activities_example$variant_id == '6_135426558_2-3',]
#' posterior_beeswarm(example_posterior,
#'                    variant_activities)
#' @export
posterior_beeswarm = function(sampler_result,
                              variant_activities,
                              color_by_sample = FALSE) {

  n_samples = dplyr::n_distinct(variant_activities$sample_id)
  n_bc = dplyr::n_distinct(variant_activities$barcode)
  n_sampler_bc = sum(grepl('dna_m', names(sampler_result)))

  if (n_bc > 2*n_sampler_bc){
    # Check to make sure the user didn't input the activities for the entire
    # assay.
    stop('Found more than twice as many barcodes in the activities input
         than the sampler result. Did you forget to restrict variant_activities
         to just this variant? Try this:
         dplyr::filter(variant_activities, variant_id == variant_id_to_plot)')
  }

  violin_dat = sampler_result %>%
    rstan::extract(pars = c('ref_act', 'alt_act')) %>%
    bind_cols %>%
    select_all(~gsub('_act', '', .)) %>%
    gather('allele', 'post_sample') %>%
    mutate(allele = factor(.data$allele, levels = c('ref', 'alt')))


  if (color_by_sample){
    message('Coloring dots by sample. Colors should be mixed i.e. identically distributed!')
    bee_plot = ggplot(aes(x = .data$allele,
                          y = .data$activity),
                      data = variant_activities) +
      ggbeeswarm::geom_beeswarm(aes(color = .data$sample_id),
                                cex = 1.5) +
      labs(color = 'Sample ID',
           x = 'Allele',
           y = 'Activity') +
      scale_color_viridis_d(end = .9) + theme_light()

  } else {
    bee_plot = ggplot(aes(x = .data$allele,
                          y = .data$activity),
                      data = variant_activities) +
      ggbeeswarm::geom_beeswarm(cex = 1.5) +
      labs(color = 'Sample ID',
           x = 'allele',
           y = 'Activity') +
      theme_light()
  }

  bee_plot +
    geom_violin(data = violin_dat,
                aes(.data$allele, .data$post_sample),
                fill = rgb(0,0,0,0),
                width = .6)

}

#' Get ratio label y
#'
#' @description helper function for plot_prior_ratios
#' @param ratio_values a vector of ratio values to get a relative position from
get_label_y = function(ratio_values){
  .3 * max(density(ratio_values)$y)
}

#' Plot prior ratios
#'
#' @description Visualize prior ratios by histogram
#' @param prior_ratios a data frame of prior ratios from \code{get_prior_ratios}
#' @param x_limits a length two vector for the x-limits of the histograms
#' @param n_bins number of bins in the histogram
#' @examples
#' # marg_prior = fit_marg_prior(umpra_example)
#' # cond_prior = fit_cond_prior(mpra_data = umpra_example,
#' #                             annotations = u_deepsea,
#' #                             min_neighbors = 10)
#' # ratios_example = get_prior_ratios(umpra_example, marg_prior, cond_prior)
#' # ratios_example is included as a data object in malacoda to make this example run quickly.
#' plot_prior_ratios(ratios_example, n_bins = 30)
#' @export
plot_prior_ratios = function(prior_ratios,
                             x_limits = c(-1,1),
                             n_bins = 100) {


  suitable = prior_ratios %>%
    filter(.data$log_prior_ratio > x_limits[1] &
             .data$log_prior_ratio < x_limits[2])

  if (nrow(suitable) < nrow(prior_ratios)){
    message(paste0('Removing ', format(100*(1- nrow(suitable) / nrow(prior_ratios)), digits = 3), '% of parameter estimates from plot x-range. Set x_limits to larger values to avoid'))
  }


  fraction_labels = prior_ratios %>%
    group_by(.data$param_type, .data$sample_id) %>%
    summarise(fraction_improved = format(sum(.data$log_prior_ratio > 0) / n(), digits = 3),
              facet_label = paste('Fraction above 0:\n', .data$fraction_improved, sep = ''),
              y = get_label_y(.data$log_prior_ratio))

  prior_ratios %>%
    ggplot(aes(.data$log_prior_ratio)) +
    geom_histogram(bins = n_bins,
                   aes(y = .data$..density..)) +
    xlim(c(x_limits[1],x_limits[2])) +
    geom_density(adjust = 2) +
    geom_vline(lty = 2, xintercept = 0, color = 'grey35') +
    facet_grid(param_type ~ sample_id, scales = 'free_y') +
    geom_text(data = fraction_labels,
              x = x_limits[1] + .75*diff(x_limits),
              aes(label = .data$facet_label,
                  y = .data$y),
              size = 3) +
    labs(title = 'log(Conditional:Marginal) prior density of maximum-likelihood estimates',
         x = 'log Prior Ratio',
         y = 'density')

}

#' Create ratio by transcription shift hexbins
#'
#' @param prior_ratios a dataframe of barcode prior_ratios
#' @param model_result a dataframe of model results
#' @param y_limits y limits of the plot
#'
#' @details prior_ratios can be produced by get_prior_ratios() and model_result
#'   can be obtained from fit_mpra_model()
#' @note The default y_limits cut off a small fraction of points where one prior
#'   or the other does vastly better. It can be set to larger values to avoid
#'   this behavior.
#' marg_prior = fit_marg_prior(umpra_example)
#' cond_prior = fit_cond_prior(mpra_data = umpra_example,
#'                             annotations = u_deepsea,
#'                             min_neighbors = 10)
#' ratio_df = get_prior_ratios(umpra_example, marg_prior, cond_prior)
#' example_variants = c("11_8839229_1-2", "15_75303554_2-3", "1_203652141_2-3")
#'
#' examples_to_evaluate = umpra_example[umpra_example$variant_id %in% example_variants,]
#' example_result = fit_mpra_model(mpra_data = examples_to_evaluate,
#'  priors = marg_prior,
#' vb_pass = FALSE,
#' tot_samp = 100,
#' n_warmup = 10) # Likewise, n_warmup should be >500
#'
#' plot_ratio_hexs(ratio_df, example_result)
#' @export
plot_ratio_hexs = function(prior_ratios,
                           model_result,
                           y_limits = c(-.5, .5)){

  mean_ts = model_result %>%
    select(.data$variant_id, .data$ts_post_mean) %>%
    rename(ts = .data$ts_post_mean)

  prior_ratios %>%
    left_join(mean_ts, by = 'variant_id') %>%
    ggplot(aes(.data$ts, .data$log_prior_ratio)) +
    ylim(y_limits) +
    geom_hex(aes(color = .data$..count..)) +
    guides(color = FALSE) +
    scale_color_viridis_c(trans = 'log') +
    scale_fill_viridis_c(trans = 'log') +
    facet_wrap("sample_id") +
    geom_hline(yintercept = 0,
               lty = 2,
               color = 'black') +
    labs(x = 'Transcription Shift',
         y = 'MLE estimate\nConditional : Marginal prior ratio')

}

#' Plot MPRA sample correlations
#'
#' @param sample_correlations a data frame of pairwise sample correlations
#' @return a ggplot showing the pairwise sample correlations
#' @details This function visualizes the pairwise correlations of count samples
#'   in MPRA data. This can be a useful QC metric - samples should correlate
#'   highly with other samples of the same type (DNA or RNA). If ALL samples
#'   correlate highly with one another, this can indicate DNA contamination in
#'   the RNA libraries.
#' @note Get \code{sample_correlations} from get_sample_correlations
#' @examples
#' example_correlations = get_sample_correlations(umpra_example)
#' plot_mpra_correlations(example_correlations)
#' @export
plot_mpra_correlations = function(sample_correlations){
  sample_correlations %>%
    mutate(corr_label = format(.data$correlation, digits = 3)) %>%
    ggplot(aes(x = .data$sample_1, y = .data$sample_2)) +
    geom_tile(aes(fill = .data$correlation)) +
    geom_text(aes(label = .data$corr_label), color = 'white') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = 'Sample', y = 'Sample',
         fill = 'Correlation')
}

#' Plot DNA representation
#'
#' @description Plot DNA representation histogram
#'
#' @param dna_df a dataframe of depth-normalized DNA abundances
#' @param rep_cutoff a representation cutoff
#'
#' @return a DNA representation ggplot
#' @note see the first command in the source code of get_well_represented() for
#'   how to generate dna_df
plot_dna_representation = function(dna_df, rep_cutoff){

  dna_df %>%
    ggplot(aes(.data$mean_depth_adj_count)) +
    geom_histogram(bins = 40,
                   color = 'black',
                   fill = 'grey50') +
    geom_vline(xintercept = quantile(dna_df$mean_depth_adj_count,
                                     probs = rep_cutoff),
               lty = 2) +
    scale_x_log10() +
    labs(x = 'Mean Depth Adjusted DNA barcode count',
         title = 'DNA barcode abundance and cutoff',
         subtitle = paste0('Using depth-adjusted DNA barcode count cutoff of ',
                           round(quantile(dna_df$mean_depth_adj_count,
                                          probs = rep_cutoff),
                                 digits = 3))) +
    geom_rug(alpha = .01)

}
