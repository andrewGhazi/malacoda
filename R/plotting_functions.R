#' Create a posterior beeswarm plot
#'
#' @description Create a combined beeswarm & activity posterior plot
#'
#' @param sampler_result a stanfit object
#' @param activities a dataframe of activities
#' @param color_by_sample optional logical indicating whether to color dots by
#'   sample
#' @return a posterior beeswarm ggplot object
#' @note sampler_result objects are written by fit_mpra_model to the out_dir argument for each variant_id
#' @export
posterior_beeswarm = function(sampler_result,
                              activities,
                              color_by_sample = FALSE) {


  violin_dat = sampler_result %>%
    rstan::extract(pars = c('ref_act', 'alt_act')) %>%
    bind_cols %>%
    select_all(~gsub('_act', '', .)) %>%
    gather('allele', 'post_sample') %>%
    mutate(allele = factor(.data$allele, levels = c('ref', 'alt')))


  if (color_by_sample){
    print('Coloring dots by sample. Colors should be mixed i.e. identically distributed!')
    bee_plot = ggplot(aes(x = .data$allele,
                          y = .data$activity),
                      data = activities) +
      ggbeeswarm::geom_beeswarm(aes(color = .data$sample_id),
                                cex = 1.5)
  } else {
    bee_plot = ggplot(aes(x = .data$allele,
                          y = .data$activity),
                      data = activities) +
      ggbeeswarm::geom_beeswarm(cex = 1.5)
  }

  bee_plot +
    geom_violin(data = violin_dat,
                aes(.data$allele, .data$post_sample),
                fill = rgb(0,0,0,0),
                width = .75)

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
plot_prior_ratios = function(prior_ratios,
                             x_limits = c(-1,1),
                             n_bins = 100) {


  suitable = prior_ratios %>%
    filter(.data$log_prior_ratio > x_limits[1] &
             .data$log_prior_ratio < x_limits[2])

  if (nrow(suitable) < nrow(prior_ratios)){
    warning(paste0('Removing ', format(100*(1- nrow(suitable) / nrow(prior_ratios)), digits = 3), '% of parameter estimates from plot x-range. Set x_limits to larger values to avoid'))
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
    geom_density() +
    geom_vline(lty = 2, xintercept = 0, color = 'grey35') +
    facet_grid(param_type ~ sample_id, scales = 'free_y') +
    geom_text(data = fraction_labels,
              x = x_limits[1] + .75*diff(x_limits),
              aes(label = .data$facet_label,
                  y = .data$y),
              size = 3) +
    labs(title = 'log(Conditional:Marginal) prior density of maximum-likelihood estimates')

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
plot_ratio_hexs = function(prior_ratios,
                           model_result,
                           y_limits = c(-.5, .5)){

  mean_ts = model_result %>%
    select(variant_id, ts_post_mean) %>%
    rename(ts = ts_post_mean)

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
