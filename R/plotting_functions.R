#' Create a posterior beeswarm plot
#'
#' @description Create a combined beeswarm & activity posterior plot
#'
#' @param sampler_result a stanfit object
#' @param activities a dataframe of activities
#' @param color_by_sample optional logical indicating whether to color dots by
#'   sample
#' @return a posterior beeswarm ggplot object
#' @export
posterior_beeswarm = function(sampler_result,
                              activities,
                              color_by_sample = FALSE) {


  violin_dat = sampler_result %>%
    rstan::extract(pars = c('ref_act', 'alt_act')) %>%
    bind_cols %>% select_all(~gsub('_act', '', .)) %>%
    gather(allele, post_sample) %>%
    mutate(allele = factor(allele, levels = c('ref', 'alt')))


  if (color_by_sample){
    print('Coloring dots by sample. Colors should be mixed i.e. identically distributed!')
    bee_plot = ggplot(aes(allele, activity),
                      data = activities) +
      ggbeeswarm::geom_beeswarm(aes(color = sample_id),
                                cex = 1.5)
  } else {
    bee_plot = ggplot(aes(allele, activity),
                      data = activities) +
      ggbeeswarm::geom_beeswarm(cex = 1.5)
  }

  bee_plot +
    geom_violin(data = violin_dat,
                aes(allele, post_sample),
                fill = rgb(0,0,0,0),
                width = .75)

}

#' Get ratio label y
#'
#' @description helper function for plot_prior_ratios
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


  suitable = prior_ratios %>% filter(log_prior_ratio > x_limits[1] & log_prior_ratio < x_limits[2])

  if (nrow(suitable) < nrow(prior_ratios)){
    warning(paste0('Removing ', format(100*(1- nrow(suitable) / nrow(prior_ratios)), digits = 3), '% of parameter estimates from plot x-range. Set x_limits to larger values to avoid'))
  }


  fraction_labels = prior_ratios %>%
    group_by(param_type, sample_id) %>%
    summarise(fraction_improved = format(sum(log_prior_ratio > 0) / n(), digits = 3),
              facet_label = paste('Fraction above 0:\n', fraction_improved, sep = ''),
              y = get_label_y(log_prior_ratio))

  prior_ratios %>%
    ggplot(aes(log_prior_ratio)) +
    geom_histogram(bins = n_bins,
                   aes(y = ..density..)) +
    xlim(c(x_limits[1],x_limits[2])) +
    geom_density() +
    geom_vline(lty = 2, xintercept = 0, color = 'grey35') +
    facet_grid(param_type ~ sample_id, scales = 'free_y') +
    geom_text(data = fraction_labels,
              x = x_limits[1] + .75*diff(x_limits),
              aes(label = facet_label,
                  y = y),
              size = 3) +
    labs(title = 'log(Conditional:Marginal) prior density of maximum-likelihood estimates')

}

#' Create ratio by transcription shift hexbins
#'
#' @param prior_ratios a dataframe of barcode prior_ratios
#' @param activities a dataframe of barcode activities
#' @param y_limits y limits of the plot
#'
#' @details prior_ratios can be produced by get_prior_ratios() and activities
#'   can be produced by compute_activities
#' @note The default y_limits cut off a small fraction of points where one prior
#'   or the other does vastly better. It can be set to larger values to avoid
#'   this behavior.
plot_ratio_hexs = function(prior_ratios,
                           activities,
                           y_limits = c(-1,1)){
  mean_ts = activities %>%
    group_by(variant_id, allele) %>%
    summarise(mean_act = mean(activity)) %>%
    filter(all(c('ref', 'alt') %in% allele)) %>%
    summarise(ts = mean_act[allele == 'alt'] - mean_act[allele == 'ref'])

  prior_ratios %>%
    left_join(mean_ts, by = 'variant_id') %>%
    ggplot(aes(ts, log_prior_ratio)) +
    ylim(y_limits) +
    geom_hex() +
    facet_wrap("sample_id")

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
#' @note Get \code{sample_correlations} from get_mpra_correlations
plot_mpra_correlations = function(sample_correlations){
  sample_correlations %>%
    mutate(corr_label = format(correlation, digits = 3)) %>%
    ggplot(aes(x = sample_1, y = sample_2)) +
    geom_tile(aes(fill = correlation)) +
    geom_text(aes(label = corr_label), color = 'white') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
