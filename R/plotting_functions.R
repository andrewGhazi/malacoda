#' Create a posterior beeswarm plot
#'
#' @description Create a combined beeswarm & activity posterior plot
#'
#' @param sampler_result a stanfit object
#' @param variant_activities a dataframe of activities (only for the same variant as sampler_result)
#' @param color_by_sample optional logical indicating whether to color dots by
#'   sample
#' @param adjust violin bandwidth adjustment. See ?ggplot2::geom_violin
#' @param verbose logical indicating whether to print messages
#' @return a posterior beeswarm ggplot object
#' @note sampler_result objects are written by fit_mpra_model to the out_dir argument for each variant_id
#' @examples
#' variant_activities = activities_example[activities_example$variant_id == '6_135426558_2-3',]
#' posterior_beeswarm(example_posterior,
#'                    variant_activities)
#' @export
posterior_beeswarm = function(sampler_result,
                              variant_activities,
                              color_by_sample = FALSE,
                              verbose = TRUE,
                              adjust = 1.5) {

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
    if (verbose) {message('Coloring dots by sample. Colors should be mixed i.e. identically distributed!')}

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
                width = .6,
                adjust = adjust)

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
#' @examples
#' plot_ratio_hexs(ratios_example, example_result)
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

format_parameter_name = function(par_name){
  if(grepl('[0-9]+', par_name)){
    return(str_extract(par_name, '[0-9]+'))
  } else{
    return(par_name)
  }
}

#' Plot a variant's posterior intervals
#'
#' @description Given a result object from the malacoda model, this function
#'   plots the posterior intervals for all the modeled parameters. The outer
#'   black bars represent the limits of the 95% interval, while the inner blue
#'   bar represents the 80% interval (by default). The black dots are the
#'   posterior means.
#'
#' @param sampler_res A stanfit object for one variant
#' @param inner_probs quantiles of the posterior for the inner bar
#' @param outer_probs quantiles of the posterior for the outer bar
#' @details The sampler_res object for a given variant is saved into the out_dir
#'   argument of fit_mpra_model. The barcode concentrations (i.e. DNA means) are
#'   numbered in the same order in which they are given to fit_mpra_model, with
#'   poorly represented barcodes removed entirely.
#' @seealso fit_mpra_model get_well_represented
#' @examples malacoda_intervals(example_posterior)
#' @export
malacoda_intervals = function(sampler_res,
                              inner_probs = c(.1, .9),
                              outer_probs = c(.025, .975)){

summary_df =   summary(sampler_res,
                       probs = c(outer_probs[1],
                                 inner_probs[1],
                                 inner_probs[2],
                                 outer_probs[2])) %>%
  .$summary %>%
  as.data.frame %>%
  rownames_to_column('parameter') %>%
  as_tibble

new_names = names(summary_df)
new_names[5:8] = c('outer_left', 'inner_left', 'inner_right', 'outer_right')

summary_df %<>%
    set_names(new_names) %>%
    filter(.data$parameter != 'lp__') %>%
    mutate(par_type = str_remove(.data$parameter, pattern = '\\[[0-9]+\\]|ref_|alt_'),
           parameter = map_chr(.data$parameter, format_parameter_name)) %>%
    mutate(parameter = case_when(grepl('rna_m|rna_p', .data$par_type) & .data$parameter == '1' ~ 'ref',
                                 grepl('rna_m|rna_p', .data$par_type) & .data$parameter == '2' ~ 'alt',
                                 .data$parameter == 'ref_act' ~ 'ref',
                                 .data$parameter == 'alt_act' ~ 'alt',
                                 .data$par_type == 'transcription_shift' ~ 'shift',
                                 .data$par_type == 'dna_p' ~ 'phi',
                                 TRUE ~ as.character(.data$parameter))) %>%
    mutate(parameter = factor(.data$parameter,
                              levels = rev(unique(.data$parameter))),
           par_type = factor(case_when(.data$par_type == 'act' ~ 'Activities',
                                       .data$par_type == 'dna_m_alt' ~ 'Alternate Allele DNA means',
                                       .data$par_type == 'dna_m_ref' ~ 'Reference Allele DNA means',
                                       .data$par_type == 'dna_p' ~ 'DNA dispersion',
                                       .data$par_type == 'rna_m' ~ 'RNA means',
                                       .data$par_type == 'rna_p' ~ 'RNA dispersions',
                                       .data$par_type == 'transcription_shift' ~ 'Transcription Shift'),
                             levels = c('Reference Allele DNA means',
                                        'Alternate Allele DNA means',
                                        'DNA dispersion',
                                        'RNA means',
                                        'RNA dispersions',
                                        'Activities',
                                        'Transcription Shift')))

summary_df %>%
    ggplot(aes(.data$mean, .data$parameter)) +
    geom_segment(aes(x = .data$outer_left, xend = .data$outer_right, yend = .data$parameter),
                 color = 'black',
                 size = 0.5, lineend = 'butt') +
    geom_segment(aes(x = .data$inner_left, xend = .data$inner_right, yend = .data$parameter),
                 color = 'deepskyblue2',
                 size = 1.5, lineend = 'butt') +
    geom_point(size = 1) +
    facet_wrap('par_type', scales = 'free') +
    theme_light() +
    labs(x = 'parameter value',
         y = NULL)  +
    geom_vline(xintercept = 0,
               lty = 2) +
    theme(strip.text = element_text(color = 'grey10', size = 9))
}


#' Plot prior samples
#' @description Visualize a malacoda prior on transcription shift
#' @examples
#' prior_draws = sample_from_prior(marg_prior_example, n_samp = 2000)
#' plot_prior_samples(prior_draws)
#' @return a ggplot visualizing the prior density on TS
#' @param prior_draws a data frame of prior samples
#' @seealso \code{\link{sample_from_prior}}
#' @export
plot_prior_samples = function(prior_draws){
  ggplot(prior_draws,
         aes(.data$sim_ts)) +
    geom_histogram(aes(y = .data$..density..)) +
    labs(x = 'Transcription Shift',
         y = 'Prior density')

}

#' Make a MPRA tile plot
#' @description Make a tile plot of the MPRA counts for a given variant
#' @param variant_id variant ID
#' @param variant_counts a data frame of counts for one variant
#' @param show_barcodes a logical indicating whether to show the barcodes on the
#'   left side of the tile plot
#' @param composite a logicial indicating whether or not to include the extra
#'   plot components on the right sides
#' @param sampler_res A stanfit object for one variant
#' @param prob if composite = TRUE, a fraction between 0 and 1 indicating the
#'   posterior probability mass to include in the inner intervals.
#' @param prob_outer if composite = TRUE, a fraction between 0 and 1 indicating
#'   the posterior probability mass to include in the outer intervals.
#' @details The barcodes in the output tile plot are sorted according the their
#'   average DNA representation.
#'
#'   The barcodes in variant_counts need to be in the same order that they were
#'   presented to the malacoda sampler. malacoda tries to preserve the ordering
#'   present in the original mpra_data dataset (not including removed poorly
#'   represented barcodes), but if you're unsure you can pull from the
#'   variant_data column in the analysis_res object. Currently it's necessary
#'   that the allele be encoded by 'alt' and 'ref'.
#'
#'   Setting composite = TRUE will include three additional components: a column
#'   of tiles whose color indicates the log(mean(RNA)/ mean(DNA)) by row, and
#'   posterior interval plots for the variant's activities and transcription
#'   shift.
#' @export
mpra_tile_plot = function(variant_id, variant_counts,
                          show_barcodes = FALSE,
                          composite = FALSE,
                          sampler_res = NULL,
                          prob = .5,
                          prob_outer = .95){

  # TODO make these robust against other encodings of allele
  nref = sum(variant_counts$allele == 'ref')
  nalt = sum(variant_counts$allele != 'ref')

  dna_rep = variant_counts %>%
    select(.data$barcode, matches('DNA')) %>%
    gather("sample_id", "count", -.data$barcode) %>%
    group_by(.data$barcode) %>%
    summarise(dna_rep = mean(.data$count))

  var_tile_input = variant_counts %>%
    left_join(dna_rep, by = 'barcode') %>%
    gather("sample_id", "count", matches('[DR]NA', ignore.case = FALSE)) %>%
    arrange(.data$allele, .data$dna_rep) %>%
    mutate(barcode = factor(.data$barcode,
                            levels = unique(.data$barcode)))

  var_color_breaks = seq(min(var_tile_input$count),
                         max(var_tile_input$count),
                         length.out = 3)

  n_dna = names(variant_counts) %>% stringr::str_count('DNA') %>% sum
  n_rna = names(variant_counts) %>% stringr::str_count('RNA') %>% sum


  leg_pos = ifelse(composite, 'bottom', 'right' )
  axis_text_y = if (show_barcodes){
    element_text(hjust = 0)
  } else {
    element_blank()
  }

  var_tile_plot = var_tile_input %>%
    ggplot(aes(.data$sample_id, .data$barcode)) +
    geom_tile(aes(fill = .data$count)) +
    geom_segment(aes(x = .5,
                     xend = n_dna + n_rna + .5,
                     y = nalt+.5,
                     yend = nalt+.5),
                 lwd = 1.3,
                 color = 'grey50') +
    geom_vline(xintercept = n_dna + .5, color = 'white',
               lwd = 2) +
    scale_fill_viridis_c(na.value = 'black',
                         trans = 'log10',
                         guide = guide_colorbar(title.position = 'top')) +
    theme(axis.text.y = axis_text_y,
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = .5),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          legend.position = leg_pos) +
    coord_cartesian(expand = FALSE) +
    labs(title = variant_id,
         x = 'sample ID',
         fill = 'barcode count',
         y = 'alt   |   ref')

  if (composite) {
    if (missing(sampler_res)){
      stop('Composite plots require a sampler_res input.')
    }

    plot_probs = c((1 - prob_outer[1])/2,
                   (1 - prob) / 2,
                   .5,
                   (1 - prob) / 2 + prob,
                   (1 - prob_outer[1])/2 + prob_outer)

    fit_summary = sampler_res %>%
      as.data.frame() %>%
      as_tibble %>%
      select(.data$ref_act:.data$transcription_shift) %>%
      map_dfc(~quantile(.x, probs = plot_probs)) %>%
      mutate(quantile = plot_probs) %>%
      gather("param", "value", 1:3) %>%
      mutate(is_act = !grepl('act', .data$param)) %>%
      spread("quantile", "value") %>%
      mutate(param = factor(c('alt activity', 'ref activity', 'transcription\nshift')))

    act_ints = fit_summary %>% .[1:2,] %>%
      mutate(param = c('alt', 'ref')) %>%
      ggplot(aes(.data$`0.5`, .data$param)) +
      geom_segment(aes(x = .data$`0.025`, xend = .data$`0.975`,
                       yend = .data$param,
                       y = .data$param)) +
      geom_segment(aes(x = .data$`0.25`, xend = .data$`0.75`,
                       yend = .data$param,
                       y = .data$param),
                   color = 'deepskyblue2',
                   lwd = 2) +
      geom_point(size = 2) +
      theme_light() +
      theme(strip.text = element_blank(),
            axis.text.y = element_text(size = 16),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            plot.margin = margin(c(35, 5.5, 5.5, 5.5), unit = 'pt')) +
      labs(x = 'activity',
           y = NULL)

    ts_int = fit_summary %>%
      .[3,] %>%
      ggplot(aes(.data$`0.5`, 1)) +
      geom_segment(aes(x = .data$`0.025`, xend = .data$`0.975`,
                       yend = 1, y = 1)) +
      geom_segment(aes(x = .data$`0.25`, xend = .data$`0.75`,
                       yend = 1, y = 1),
                   color = 'deepskyblue2',
                   lwd = 2) +
      geom_point() +
      geom_vline(xintercept = 0,
                 lty = 2,
                 color = 'grey50') +
      labs(y = '   ',
           x = 'transcription shift') +
      theme_light() +
      theme(axis.ticks.y = element_blank(),
            axis.title.y = element_text(color = 'white'),
            axis.text.y = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            plot.margin = margin(c(25, 5.5, 35, 5.5), unit = 'pt'))

    int_grob = gridExtra::arrangeGrob(grobs = list(act_ints, ts_int), ncol = 1, heights = c(1.5,1))

    v_rna = variant_counts %>%
      select(matches('RNA')) %>%
      rowMeans()
    v_dna = variant_counts %>%
      select(matches('DNA')) %>%
      rowMeans()

    summary_input = variant_counts %>%
      mutate(avg_rna = v_rna,
             avg_dna = v_dna) %>%
      mutate(new_color = log10(.data$avg_rna / .data$avg_dna)) %>%
      mutate(barcode = factor(.data$barcode,
                              levels = levels(var_tile_input$barcode)))

    newcols = summary_input$new_color[is.finite(summary_input$new_color)]
    color_breaks = seq(min(newcols), max(newcols), length.out = 7)[c(2,4,6)]
    summary_color = summary_input %>%
      ggplot(aes('RNA1', .data$barcode)) +
      geom_tile(aes(color = .data$new_color, fill = .data$new_color)) +
      geom_segment(aes(x = .5,
                       xend = 1.5,
                       y = nalt+.5,
                       yend = nalt+.5),
                   lwd = 1.3,
                   color = 'grey40') +
      scale_fill_viridis_c(option = 'plasma',
                           guide = guide_colorbar(title.position = 'top'),
                           breaks = color_breaks,
                           labels = round(color_breaks, digits = 2),
                           na.value = 'black') +
      scale_color_viridis_c(option = 'plasma',
                            guide = guide_colorbar(title.position = 'top'),
                            breaks = color_breaks,
                            labels = round(color_breaks, digits = 2),
                            na.value = 'black') +
      theme_light() +
      coord_cartesian(expand = FALSE,
                      xlim = c(0,2)) +
      labs(title = '',
           y = NULL,
           x = 'sample ID',
           color = 'log(RNA/DNA)',
           fill = 'log(RNA/DNA)') +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.x = element_text(angle = 90, color = 'white'),
            axis.title.x = element_text(color = 'white'),
            legend.position = 'bottom',
            panel.grid = element_blank(),
            panel.border = element_blank(),
            axis.ticks.x = element_line(color = 'white'),
            plot.margin = margin(c(5.5, 11, 5.5, 0), unit = 'pt'),
            legend.box.margin = margin(c(0,0,0,0), unit = 'cm'))

    result_plot = gridExtra::grid.arrange(grobs = list(var_tile_plot, summary_color, int_grob), nrow = 1, widths = c(9,1.5, 3),
                                             newpage = FALSE)
  } else {
    result_plot = var_tile_plot
  }

  return(result_plot)
}

