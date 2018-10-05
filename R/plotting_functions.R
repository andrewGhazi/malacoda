#' Create a posterior beeswarm plot
#'
#' @description Create a combined beeswarm & activity posterior plot
#'
#' @param sampler_result a stanfit object
#' @param activities a dataframe of activities
#' @param color_by_sample optional logical indicating whether to color dots by sample
#' @return a posterior beeswarm ggplot object
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
