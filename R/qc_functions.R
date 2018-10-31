
#' Get MPRA sample correlations
#'
#' @description Compute the correlations of MPRA sequencing samples. This can be
#'   a helpful as a QC diagnostic metric. Libraries of the same type (DNA or RNA) should
#'   correlate highly with one another
#' @param mpra_data a data frame of mpra data
#' @return a data frame of pairwise correlations
#' @details the output is ready to be presented to plot_mpra_correlations
get_sample_correlations = function(mpra_data){
  mpra_data %>%
    select(matches('[DR]NA')) %>%
    cor() %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = 'sample_1') %>%
    gather("sample_2", 'correlation', -.data$sample_1)
}
