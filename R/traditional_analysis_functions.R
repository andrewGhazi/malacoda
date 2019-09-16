#### Functions for things like computing activity levels, performing t-tests, etc.

#' Compute MPRA activities
#'
#' @description This function computes traditional MPRA activities, i.e. the
#'   log-ratio of RNA to DNA, for each barcode-RNA sample
#'
#' @param mpra_data a data frame of MPRA data
#' @param rep_cutoff a DNA representation cutoff
#' @param plot_rep_cutoff a logical indicating whether or not to plot the DNA
#'   barcode abundance and the applied representation cutoff
#' @param verbose logical indicating whether to print messages
#' @return a data frame of activities
#' @note the output is returned in a "tall" format, with sample_id's gathered
#'   into one column
#' @examples compute_activities(umpra_example)
#' @export
compute_activities = function(mpra_data,
                              rep_cutoff = .15,
                              plot_rep_cutoff = TRUE,
                              verbose = TRUE){

  sample_depths = get_sample_depths(mpra_data)

  if (verbose) {
    message('Determining well-represented variants, see plot...')
  }

  well_represented = get_well_represented(mpra_data,
                                          sample_depths,
                                          rep_cutoff = rep_cutoff,
                                          plot_rep_cutoff = plot_rep_cutoff,
                                          verbose = verbose)

  mean_dna_abundance = mpra_data %>%
    filter(.data$barcode %in% well_represented$barcode) %>%
    select(.data$barcode, matches('DNA')) %>%
    gather("sample_id", "count", -.data$barcode) %>%
    left_join(sample_depths, by = 'sample_id') %>%
    mutate(depth_adj_count = .data$count / .data$depth_factor) %>%
    group_by(.data$barcode) %>%
    summarise(mean_depth_adj_dna = mean(.data$depth_adj_count))

  activities_raw = mpra_data %>%
    select(-matches('DNA')) %>%
    gather("sample_id", "count", matches('RNA')) %>%
    left_join(sample_depths, by = 'sample_id') %>%
    left_join(mean_dna_abundance, by = 'barcode') %>%
    mutate(depth_adj_count = .data$count / .data$depth_factor,
           activity = log(.data$depth_adj_count / .data$mean_depth_adj_dna))

  if (any(is.infinite(activities_raw$activity))) {
    if (verbose) {
      message(paste0('Removing ', sum(is.infinite(activities_raw$activity)),
                     ' infinite activity measurements (i.e. 0 RNA counts) out of ',
                     nrow(activities_raw),
                     ' (', round(100*sum(is.infinite(activities_raw$activity)) / nrow(activities_raw), digits = 2), '%)'))
    }


    activities = activities_raw %>%
      filter(is.finite(.data$activity))
  }

  activities %<>%
    mutate(allele = factor(case_when(tolower(.data$allele) == 'ref' ~ 'ref',
                                     tolower(.data$allele) != 'ref' ~ 'alt'), levels = c('ref', 'alt')))


  return(activities)

}


test_one_variant = function(variant_activities,
                            test_type = 't') {

  if (n_distinct(variant_activities$allele) != 2 | any(table(variant_activities$allele) < 2)) {
    return(NA)
  }

  if (test_type %in% c('t', 't-test', 't.test')) {
    variant_activities %>%
      summarise(test_res = list(broom::tidy(t.test(x = .data$activity[tolower(.data$allele) == 'ref'],
                                                   y = .data$activity[tolower(.data$allele) != 'ref'])))) %>%
      unnest(c(.data$test_res)) %>%
      mutate(estimate = -.data$estimate) %>%
      dplyr::rename(ts_estimate = .data$estimate,
                    ref_mean_estimate = .data$estimate1,
                    alt_mean_estimate = .data$estimate2)
  } else if (test_type %in% c('u', 'U', 'wilcox.test', 'Mann.Whitney', 'Mann-Whitney')) {
    variant_activities %>%
      summarise(test_res = list(broom::tidy(wilcox.test(x = .data$activity[tolower(.data$allele) == 'ref'],
                                                        y = .data$activity[tolower(.data$allele) != 'ref'])))) %>%
      unnest(c(.data$test_res))
  } else {
    stop('test_type input given not supported')
  }
}

#' Get Sample Depths
#'
#' @description Computes the sum of all barcode counts
#' @param mpra_data a dataframe of MPRA data
#' @param depth_multiplier a numeric to divide through the depths to make the
#'   numbers more easily interpretable
#' @note The \code{depth_multiplier} input has no effect on the downstream
#'   analysis as it divides out because both RNA and DNA counts are normalized
#'   by this same factor. It simply sets the scale of depth factors at an
#'   easily-readable range.
#' @examples
#' get_sample_depths(umpra_example)
#' # Data from non-subsampled datasets with typical sequencing depths will typically
#' # show depth factors larger than this example, in the 10-100 range.
#' @export
get_sample_depths = function(mpra_data,
                             depth_multiplier = 1e6){
  mpra_data %>%
    gather('sample_id', 'counts', matches('DNA|RNA')) %>%
    group_by(.data$sample_id) %>%
    summarise(depth_factor = sum(.data$counts) / depth_multiplier)
}

#' Run MPRA activity tests
#'
#' @description Given a set of MPRA activity measurements, run significance
#'   tests to see if there's a difference in activity tests between alleles
#' @param mpra_activities a "tall" data frame of MPRA activities
#' @param test_type a string indicating the type of test to use to compare
#'   activities, either 't' or 'U'
#' @param compute_FDR logical indicating whether or not to add on a column of
#'   FDR q-values
#' @param plot_p_hist logical indicating whether or not to plot a histogram of p-values
#'
#' @details Compute activity tests by MPRA variant, comparing the activity
#'   measurements between alleles across all sequencing samples.
#' @return a data frame of test statistics
#' @note These tests assume that after depth correction, MPRA activities are
#'   identically distributed across sequencing samples. Thus activities for the
#'   same barcodes across different samples are considered as independent
#'   measurements.
#' @examples run_activity_tests(activities_example)
#' @export
run_activity_tests = function(mpra_activities,
                              test_type = 't',
                              compute_FDR = TRUE,
                              plot_p_hist = TRUE){

  activity_tests = mpra_activities %>%
    group_by(.data$variant_id) %>%
    nest() %>%
    dplyr::rename('act_dat' = 'data') %>%
    mutate(test_result = map(.data$act_dat, test_one_variant, test_type = test_type),
           note = map_chr(.data$test_result, ~ifelse(all(is.na(.x)),
                                                     'Not enough activity measurements present for both alleles, returning NA',
                                                     '')))

  na_variants = activity_tests %>% # lol I have no idea how to do this in a tidy way
    filter(.data$note != '') %>%
    select(.data$variant_id, .data$note)

  activity_tests = activity_tests %>%
    filter(map_lgl(.data$test_result, ~!all(is.na(.x)))) %>%
    select(-.data$act_dat) %>%
    unnest(c(.data$test_result)) %>%
    full_join(na_variants,
              by = c('variant_id', 'note')) %>%
    select(.data$variant_id, .data$ts_estimate:.data$alternative, .data$note) %>%
    ungroup

  if(compute_FDR){
    activity_tests = activity_tests %>%
      mutate(q_value = p.adjust(.data$p.value, method = 'fdr'))
  }

  if (plot_p_hist) {
    p_hist = activity_tests %>%
      ggplot(aes(.data$p.value)) +
      geom_histogram(boundary = 0,
                     binwidth = 1/40,
                     bins = 40,
                     color = 'black',
                     fill = 'grey50') +
      labs(title = 'Activity tests p-value distribution',
           x = 'p-value')
    print(p_hist)
  }

  return(activity_tests)

}
