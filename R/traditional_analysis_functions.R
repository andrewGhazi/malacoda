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
#'
#' @return a data frame of activities
#' @note the output is returned in a "tall" format, with sample_id's gathered
#'   into one column
#' @export
compute_activities = function(mpra_data,
                             rep_cutoff = .15,
                             plot_rep_cutoff = TRUE){

  sample_depths = mpra_data %>%
    gather(sample_id, counts, matches('DNA|RNA')) %>%
    group_by(sample_id) %>%
    summarise(depth_factor = sum(counts) / 1e6)

  print('Determining well-represented variants, see plot...')
  well_represented = get_well_represented(mpra_data,
                                          sample_depths,
                                          rep_cutoff = rep_cutoff,
                                          plot_rep_cutoff = plot_rep_cutoff,
                                          verbose = TRUE)

  mean_dna_abundance = mpra_data %>%
    filter(barcode %in% well_represented$barcode) %>%
    select(barcode, matches('DNA')) %>%
    gather(sample_id, count, -barcode) %>%
    left_join(sample_depths, by = 'sample_id') %>%
    mutate(depth_adj_count = count / depth_factor) %>%
    group_by(barcode) %>%
    summarise(mean_depth_adj_dna = mean(depth_adj_count))

  activities_raw = mpra_data %>%
    select(-matches('DNA')) %>%
    gather(sample_id, count, matches('RNA')) %>%
    left_join(sample_depths, by = 'sample_id') %>%
    left_join(mean_dna_abundance, by = 'barcode') %>%
    mutate(depth_adj_count = count / depth_factor,
           activity = log(depth_adj_count / mean_depth_adj_dna))

  if (any(is.infinite(activities_raw$activity))) {
    print(paste0('Removing ', sum(is.infinite(activities_raw$activity)),
                 ' infinite activity measurements (i.e. 0 RNA counts) out of ',
                 nrow(activities_raw),
                 ' (', round(100*sum(is.infinite(activities_raw$activity)) / nrow(activities_raw), digits = 2), '%)'))

    activities = activities_raw %>%
      filter(is.finite(activity))
  }

  activities %<>%
    mutate(allele = factor(case_when(tolower(allele) == 'ref' ~ 'ref',
                                     tolower(allele) != 'ref' ~ 'alt'), levels = c('ref', 'alt')))


  return(activities)

}


test_one_variant = function(variant_activities,
                            test_type = 't') {

  if (n_distinct(variant_activities$allele) != 2) {

    return(NA)
  }

  if (test_type %in% c('t', 't-test', 't.test')) {
    variant_activities %>%
      summarise(test_res = list(broom::tidy(t.test(x = activity[tolower(allele) == 'ref'],
                                                   y = activity[tolower(allele) != 'ref'])))) %>%
      unnest %>%
      rename(ts_estimate = estimate,
             ref_mean_estimate = estimate1,
             alt_mean_estimate = estimate2)
  } else if (test_type %in% c('u', 'U', 'wilcox.test', 'Mann.Whitney', 'Mann-Whitney')) {
    variant_activities %>%
      summarise(test_res = list(broom::tidy(wilcox.test(x = activity[tolower(allele) == 'ref'],
                                                        y = activity[tolower(allele) != 'ref'])))) %>%
      unnest
  } else {
    stop('test_type input given not supported')
  }
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
#' @export
run_activity_tests = function(mpra_activities,
                              test_type = 't',
                              compute_FDR = TRUE,
                              plot_p_hist = TRUE){

  activity_tests = mpra_activities %>%
    group_by(variant_id) %>%
    nest(.key = act_dat) %>%
    mutate(test_result = map(act_dat, test_one_variant, test_type = test_type),
           note = map_chr(test_result, ~ifelse(all(is.na(.x)),
                                               'Activity measurements not present for both alleles, returning NA',
                                               '')))

  na_variants = activity_tests %>% # lol I have no idea how to do this in a tidy way
    filter(note != '') %>%
    select(variant_id, note)

  activity_tests = activity_tests %>%
    filter(map_lgl(test_result, ~!all(is.na(.x)))) %>%
    select(-act_dat) %>%
    unnest() %>%
    full_join(na_variants,
              by = c('variant_id', 'note')) %>%
    select(variant_id, ts_estimate:alternative, note)

  if(compute_FDR){
    activity_tests = activity_tests %>%
      mutate(q_value = p.adjust(p.value, method = 'fdr'))
  }

  if (plot_p_hist) {
    p_hist = activity_tests %>%
      ggplot(aes(p.value)) +
      geom_histogram(boundary = 0,
                     bins = 30,
                     color = 'black',
                     fill = 'grey50') +
      labs(title = 'Activity tests p-value distribution')
    print(p_hist)
  }

  return(activity_tests)

}
