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
compute_activites = function(mpra_data,
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
                                               plot_rep_cutoff = plot_rep_cutoff)



  mean_dna_abundance = mpra_data %>%
    filter(barcode %in% well_represented$barcode) %>%
    select(barcode, matches('DNA')) %>%
    gather(sample_id, count, -barcode) %>%
    left_join(sample_depths, by = 'sample_id') %>%
    mutate(depth_adj_count = count / depth_factor) %>%
    group_by(barcode) %>%
    summarise(mean_depth_adj_dna = mean(depth_adj_count))

  activities = mpra_data %>%
    select(-matches('DNA')) %>%
    gather(sample_id, count, matches('RNA')) %>%
    left_join(sample_depths, by = 'sample_id') %>%
    left_join(mean_dna_abundance, by = 'barcode') %>%
    mutate(depth_adj_count = count / depth_factor,
           activity = log(depth_adj_count / mean_depth_adj_dna))


  return(activities)

}


test_one_variant = function(variant_activities,
                            test_type = 't'){

  variant_activities %>%
    summarise(test_res = list(t.test(x = activity[tolower(allele) == 'ref'],
                                     y = activity[tolower(allele) != 'ref'])))

}

run_activity_tests = function(mpra_activities,
                              test_type = 't'){

  activity_tests = mpra_activities %>%
    group_by(variant_id) %>%
    nest(.key = act_dat) %>%
    mutate(test_result = map(act_dat, test_one_variant))

}
