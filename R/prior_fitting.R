#' Generate a distance matrix from a matrix of annotations
#'
#' Given an nxd matrix of variant annotations, produce an nxn distance matrix
#' describing the inter-variant distances in annotation space
#'
#' @param annotations an n x d data frame of annotations
#' @param log_distance a logical indicating to use the log1p of the distances (TRUE) or the raw euclidean distances (FALSE)
#' @param scale_annotations logical indicating whether to base::scale to center and scale annotations
#' @param verbose logical indicating whether to print messages
#' @importFrom magrittr %>%
generate_distance_matrix = function(annotations,
                                    log_distance = FALSE,
                                    scale_annotations = TRUE,
                                    verbose = TRUE){

  if(nrow(annotations) > 10000 & verbose){
    message('Computing distance matrix for more than 10000 variants, I hope you have enough memory!')
  }

  if (log_distance & scale_annotations) {
    annotations %>%
      dplyr::mutate_at(.vars = vars(-.data$variant_id),
                       .funs = scale) %>%
      as.data.frame() %>%
      tibble::column_to_rownames('variant_id') %>%
      as.matrix %>%
      scale %>%
      dist %>%
      as.matrix %>%
      log1p
  } else if (log_distance & !scale_annotations) {
    annotations %>%
      as.data.frame() %>%
      tibble::column_to_rownames('variant_id') %>%
      as.matrix %>%
      scale %>%
      dist %>%
      as.matrix %>%
      log1p
  } else if (!log_distance & scale_annotations){
    annotations %>%
      dplyr::mutate_at(.vars = vars(-.data$variant_id),
                       .funs = scale) %>%
      as.data.frame() %>%
      tibble::column_to_rownames('variant_id') %>%
      as.matrix %>%
      scale %>%
      dist %>%
      as.matrix
  } else {
    annotations %>%
      as.data.frame() %>%
      tibble::column_to_rownames('variant_id') %>%
      as.matrix %>%
      scale %>%
      dist %>%
      as.matrix
  }
}

#' Find prior weights
#'
#' @description For a given variant and annotation set (scaled) and distance
#'   matrix, get the weights of all other variants
#' @param given_id a variant_id to get weights for
#' @param scaled_annotations a tall annotation data frame where the annotations
#'   have been set to the same scale
#' @param dist_mat a distance matrix of euclidean distances between variants in
#'   scaled annotation space
#' @param min_dist_kernel an initialization of the distance kernel at some tiny
#'   value
#' @param kernel_fold_change the multiplier by which to iteratively increase the
#'   distance kernel
#' @param min_num_neighbors the minimum number of neighbors which must make a
#'   meaningful contribution to the weights before stopping
#' @param verbose logical indicating whether to print messages
#' @details The "meaningful contribution" is defined in this way: The variants
#'   are sorted by weight. The min_num_neighbors-th (e.g. 100th) variant will be
#'   weighted to at least 1% of the highest most strongly weighted variant. This
#'   prevents some small number of extremely close neighbors from dominating the
#'   prior estimation later on.
find_prior_weights = function(given_id,
                              scaled_annotations,
                              dist_mat,
                              min_dist_kernel,
                              kernel_fold_change = 1.3,
                              min_num_neighbors = 100,
                              verbose = TRUE){

  n_annotations = dplyr::n_distinct(scaled_annotations$annotation)

  given_annotations = scaled_annotations %>%
    filter(.data$variant_id == given_id)

  pos_vec = given_annotations$value

  same_annotation_pos = scaled_annotations %>%
    filter(.data$variant_id != given_id) %>%
    group_by(.data$variant_id) %>%
    summarise(same_pos = all(.data$value == pos_vec)) %>%
    filter(.data$same_pos)

  if(nrow(same_annotation_pos) >= min_num_neighbors){
    if (verbose) {
      message('>= min_num_neighbors at the exact same annotation position. Using these evenly for prior estimation while not using others.')
    }

    weight_res = scaled_annotations %>%
      filter(.data$variant_id != given_id) %>%
      mutate(same_pos = .data$variant_id %in% same_annotation_pos$variant_id,
             weight = case_when(.data$same_pos ~ 1,
                                !.data$same_pos ~ 0)) %>%
      select(.data$variant_id, .data$weight)
    return(weight_res)

  }

  dist_to_others = scaled_annotations %>%
    filter(.data$variant_id != given_id) %>%
    group_by(.data$variant_id) %>%
    mutate(dist = .data$value - pos_vec)

  if (n_annotations == 1){
    weight_df = dist_to_others %>%
      ungroup %>%
      select(-.data$value) %>%
      mutate(mv_dens = mvtnorm::dmvt(matrix(.data$dist),
                                     sigma = diag(min_dist_kernel, n_annotations), log = FALSE)) %>% # Using a t kernel
      mutate(frac_weight = .data$mv_dens / sum(.data$mv_dens)) %>%
      arrange(desc(.data$frac_weight)) %>%
      mutate(cs = cumsum(.data$frac_weight),
             n = 1:dplyr::n())
  } else if (n_annotations > 1){
    dist_by_variant = dist_to_others %>%
      ungroup %>%
      select(-.data$value) %>%
      spread(.data$annotation, .data$dist) %>%
      select(-.data$variant_id) %>%
      as.matrix()

    variant_df = dist_to_others %>%
      ungroup %>%
      select(.data$variant_id) %>%
      unique

    weight_df = variant_df %>%
      mutate(mv_dens = mvtnorm::dmvt(dist_by_variant,
                                     sigma = diag(min_dist_kernel, n_annotations), log = FALSE)) %>% # Using a t kernel
      mutate(frac_weight = .data$mv_dens / sum(.data$mv_dens)) %>%
      arrange(desc(.data$frac_weight)) %>%
      mutate(cs = cumsum(.data$frac_weight),
             n = 1:dplyr::n())
  }

  if (weight_df$frac_weight[min_num_neighbors] / weight_df$frac_weight[1] < .01){
    # If the 100th variant has less than 1% the weight of the top variant,
    # expand the kernel width and re-weight.

    while (weight_df$frac_weight[min_num_neighbors] / weight_df$frac_weight[1] < .01) {
      min_dist_kernel = kernel_fold_change * min_dist_kernel

      if (n_annotations == 1){
        weight_df = dist_to_others %>%
          ungroup %>%
          select(-.data$value) %>%
          mutate(mv_dens = mvtnorm::dmvt(matrix(.data$dist),
                                         sigma = diag(min_dist_kernel, n_annotations), log = FALSE)) %>% # Using a t kernel
          mutate(frac_weight = .data$mv_dens / sum(.data$mv_dens)) %>%
          arrange(desc(.data$frac_weight)) %>%
          mutate(cs = cumsum(.data$frac_weight),
                 n = 1:dplyr::n())
      } else if (n_annotations > 1){

        weight_df = variant_df %>%
          mutate(mv_dens = mvtnorm::dmvt(dist_by_variant,
                                         sigma = diag(min_dist_kernel, n_annotations), log = FALSE)) %>% # Using a t kernel
          mutate(frac_weight = .data$mv_dens / sum(.data$mv_dens)) %>%
          arrange(desc(.data$frac_weight)) %>%
          mutate(cs = cumsum(.data$frac_weight),
                 n = 1:dplyr::n())
      }
    }
  }

  weight_res = weight_df %>%
    select(.data$variant_id, weight = .data$mv_dens)

  return(weight_res)

}

#' Get well represented barcodes
#'
#' @description Identify barcodes well-represented in DNA samples from input
#'   MPRA data.
#' @param mpra_data a data frame of MPRA data
#' @param sample_depths a data frame of sample depths
#' @param rep_cutoff a representation cutoff
#' @param plot_rep_cutoff logical indicating whether to plot the DNA
#'   representation distribution with the input rep_cutoff indicated
#' @param verbose logical indicating to print messages about DNA removal
#'   statistics
#' @details Use this function to tune the representation cutoff shown on the
#'   resulting histogram in order to discard failed and poorly-represented
#'   barcodes. These will sometimes be visible as a noticeable bump on the left
#'   side of the plot, though carefully prepared oligo libraries my not exhibit
#'   this.
#'
#'   After turning, plot_rep_cutoff may be set to FALSE to suppress the
#'   plotting.
#' @note properly formatted sample_depths can be obtained from
#'   \code{get_sample_depths()}
#' @examples
#' example_depths = get_sample_depths(umpra_example)
#' get_well_represented(umpra_example,
#'     sample_depths = example_depths,
#'     rep_cutoff = .15,
#'     plot_rep_cutoff = TRUE,
#'     verbose = TRUE)
#' # The final depth adjusted cutoff value will be lower for non-subsampled
#' # datasets that have higher total sequencing depth.
#' @export
get_well_represented = function(mpra_data,
                                sample_depths,
                                rep_cutoff,
                                plot_rep_cutoff = FALSE,
                                verbose = TRUE){

  if (verbose & plot_rep_cutoff) {
    message('Determining well-represented variants, see plot...')
  } else if (verbose) {
    message('Determining well-represented variants ...')
  }

  all_dna = mpra_data %>%
    select(.data$variant_id, .data$allele, .data$barcode, matches('DNA')) %>%
    gather('sample_id', 'counts', matches('DNA|RNA')) %>%
    left_join(sample_depths,
              by = 'sample_id') %>%
    mutate(depth_adj_count = .data$counts / .data$depth_factor) %>%
    group_by(.data$barcode) %>%
    summarise(mean_depth_adj_count = mean(.data$depth_adj_count))

  zero_num = all_dna %>%
    filter(.data$mean_depth_adj_count == 0) %>%
    nrow

  zero_pct = round(zero_num / nrow(all_dna) * 100, digits = 3)

  if (verbose){
    message(paste0('* ', zero_pct, '% of all barcodes have exactly 0 representation in all samples.'))
  }

  if (plot_rep_cutoff) {
     rep_cutoff_plot = all_dna %>%
       plot_dna_representation(rep_cutoff = rep_cutoff)

     print(rep_cutoff_plot)
  }

  well_represented = all_dna %>%
    filter(.data$mean_depth_adj_count > quantile(all_dna$mean_depth_adj_count,
                                                 probs = rep_cutoff)) %>%
    select(.data$barcode) %>%
    unique

  if (verbose) {
    message(paste0('* ', nrow(well_represented) , ' out of ', n_distinct(mpra_data$barcode),
                 ' (', round(100* nrow(well_represented) / n_distinct(mpra_data$barcode), digits = 2),'%)',
                 ' barcodes in input are well represented in the DNA pools.'))
  }

  return(well_represented)

}

#' Fit a marginal prior
#'
#' @param mpra_data a data frame of mpra data
#' @param n_cores number of cores to parallelize across
#' @param plot_rep_cutoff logical indicating whether to plot the representation
#'   cutoff used
#' @param rep_cutoff fraction indicating the depth-adjusted DNA count quantile
#'   to use as the cutoff
#' @param sample_depths optional inputs to allow passing in sample_depths and
#'   well_represented objects
#' @param well_represented optional inputs to allow passing in sample_depths and
#'   well_represented objects
#' @param verbose logical indicating whether to print messages
#' @examples marg_prior = fit_marg_prior(umpra_example,
#'  rep_cutoff = .15,
#'  plot_rep_cutoff = TRUE,
#'  n_cores = 1)
#' @export
fit_marg_prior = function(mpra_data,
                          n_cores = 1,
                          plot_rep_cutoff = TRUE,
                          rep_cutoff = .15,
                          sample_depths,
                          well_represented,
                          verbose = TRUE){

  if (missing(sample_depths) | missing(well_represented)){
    sample_depths = get_sample_depths(mpra_data)

    well_represented = get_well_represented(mpra_data,
                                            sample_depths,
                                            rep_cutoff = rep_cutoff,
                                            plot_rep_cutoff = plot_rep_cutoff,
                                            verbose = verbose)
  }

  if (!any(tolower(mpra_data$allele) == 'ref')) {
    stop('Cannot identify which alleles are reference or alternate. Map existing values onto "ref" or "alt" and try again. The function stringr::str_replace_all might help with this. Try str_replace_all(allele, c("A" = "ref", "B" = "alt")) where allele is the existing allele column and "A" and "B" are the current allele indicators.')
  }

  if (verbose) {
    message('Fitting negative binomial MLEs...')
  }



  # All Negative Binomial Maximum Likelihood Estimates
  all_nb_mle = fit_all_nb_mle(mpra_data,
                              well_represented = well_represented,
                              sample_depths = sample_depths,
                              n_cores = n_cores)


  # MLEs of dispersion parameters are known to be biased to be too large. The
  # exact amount of bias depends on number of data points and the true, values,
  # but as a heuristic we simply remove the top 5% of largest dispersion values.
  estimates_for_priors = all_nb_mle %>%
    select(matches('_p|_m')) %>%
    gather('par', 'value') %>%
    group_by(.data$par) %>%
    mutate(to_keep = grepl('_m', .data$par) | .data$value < quantile(.data$value, .95)) %>%
    ungroup %>%
    filter(.data$to_keep) %>%
    select(-.data$to_keep)

  gamma_estimates = estimates_for_priors %>%
    group_by(.data$par) %>%
    summarise(prior = list(fit_gamma_stan(.data$value))) %>%
    mutate(alpha_est = map_dbl(.data$prior, ~.x$par[['alpha']]),
           beta_est = map_dbl(.data$prior, ~.x$par[['beta']]),
           acid_type = toupper(str_extract(string = .data$par, pattern = '[dr]na')),
           allele = c('[1]' = 'ref', '[2]' = 'alt')[str_extract(pattern = '\\[[12]\\]',
                                                                string = .data$par)]) %>%  # kill me
    dplyr::rename('prior_type' = 'par') %>%
    mutate(prior_type = str_replace_all(str_extract(.data$prior_type,
                                                    'm|p'),
                                        c('p' = 'phi_prior',
                                          'm' = 'mu_prior')))

  # There is room for improvement with this phi prior estimation. It disregards
  # sequencing sample depth V and cuts out observed larged phi values.

  return(gamma_estimates)


  # mpra_data %>%
  #   select(variant_id, allele, barcode, matches('RNA')) %>%
  #   gather(sample_id, counts, matches('DNA|RNA')) %>%
  #   left_join(sample_depths , by = 'sample_id') %>%
  #   left_join(mean_dna_abundance, by = 'barcode') %>%
  #   mutate(count_remnant = .1 + counts / depth_factor / mean_depth_adj_count) %>%
  #   ggplot(aes(count_remnant)) +
  #   geom_density(aes(color = sample_id)) +
  #   scale_x_log10() +
  #   stat_function(fun = dgamma, args = list(shape = 1.04, rate = 1.128))

  # ^ per-sample priors might be worthwhile

}

#' Fit priors by group
#'
#' @description If you have several distinct categories of variants, one may
#'   want to fit priors for them separately. Categories could be genomic region:
#'   5'UTR vs intronic vs 3'UTR vs upstream vs downstream. Perhaps you want to
#'   quantify the difference by some prediction outputs: up vs down vs
#'   no-effect.
#'
#'   This yields a pseudo-hierarchical model without the computational problems
#'   associated with fitting a joint model on thousands of variants at once.
#' @inheritParams fit_marg_prior
#' @param group_df a data frame giving group identity by variant_id in mpra_data
#'
#' @details group_df should have two columns: variant_id and group_id. This
#'   function checks that there are >100 variants per group and that there
#'   aren't more than 20 groups. These are somewhat arbitrary magic numbers, but
#'   having loads of tiny groups is a recipe for over-fitting.
#'
#' @return a grouped prior list
fit_grouped_prior = function(mpra_data,
                             group_df,
                             n_cores,
                             plot_rep_cutoff = TRUE,
                             rep_cutoff = .15,
                             verbose = TRUE) {



  #### Input checks ----

  # Enough per group?
  group_sizes = group_df %>%
    dplyr::count(.data$group_id,
                 name = 'n')

  enough_per_group = all(group_sizes$n > 100)
  # 100 is just a magic number. Should I make that a user input? Or make it just
  # completely override-able?

  if (!enough_per_group) {
    stop('Not enough variants in all the groups!')
  }

  if (!any(tolower(mpra_data$allele) == 'ref')) {
    stop('Cannot identify which alleles are reference or alternate. Map existing values onto "ref" or "alt" and try again. The function stringr::str_replace_all might help with this. Try str_replace_all(allele, c("A" = "ref", "B" = "alt")) where allele is the existing allele column and "A" and "B" are the current allele indicators.')
  }

  # Not too many groups?
  num_groups = group_df$group_id %>%
    dplyr::n_distinct()

  too_many_groups = num_groups > 20

  if (too_many_groups) {
    stop('Too many groups!')
  }

  #### Preparation ----

  # we want to avoid computing sample depths and well represented barcodes
  # separately for each group, so we do it here globally first. Then pass these
  # results to fit_marg_prior.

  sample_depths = get_sample_depths(mpra_data)

  well_represented = get_well_represented(mpra_data,
                                          sample_depths,
                                          rep_cutoff = rep_cutoff,
                                          plot_rep_cutoff = plot_rep_cutoff,
                                          verbose = verbose)

  #### Fit grouped prior ----
  grouped_prior = mpra_data %>%
    left_join(group_df, by = 'variant_id') %>%
    dplyr::group_by(.data$group_id) %>%
    nest() %>%
    dplyr::rename('group_data' = 'data') %>%
    ungroup %>%
    dplyr::mutate(group_prior = purrr::map(group_data, fit_marg_prior,
                                           sample_depths = sample_depths,
                                           well_represented = well_represented,
                                           verbose = verbose,
                                           n_cores = n_cores, plot_rep_cutoff = FALSE, rep_cutoff = .15)) %>%
    dplyr::select(-'group_data')

  return(grouped_prior)
}

#' Fit a informative conditional prior
#'
#' @description Use informative annotations to bias prior estimation towards
#'   alleles that show similar annotations in the provided annotation space.
#'
#' @details The empirical prior returned by this object is "conditional" in the
#'   sense that the prior estimation weights are conditional on the annotations.
#'
#'   The DNA prior is still estimated marginally because the annotations should
#'   not be able to provide any information on the DNA inputs (which are
#'   presumably only affected by the preparation of the oligonucleotide library
#'   at the vendor).
#'
#'   The RNA prior is estimated from the RNA observations of
#'   other variants in the assay that are nearby in annotation space. A
#'   multivariate t distribution centered on the variant in question is used to
#'   weight all other variants in the assay. It is initialized with a very small
#'   width, and if there are fewer than \code{min_neighbors} that provide
#'   substantial input to the prior, the width is iteratively increased by a
#'   factor of \code{kernel_fold_increase} until that condition is satisfied.
#'   This prevents the prior estimation for variants in sparse regions of
#'   annotation space from being influenced too heavily by their nearest
#'   neighbors.
#'
#' @return A list of two data frames. The first is for the DNA and the second is
#'   by-variant RNA priors.
#'
#' @param mpra_data a data frame of mpra data
#' @param annotations a data frame of annotations for the same variants in
#'   mpra_data
#' @param n_cores number of cores to parallelize across
#' @param plot_rep_cutoff logical indicating whether to plot the representation
#'   cutoff used
#' @param rep_cutoff fraction indicating the depth-adjusted DNA count quantile
#'   to use as the cutoff
#' @param min_neighbors The minimum number of neighbors in annotation space that
#'   must contribute to prior estimation
#' @param kernel_fold_increase The amount to iteratively increase kernel width
#'   by when estimating conditional priors. Smaller values (closer to 1) will
#'   yield more refined priors but take longer.
#' @param verbose logical indicating whether to print messages
#' @export
#' @examples
#' cond_prior = fit_cond_prior(mpra_data = umpra_example,
#'                             annotations = u_deepsea,
#'                             n_cores = 1,
#'                             rep_cutoff = .15,
#'                             plot_rep_cutoff = TRUE,
#'                             min_neighbors = 5)
fit_cond_prior = function(mpra_data,
                          annotations,
                          n_cores = 1,
                          plot_rep_cutoff = TRUE,
                          rep_cutoff = .15,
                          min_neighbors = 100,
                          kernel_fold_increase = 1.4142,
                          verbose = TRUE){

  # Input checks ----
  if(!all(mpra_data$variant_id %in% annotations$variant_id)){
    stop('There aren\'t annotations for each variant!')
  }

  if (!any(tolower(mpra_data$allele) == 'ref')) {
    stop('Cannot identify which alleles are reference or alternate. Map existing values onto "ref" or "alt" and try again. The function stringr::str_replace_all might help with this. Try str_replace_all(allele, c("A" = "ref", "B" = "alt")) where allele is the existing allele column and "A" and "B" are the current allele indicators.')
  }

  if(any(duplicated(annotations$variant_id))){
    stop('There are duplicate annotations for some variants!')
  }

  if (n_distinct(annotations$variant_id) > n_distinct(mpra_data$variant_id)){
    n_extra = n_distinct(annotations$variant_id) - n_distinct(mpra_data$variant_id)

    if (verbose) {
      message(paste0('Annotations provided for variants not included in mpra_data. Removing ', n_extra, ' unneeded annotations.'))
    }

    annotations = dplyr::filter(annotations,
                                .data$variant_id %in% mpra_data$variant_id)
  }

  if (nrow(unique(select(annotations, -.data$variant_id))) < 10) {
    warning('Found fewer than 10 unique annotations values. Were you looking for fit_grouped_prior()?')
  }

  #### DNA prior fitting ----

  if (verbose) {
    message('Evaluating data depth/DNA representation properties...')
  }

  sample_depths = get_sample_depths(mpra_data)

  well_represented = get_well_represented(mpra_data,
                                          sample_depths,
                                          rep_cutoff = rep_cutoff,
                                          plot_rep_cutoff = plot_rep_cutoff,
                                          verbose = verbose)

  if (verbose) {
    message('Fitting per-variant maximum likelihood estimates...')
  }

  all_nb_mle = fit_all_nb_mle(mpra_data,
                              well_represented = well_represented,
                              sample_depths = sample_depths,
                              n_cores = n_cores)

  if (!all(all_nb_mle$converged) & verbose) {
    message(paste0(sum(!all_nb_mle$converged),
                   ' out of ',
                   nrow(all_nb_mle),
                   ' (', round(sum(!all_nb_mle$converged) / nrow(all_nb_mle) * 100, digits = 3), '%)',
                   ' negative binomial maximum likelihood estimates failed to converge. A small fraction (<1%) failing is acceptable.'))
  }


  # MLEs of dispersion parameters are known to be biased to be too large. The
  # exact amount of bias depends on number of data points and the true, values,
  # but as a heuristic we simply remove the top 5% of largest dispersion values.
  estimates_for_marg = all_nb_mle %>%
    select(matches('_p|_m')) %>%
    gather('par', 'value') %>%
    group_by(.data$par) %>%
    mutate(to_keep = grepl('_m', .data$par) | .data$value < quantile(.data$value, .95)) %>%
    ungroup %>%
    filter(.data$to_keep) %>%
    select(-.data$to_keep) %>%
    filter(grepl('dna', .data$par))

  dna_gamma_prior = estimates_for_marg %>%
    group_by(.data$par) %>%
    summarise(prior = list(fit_gamma_stan(.data$value))) %>%
    mutate(alpha_est = map_dbl(.data$prior, ~.x$par[['alpha']]),
           beta_est = map_dbl(.data$prior, ~.x$par[['beta']]),
           acid_type = toupper(str_extract(string = .data$par, pattern = '[dr]na')),
           allele = c('[1]' = 'ref', '[2]' = 'alt')[str_extract(pattern = '\\[[12]\\]',
                                                                string = .data$par)]) %>%  # kill me
    dplyr::rename('prior_type' = 'par') %>%
    mutate(prior_type = str_replace_all(str_extract(.data$prior_type,
                                                    'm|p'),
                                        c('p' = 'phi_prior',
                                          'm' = 'mu_prior')))

  #### Fit conditional RNA prior ----

  # generate annotation distance matrix
  dist_mat = generate_distance_matrix(annotations = annotations,
                                      verbose = verbose)

  scaled_annotations = annotations %>%
    mutate_at(.vars = vars(-.data$variant_id),
              .funs = scale) %>%
    gather('annotation', 'value', -.data$variant_id) %>%
    arrange(.data$variant_id, .data$annotation)

  n_annotations = dplyr::n_distinct(scaled_annotations$annotation)

  min_dist_kernel = dist_mat[upper.tri(dist_mat)] %>%
    unlist() %>%
    sort() %>% #sort all observed distances
    .[. > 0] %>%
    quantile(probs = .001) %>%
    {. ^ (1 / n_annotations)} # adjustment for the number of annotations provided

  # For each variant, get a vector of weights for all other variants in the assay
  if (verbose) {
    message('Weighting variants in annotation space')
  }

  annotation_vectors = annotations %>%
    gather('anno_id', 'anno_value', -.data$variant_id) %>%
    arrange(.data$anno_id) %>%
    group_by(.data$variant_id) %>%
    summarise(anno_labels = list(c(.data$anno_id)),
              anno_vec = list(c(.data$anno_value)))

  unique_annotations = annotation_vectors %>% # You only need to fit the priors once for each unique vector of annotations
    filter(!duplicated(.data$anno_vec))

  prior_weights = unique_annotations %>%
    mutate(annotation_weights = parallel::mclapply(.data$variant_id, find_prior_weights,
                                                   scaled_annotations = scaled_annotations,
                                                   dist_mat = dist_mat,
                                                   min_dist_kernel = min_dist_kernel,
                                                   min_num_neighbors = min_neighbors,
                                                   kernel_fold_change = kernel_fold_increase,
                                                   verbose = verbose,
                                                   mc.cores = n_cores))


  # Perform weighted density estimation for each variant
  if (verbose) {
    message('Fitting annotation-weighted distributions...')
  }

  estimates_for_weighted = all_nb_mle %>%
    select(.data$variant_id, matches('rna_p|rna_m')) %>%
    gather('par', 'value', -.data$variant_id) %>%
    group_by(.data$par) %>%
    mutate(to_keep = grepl('_m', .data$par) | .data$value < quantile(.data$value, .95)) %>%
    ungroup %>%
    filter(.data$to_keep) %>%
    select(-.data$to_keep)

  if (n_cores == 1) {

    rna_priors = prior_weights %>%
      mutate(gamma_priors = map(.data$annotation_weights,
                                fit_weighted_gammas_stan,
                                estimates_to_weight = estimates_for_weighted)) %>%
      mutate(variant_m_prior = map(.data$gamma_priors, 1), # This part is necessary to make it work with the existing code
             variant_p_prior = map(.data$gamma_priors, 2)) %>%
      select(-.data$gamma_priors)

  } else if (n_cores > 1) {

    rna_priors = prior_weights %>%
      mutate(gamma_priors = parallel::mclapply(.data$annotation_weights,
                                               fit_weighted_gammas_stan,
                                               estimates_to_weight = estimates_for_weighted,
                                               mc.cores = n_cores,
                                               mc.preschedule = FALSE)) %>%
      mutate(variant_m_prior = map(.data$gamma_priors, 1), # This part is necessary to make it work with the existing code
             variant_p_prior = map(.data$gamma_priors, 2)) %>%
      select(-.data$gamma_priors)
  }

  #TODO tack the appropriate priors onto annotation_vectors according to the value of anno_vec
  # Normal join operations don't work on list columns, so you'll need to match them up manually

  annotation_vectors$weight_row = map_int(annotation_vectors$anno_vec, # kill me
                                          ~which(sapply(prior_weights$anno_vec, function(x){all(.x == x)}))) # later

  annotation_vectors$annotation_weights = prior_weights$annotation_weights[annotation_vectors$weight_row]
  annotation_vectors$variant_m_prior = rna_priors$variant_m_prior[annotation_vectors$weight_row]
  annotation_vectors$variant_p_prior = rna_priors$variant_p_prior[annotation_vectors$weight_row]

  both_rna = annotation_vectors %>%
    dplyr::select(.data$variant_id,
                  .data$annotation_weights,
                  .data$variant_m_prior,
                  .data$variant_p_prior)

  res_list = list(dna_prior = dna_gamma_prior,
                  rna_priors = both_rna)

  return(res_list)
}


fit_one_m_prior = function(given_id,
                           annotation_weights,
                           mpra_data,
                           sample_depths,
                           well_represented,
                           mean_dna_abundance){

mpra_data %>%
    filter(.data$barcode %in% well_represented$barcode) %>%
    filter(.data$variant_id != given_id) %>%
    select(-matches('DNA')) %>%
    left_join(annotation_weights, by = 'variant_id') %>%
    gather('sample_id', 'counts', matches('RNA')) %>%
    left_join(sample_depths, by = 'sample_id') %>%
    left_join(mean_dna_abundance, by = 'barcode') %>%
    mutate(count_remnant = .1 + .data$counts / .data$depth_factor / .data$mean_depth_adj_count) %>% # the variability of the count after accounting for depth and DNA input
    group_by(.data$allele) %>%
    summarise(mu_prior = list(fit_gamma(.data$count_remnant,
                                        weights = .data$weight))) %>% # WEIGHTS! A CONDITIONALLY WEIGHTED PRIOR!
    gather('prior_type', 'prior', matches('prior')) %>%
    mutate(alpha_est = map_dbl(.data$prior, ~.x$par[1]),
           beta_est = map_dbl(.data$prior, ~.x$par[2])) %>%
    mutate(acid_type = 'RNA')
}

fit_one_p_prior = function(given_id,
                           annotation_weights,
                           mpra_data,
                           sample_depths,
                           well_represented,
                           mean_dna_abundance){

  mpra_data %>%
    filter(.data$barcode %in% well_represented$barcode) %>%
    filter(.data$variant_id != given_id) %>%
    select(-matches('DNA')) %>%
    left_join(annotation_weights, by = 'variant_id') %>%
    gather('sample_id', 'counts', matches('RNA')) %>%
    group_by(.data$allele, .data$barcode) %>%
    summarise(mean_est = mean(.data$counts),
              var_est = var(.data$counts),
              size_guess = .data$mean_est^2 / (.data$var_est - .data$mean_est),
              weight = .data$weight[1]) %>%
    filter(.data$size_guess > 0 & is.finite(.data$size_guess)) %>% # negative size guess = var < mean = underdispersed
    filter(.data$size_guess < quantile(.data$size_guess, probs = .99)) %>%
    summarise(phi_prior = list(fit_gamma(.data$size_guess,
                                         weights = .data$weight))) %>% # WEIGHTS!
    gather('prior_type', 'prior', matches('prior')) %>%
    mutate(alpha_est = map_dbl(.data$prior, ~.x$par[1]),
           beta_est = map_dbl(.data$prior, ~.x$par[2])) %>%
    mutate(acid_type = 'RNA')
}

format_conditional_prior = function(given_id, cond_priors){
  dna_prior = cond_priors$dna_prior
  rna_prior = cond_priors$rna_priors %>%
    filter(.data$variant_id == given_id) %>%
    select(-.data$annotation_weights) %>%
    gather('prior_name', 'prior', matches('prior')) %>%
    unnest(c(.data$prior)) %>%
    select(-.data$variant_id, -.data$prior_name)

  return(bind_rows(dna_prior, rna_prior))

}

#' Increase a malacoda prior object's regularization
#'
#' @description This function takes a malacoda prior object and increases the
#'   regularization on transcription shift by a set factor.
#' @details This function multiplies the alpha and beta parameters of the gamma
#'   priors on the mean for the RNA counts in order to keep the mean the same
#'   while ascribing less prior density to extreme values, thus effectively
#'   increasing the shrinkage applied to transcription shift while evaluating
#'   the posterior.
#' @param prior_object a marginal prior object returned by fit_marg_prior()
#' @param increase_factor a number indicating the factor by which to increase
#' @return A malacoda prior object with increase
#' @examples increase_regularization(marg_prior_example)
#' @seealso \code{\link{fit_marg_prior}} \code{\link{sample_from_prior}}
#'   \code{\link{summarise_prior_samples}}
#' @export
increase_regularization = function(prior_object,
                                   increase_factor = 2){
  if (class(prior_object)[1] == 'list'){
    stop('This function currently only works for marginal priors')
  }

  row_id = prior_object$prior_type == 'mu_prior' & prior_object$acid_type == 'RNA'
  prior_object$alpha_est[row_id] = increase_factor * prior_object$alpha_est[row_id]
  prior_object$beta_est[row_id] = increase_factor * prior_object$beta_est[row_id]

  return(prior_object)
}
