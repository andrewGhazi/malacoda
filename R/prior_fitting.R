#' Generate a distance matrix from a matrix of annotations
#'
#' Given an nxd matrix of variant annotations, produce an nxn distance matrix
#' describing the inter-variant distances in annotation space
#'
#' @param annotations an n x d data frame of annotations
#' @param log_distance a logical indicating to use the log1p of the distances (TRUE) or the raw euclidean distances (FALSE)
#' @param scale_annotations logical indicating whether to base::scale to center and scale annotations
#'
#'
#' @importFrom magrittr %>%
generate_distance_matrix = function(annotations,
                                    log_distance = FALSE,
                                    scale_annotations = TRUE){

  if(nrow(annotations) > 10000){
    warning('Computing distance matrix for more than 10000 variants, I hope you have enough memory!')
  }

  if (log_distance & scale_annotations) {
    annotations %>%
      dplyr::mutate_at(.vars = vars(-variant_id),
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
      dplyr::mutate_at(.vars = vars(-variant_id),
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
#' @details "Meaningful contribution" to the weights means that the sum of the
#'   weights of the first min_num_neighbors of other variants must account for
#'   <= 99% of the total weight. This prevents some small number of extremely
#'   close neighbors from dominating the prior estimation later on.
#' @export
find_prior_weights = function(given_id,
                              scaled_annotations,
                              dist_mat,
                              min_dist_kernel,
                              kernel_fold_change = 1.3,
                              min_num_neighbors = 30){

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
    print('>= min_num_neighbors at the exact same annotation position. Using these evenly for prior estimation while not using others.')

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

  weight_df = dist_to_others %>%
    select(-.data$value) %>%
    summarise(mv_dens = mvtnorm::dmvt(.data$dist,
                                      sigma = diag(min_dist_kernel, n_annotations), log = FALSE)) %>% # Using a t kernel
    mutate(frac_weight = .data$mv_dens / sum(.data$mv_dens)) %>%
    arrange(desc(.data$frac_weight)) %>%
    mutate(cs = cumsum(.data$frac_weight),
           n = 1:n())

  if (weight_df$cs[min_num_neighbors] > .99){
    # If the first 30 (min_num_neighbors) weights account for more than 99% of
    # all weight, we need to increase the kernel and try again

    while (weight_df$cs[min_num_neighbors] > .99) {
      min_dist_kernel = kernel_fold_change * min_dist_kernel
      weight_df = dist_to_others %>%
        select(-.data$value) %>%
        summarise(mv_dens = mvtnorm::dmvt(.data$dist, sigma = diag(min_dist_kernel, n_annotations), log = FALSE)) %>% # Using a t kernel
        mutate(frac_weight = .data$mv_dens / sum(.data$mv_dens)) %>%
        arrange(desc(.data$frac_weight)) %>%
        mutate(cs = cumsum(.data$frac_weight),
               n = 1:n())
    }
  }

  weight_res = weight_df %>%
    select(.data$variant_id, weight = .data$mv_dens)

  return(weight_res)

}

get_well_represented = function(mpra_data,
                                     sample_depths,
                                     rep_cutoff,
                                     plot_rep_cutoff = FALSE,
                                verbose = TRUE){

  all_dna = mpra_data %>%
    select(.data$variant_id, .data$allele, .data$barcode, matches('DNA')) %>%
    gather('sample_id', 'counts', matches('DNA|RNA')) %>%
    left_join(sample_depths,
              by = 'sample_id') %>%
    mutate(depth_adj_count = .data$counts / .data$depth_factor) %>%
    group_by(.data$barcode) %>%
    summarise(mean_depth_adj_count = mean(.data$depth_adj_count))


  if (plot_rep_cutoff) {
     rep_cutoff_plot = all_dna %>%
      ggplot(aes(.data$mean_depth_adj_count)) +
      geom_histogram(bins = 40,
                     color = 'black',
                     fill = 'grey50') +
      geom_vline(xintercept = quantile(all_dna$mean_depth_adj_count,
                                       probs = rep_cutoff),
                 lty = 2) +
      scale_x_log10() +
      labs(x = 'Mean Depth Adjusted DNA barcode count',
           title = 'DNA barcode abundance and cutoff',
           subtitle = paste0('Using depth-adjusted DNA barcode count cutoff of ',
                             round(quantile(all_dna$mean_depth_adj_count,
                                            probs = rep_cutoff),
                                   digits = 3))) +
       geom_rug(alpha = .01)

     print(rep_cutoff_plot)
  }

  well_represented = all_dna %>%
    filter(mean_depth_adj_count > quantile(all_dna$mean_depth_adj_count,
                                           probs = rep_cutoff)) %>%
    select(.data$barcode) %>%
    unique

  if (verbose) {
    print(paste0(nrow(well_represented) , ' out of ', n_distinct(mpra_data$barcode),
                 ' (', round(100* nrow(well_represented) / n_distinct(mpra_data$barcode), digits = 2),'%)',
                 ' barcodes in input are well represented in the DNA pools.'))
  }

  return(well_represented)

}

#' Fit a marginal prior
#'
#' @param mpra_data a data frame of mpra data
#' @param n_cores number of cores to parallelize across
#' @param plot_rep_cutoff logical indicating whether to plot the representation cutoff used
#' @param rep_cutoff fraction indicating the depth-adjusted DNA count quantile to use as the cutoff
#' @export
fit_marg_prior = function(mpra_data,
                          n_cores = 1,
                          plot_rep_cutoff = TRUE,
                          rep_cutoff = .15){

  sample_depths = get_sample_depths(mpra_data)

  print('Determining well-represented variants, see plot...')
  well_represented = get_well_represented(mpra_data,
                                          sample_depths,
                                          rep_cutoff = rep_cutoff,
                                          plot_rep_cutoff = plot_rep_cutoff)


  print('Fitting marginal DNA prior...')

  dna_nb_fits = mpra_data %>%
    filter(.data$barcode %in% well_represented$barcode) %>%
    select(.data$variant_id, .data$allele, matches('DNA')) %>%
    gather('sample_id', 'counts', matches('DNA|RNA')) %>%
    group_by(.data$variant_id, .data$allele, .data$sample_id) %>%
    nest(.key = 'count_dat') %>%
    filter(map_lgl(.data$count_dat, ~!all(.x$counts == 0))) %>% # some borderline barcodes are 0 in some samples
    mutate(nb_fit = parallel::mclapply(.data$count_dat, fit_nb, mc.cores = n_cores),
           converged = map_lgl(.data$nb_fit, ~.x$convergence == 0))

  if (!all(dna_nb_fits$converged)) {
    warning(paste0(sum(!dna_nb_fits$converged),
                   ' out of ',
                   nrow(dna_nb_fits),
                   ' (', round(sum(!dna_nb_fits$converged) / nrow(dna_nb_fits) * 100, digits = 3), '%)',
                   ' DNA-allele-samples failed to converge when fitting negative binomial parameters. A small fraction (<5%) failing is acceptable.'))
  }

  dna_nb_fits %<>%
    filter(.data$converged) %>%
    left_join(sample_depths, by = 'sample_id') %>%
    mutate(depth_adj_mu_est = map2_dbl(.data$nb_fit, .data$depth_factor, ~.x$par[1] / .y),
           phi_est = map_dbl(.data$nb_fit, ~.x$par[2]),
           acid_type = factor(stringr::str_extract(.data$sample_id, 'DNA|RNA'))) %>%
    filter(.data$phi_est < quantile(.data$phi_est, probs = .995)) # cut out severely underdispersed alleles

  dna_gamma_prior = dna_nb_fits %>%
    summarise(mu_prior = list(fit_gamma(.data$depth_adj_mu_est)),
              phi_prior = list(fit_gamma(.data$phi_est))) %>%
    ungroup %>%
    gather('prior_type', 'prior', matches('prior')) %>%
    mutate(alpha_est = map_dbl(.data$prior, ~.x$par[1]),
           beta_est = map_dbl(.data$prior, ~.x$par[2]),
           acid_type = 'DNA') # doesn't line up :(

  #### Fit RNA prior ----
  mean_dna_abundance = mpra_data %>%
    select(.data$variant_id, .data$allele, .data$barcode, matches('DNA')) %>%
    gather('sample_id', 'counts', matches('DNA|RNA')) %>%
    left_join(sample_depths,
              by = 'sample_id') %>%
    mutate(depth_adj_count = .data$counts / .data$depth_factor) %>%
    group_by(.data$barcode) %>%
    summarise(mean_depth_adj_count = mean(.data$depth_adj_count))

  print('Fitting marginal RNA mean priors...')
  rna_m_prior = mpra_data %>%
    select(.data$variant_id, .data$allele, .data$barcode, matches('RNA')) %>%
    filter(.data$barcode %in% well_represented$barcode) %>%
    gather('sample_id', 'counts', matches('DNA|RNA')) %>%
    left_join(sample_depths , by = 'sample_id') %>%
    left_join(mean_dna_abundance, by = 'barcode') %>%
    mutate(count_remnant = .1 + .data$counts / .data$depth_factor / .data$mean_depth_adj_count) %>% # the variability of the count after accounting for depth and DNA input
    group_by(.data$allele) %>%
    summarise(mu_prior = list(fit_gamma(.data$count_remnant))) %>%
    gather('prior_type', 'prior', matches('prior')) %>%
    mutate(alpha_est = map_dbl(.data$prior, ~.x$par[1]),
           beta_est = map_dbl(.data$prior, ~.x$par[2])) %>%
    mutate(acid_type = 'RNA')

  # the +.1 is to give some non-infinite log-density to 0's. See plot below. It seems to work well.

  print('Fitting marginal RNA dispersion priors...')
  rna_p_prior = mpra_data %>%
    select(.data$variant_id, .data$allele, .data$barcode, matches('RNA')) %>%
    filter(.data$barcode %in% well_represented$barcode) %>%
    gather('sample_id', 'counts', matches('DNA|RNA')) %>%
    group_by(.data$allele, .data$barcode) %>%
    summarise(mean_est = mean(.data$counts),
              var_est = var(.data$counts),
              size_guess = .data$mean_est^2 / (.data$var_est - .data$mean_est)) %>%
    filter(.data$size_guess > 0 & is.finite(.data$size_guess)) %>% # negative size guess = var < mean = underdispersed
    filter(.data$size_guess < quantile(.data$size_guess,
                                       probs = .99)) %>% # HUGE size guess = underdispersed, cut out barcodes that are TOO consistent i.e. underdispersed
    summarise(phi_prior = list(fit_gamma(.data$size_guess))) %>%
    gather('prior_type', 'prior', matches('prior')) %>%
    mutate(alpha_est = map_dbl(.data$prior, ~.x$par[1]),
           beta_est = map_dbl(.data$prior, ~.x$par[2])) %>%
    mutate(acid_type = 'RNA')

  # There is room for improvement with this phi prior estimation. It disregards
  # sequencing sample depth V and cuts out observed larged phi values.
  # It still covers the vast majority of alleles (98%)

  tot_prior = bind_rows(dna_gamma_prior,
                        rna_m_prior,
                        rna_p_prior)

  return(tot_prior)

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
#' @param min_neighbors The minimum number of neighbors in annotation spcae that
#'   must contribute to prior estimation
#' @param kernel_fold_increase The amount to iteratively increase kernel width
#'   by when estimating conditional priors. Smaller values (closer to 1) will
#'   yield more refined priors but take longer.
#' @export
fit_cond_prior = function(mpra_data,
                          annotations,
                          n_cores = 1,
                          plot_rep_cutoff = TRUE,
                          rep_cutoff = .15,
                          min_neighbors = 30,
                          kernel_fold_increase = 1.3){

  # Input checks ----
  if(!all(mpra_data$variant_id %in% annotations$variant_id)){
    stop('There aren\'t annotations for each variant!')
  }

  if(any(duplicated(annotations$variant_id))){
    stop('There are duplicate annotations for some variants!')
  }

  #### DNA prior fitting ----

  print('Evaluating data depth/DNA representation properties...')
  sample_depths = get_sample_depths(mpra_data)

  print('Determining well-represented variants, see plot...')
  well_represented = get_well_represented(mpra_data,
                                          sample_depths,
                                          rep_cutoff = rep_cutoff,
                                          plot_rep_cutoff = plot_rep_cutoff)


  print('Fitting marginal DNA prior...')

  dna_nb_fits = mpra_data %>%
    filter(.data$barcode %in% well_represented$barcode) %>%
    select(.data$variant_id, .data$allele, matches('DNA')) %>%
    gather('sample_id', 'counts', matches('DNA|RNA')) %>%
    group_by(.data$variant_id, .data$allele, .data$sample_id) %>%
    nest(.key = 'count_dat') %>%
    filter(map_lgl(.data$count_dat, ~!all(.x$counts == 0))) %>% # some borderline barcodes are 0 in some samples
    mutate(nb_fit = parallel::mclapply(.data$count_dat, fit_nb, mc.cores = n_cores),
           converged = map_lgl(.data$nb_fit, ~.x$convergence == 0))

  if (!all(dna_nb_fits$converged)) {
    warning(paste0(sum(!dna_nb_fits$converged),
                   ' out of ',
                   nrow(dna_nb_fits),
                   ' (', round(sum(!dna_nb_fits$converged) / nrow(dna_nb_fits) * 100, digits = 3), '%)',
                   ' DNA-allele-samples failed to converge when fitting negative binomial parameters. A small fraction (<5%) failing is acceptable.'))
  }

  dna_nb_fits %<>%
    filter(.data$converged) %>%
    left_join(sample_depths, by = 'sample_id') %>%
    mutate(depth_adj_mu_est = map2_dbl(.data$nb_fit, .data$depth_factor, ~.x$par[1] / .y),
           phi_est = map_dbl(.data$nb_fit, ~.x$par[2]),
           acid_type = factor(stringr::str_extract(.data$sample_id, 'DNA|RNA'))) %>%
    filter(.data$phi_est < quantile(.data$phi_est,
                                    probs = .995)) # cut out severely underdispersed alleles

  dna_gamma_prior = dna_nb_fits %>%
    summarise(mu_prior = list(fit_gamma(.data$depth_adj_mu_est)),
              phi_prior = list(fit_gamma(.data$phi_est))) %>%
    ungroup %>%
    gather('prior_type', 'prior', matches('prior')) %>%
    mutate(alpha_est = map_dbl(.data$prior, ~.x$par[1]),
           beta_est = map_dbl(.data$prior, ~.x$par[2]),
           acid_type = 'DNA') # doesn't line up :(

  #### Fit conditional RNA prior ----

  mean_dna_abundance = mpra_data %>%
    select(.data$variant_id, .data$allele, .data$barcode, matches('DNA')) %>%
    gather('sample_id', 'counts', matches('DNA|RNA')) %>%
    left_join(sample_depths,
              by = 'sample_id') %>%
    mutate(depth_adj_count = .data$counts / .data$depth_factor) %>%
    group_by(.data$barcode) %>%
    summarise(mean_depth_adj_count = mean(.data$depth_adj_count))

  # generate annotation distance matrix
  dist_mat = generate_distance_matrix(annotations = annotations)
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
  print('Weighting variants in annotation space')
  prior_weights = mpra_data %>%
    select(.data$variant_id) %>%
    unique() %>%
    mutate(annotation_weights = parallel::mclapply(.data$variant_id, find_prior_weights,
                                                   scaled_annotations = scaled_annotations,
                                                   dist_mat = dist_mat,
                                                   min_dist_kernel = min_dist_kernel,
                                                   mc.cores = n_cores))


  # Perform weighted density estimation for each variant
  print('Fitting annotation-weighted distributions...')
  if (n_cores == 1) {
    rna_m_priors = prior_weights %>%
      mutate(variant_m_prior = map2(.data$variant_id, .data$annotation_weights,
                                    fit_one_m_prior,
                                    mpra_data = mpra_data,
                                    sample_depths = sample_depths,
                                    well_represented = well_represented,
                                    mean_dna_abundance = mean_dna_abundance))

    rna_p_priors = prior_weights %>%
      mutate(variant_p_prior = map2(.data$variant_id, .data$annotation_weights,
                                    fit_one_p_prior,
                                    mpra_data = mpra_data,
                                    sample_depths = sample_depths,
                                    well_represented = well_represented,
                                    mean_dna_abundance = mean_dna_abundance))
  } else if (n_cores > 1) {
    rna_m_priors = prior_weights %>%
      mutate(variant_m_prior = parallel::mcmapply(fit_one_m_prior,
                                                  .data$variant_id, .data$annotation_weights,
                                                  MoreArgs = list(mpra_data = mpra_data,
                                                                  sample_depths = sample_depths,
                                                                  well_represented = well_represented,
                                                                  mean_dna_abundance = mean_dna_abundance),
                                                  mc.cores = n_cores,
                                                  SIMPLIFY = FALSE))
    rna_p_priors = prior_weights %>%
      mutate(variant_p_prior = parallel::mcmapply(fit_one_p_prior,
                                                  .data$variant_id, .data$annotation_weights,
                                                  MoreArgs = list(mpra_data = mpra_data,
                                                                  sample_depths = sample_depths,
                                                                  well_represented = well_represented,
                                                                  mean_dna_abundance = mean_dna_abundance),
                                                  mc.cores = n_cores,
                                                  SIMPLIFY = FALSE))
  }

  rna_p_no_weights = rna_p_priors %>% select(-.data$annotation_weights)

  both_rna = rna_m_priors %>%
    left_join(rna_p_no_weights, by = 'variant_id')

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


#' Score an annotation source
#'
#' @description This function takes MPRA data and an annotation source (or
#'   alternatively pre-fit priors) and checks how well the conditional prior
#'   improves upon the marginal prior by a prior ratio at the maximum likelihood
#'   estimates. A conditional:marginal ratio > 1 indicates that the conditional
#'   prior does better. A higher fraction of alleles for which this is true
#'   indicates an improved conditional prior.
#' @param mpra_data a data frame of MPRA data
#' @param annotations a data frame containing annotations of the same alleles in mpra_data
#' @param nb_params an optional data frame of pre-fit maximum likelihood estimates for allele parameters
#' @param conditional_prior an optional data frame containing pre-fit conditional priors
#' @param marginal_prior an optional data frame containing a pre-fit marginal prior
#' @return A fraction between 0 and 1 indicating the number of parameter
#'   estimates that have higher density under the conditional prior than under
#'   the marginal prior.
#' @note Note that even with near-perfect annotations or conditional priors, the
#'   fraction of parameter estimate ratios above 1 will not approach 1. It will
#'   instead approach some difficult-to-determine limit defined by the noise in
#'   the assay. Thus this functionality should only be used to COMPARE
#'   annotation sources, rather than make discreate claims / measurements about
#'   individual annotation sources.
score_annotation = function(mpra_data,
                            annotations,
                            nb_params = NULL,
                            conditional_prior = NULL,
                            marginal_prior = NULL){

}

format_conditional_prior = function(given_id, cond_priors){
  dna_prior = cond_priors$dna_prior
  rna_prior = cond_priors$rna_priors %>%
    filter(variant_id == given_id) %>%
    select(-annotation_weights) %>%
    gather(prior_name, prior, matches('prior')) %>%
    unnest() %>%
    select(-variant_id, -prior_name)

  return(bind_rows(dna_prior, rna_prior))

}
