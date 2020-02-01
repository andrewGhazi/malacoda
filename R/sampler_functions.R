run_mpra_sampler = function(variant_id, variant_data, variant_prior,
                            n_chains,
                            tot_samp,
                            n_warmup,
                            n_rna,
                            n_dna,
                            depth_factors,
                            out_dir = NULL,
                            save_nonfunctional,
                            ts_hdi_prob,
                            vb_pass = TRUE,
                            vb_prob = .8,
                            ts_rope = NULL,
                            verbose = TRUE) {

  # This fits the malacoda biallelic MPRA model (i.e. the main one) for ONE variant.


  priors = variant_prior
  n_per_chain = ceiling((tot_samp + n_chains * n_warmup) / n_chains)

  if (verbose){
    refresh_setting = max(n_per_chain / 10, 1)
  } else {
    refresh_setting = 0
  }

  ref_data = variant_data %>% filter(tolower(.data$allele) == 'ref')
  alt_data = variant_data %>% filter(tolower(.data$allele) != 'ref')

  data_list = list(n_rna_samples = n_rna,
                   n_dna_samples = n_dna,
                   n_ref = ref_data %>% nrow,
                   n_alt = alt_data %>% nrow,
                   ref_dna_counts = ref_data %>% select(matches('DNA')) %>% as.matrix,
                   alt_dna_counts = alt_data %>% select(matches('DNA')) %>% as.matrix,
                   ref_rna_counts = ref_data %>% select(matches('RNA')) %>% as.matrix,
                   alt_rna_counts = alt_data %>% select(matches('RNA')) %>% as.matrix,
                   rna_depths = depth_factors %>% filter(grepl('RNA', .data$sample_id)) %>% pull(.data$depth_factor),
                   dna_depths = depth_factors %>% filter(grepl('DNA', .data$sample_id)) %>% pull(.data$depth_factor),
                   dna_m_a = priors %>% filter(.data$acid_type == 'DNA', grepl('mu', .data$prior_type)) %>% pull(.data$alpha_est),
                   dna_m_b = priors %>% filter(.data$acid_type == 'DNA', grepl('mu', .data$prior_type)) %>% pull(.data$beta_est),
                   dna_p_a = priors %>% filter(.data$acid_type == 'DNA', !grepl('mu', .data$prior_type)) %>% pull(.data$alpha_est),
                   dna_p_b = priors %>% filter(.data$acid_type == 'DNA', !grepl('mu', .data$prior_type)) %>% pull(.data$beta_est),
                   rna_m_a = priors %>% arrange(tolower(.data$allele) != 'ref') %>% filter(.data$acid_type == 'RNA', grepl('mu', .data$prior_type)) %>% pull(.data$alpha_est),
                   rna_m_b = priors %>% arrange(tolower(.data$allele) != 'ref') %>% filter(.data$acid_type == 'RNA', grepl('mu', .data$prior_type)) %>% pull(.data$beta_est),
                   rna_p_a = priors %>% arrange(tolower(.data$allele) != 'ref') %>% filter(.data$acid_type == 'RNA', !grepl('mu', .data$prior_type)) %>% pull(.data$alpha_est),
                   rna_p_b = priors %>% arrange(tolower(.data$allele) != 'ref') %>% filter(.data$acid_type == 'RNA', !grepl('mu', .data$prior_type)) %>% pull(.data$beta_est))

  if (vb_pass) {
    # If vb_pass is TRUE, run a variational first pass for the sake of a quick
    # check. In practice the distributions seem pretty well fit by the VB
    # approximations. I have Rmd on this, I should put it together for a
    # supplement to the manuscript.

    vb_res = rstan::vb(stanmodels$bc_mpra_model,
                       data = data_list,
                       tol_rel_obj = .0005)

    vb_hdi = (vb_res@sim$samples[[1]][['alt_act']] - vb_res@sim$samples[[1]][['ref_act']]) %>%
      as.matrix %>%
      coda::mcmc() %>%
      coda::HPDinterval(prob = vb_prob)

    if(!between(0, vb_hdi[1], vb_hdi[2])){
      # Check if the VB result indicates the variant is "worthy" of MCMC

      sampler_res = rstan::sampling(stanmodels$bc_mpra_model,
                                    data = data_list,
                                    chains = n_chains,
                                    warmup = n_warmup,
                                    iter = n_per_chain,
                                    cores = 1,
                                    verbose = verbose,
                                    refresh = refresh_setting)
      note = 'mcmc used for posterior evaluation'
    } else{
      sampler_res = vb_res
      note = 'VB used for posterior evaluation'
    }

  } else {
    sampler_res = rstan::sampling(stanmodels$bc_mpra_model,
                                  data = data_list,
                                  chains = n_chains,
                                  warmup = n_warmup,
                                  iter = n_per_chain,
                                  cores = 1,
                                  verbose = verbose,
                                  refresh = refresh_setting)
    note = 'mcmc used for posterior evaluation'
  }

  #### Compute some output quantities ----
  ts_vec = rstan::extract(sampler_res,
                          pars = 'transcription_shift')$transcription_shift
  ref_act = rstan::extract(sampler_res,
                           pars = 'ref_act')$ref_act
  alt_act = rstan::extract(sampler_res,
                           pars = 'alt_act')$alt_act

  ts_hdi_obj = coda::HPDinterval(coda::mcmc(matrix(ts_vec)),
                                 prob = ts_hdi_prob)

  # rope_alpha = 1 - ts_hdi_prob

  ts_rope_mass = sum(ts_vec < ts_rope[2] & ts_vec > ts_rope[1]) / tot_samp
  excludes_zero = !between(0, ts_hdi_obj[1], ts_hdi_obj[2])

  is_functional = excludes_zero # & (ts_rope_mass < rope_alpha)
  # ^ Whether or not to include the rope mass cutoff part could be made into a
  # user-accessible argument.

  #### Save the output ----
  if(!is.null(out_dir)){
    dir_ends_in_slash = grepl('/$', out_dir)
    if (!dir_ends_in_slash){
      out_dir = paste0(out_dir, '/')
    }
  }

  if ((save_nonfunctional | is_functional) && !is.null(out_dir)){
    save(sampler_res,
         file = paste0(out_dir, variant_id, '.RData'))
  }

  #### Compile summary data frame ----
  res_df = tibble(ts_post_mean = mean(ts_vec),
                  ref_act_post_mean = mean(ref_act),
                  alt_act_post_mean = mean(alt_act),
                  is_functional = is_functional,
                  ts_hdi = list(ts_hdi_obj),
                  hdi_lower = ts_hdi_obj[1],
                  hdi_upper = ts_hdi_obj[2],
                  ts_rope_mass = ts_rope_mass,
                  note = note)

  return(res_df)

}

#' Sample from prior
#'
#' @description Randomly draw activities levels and transcription shifts from
#'   the prior
#' @param prior_df a data frame of the prior for one variant
#' @param n_samp number of prior simulation draws to pull
#' @details This function draws \code{n_samp} draws from the prior for one
#'   variant.
#' @return a data data frame of activity draws and the simulated transcription
#'   shift at each draw
#' @examples
#' sample_from_prior(marg_prior_example, n_samp = 1000)
#' @export
sample_from_prior = function(prior_df, n_samp){
  sim_df = prior_df %>% filter(.data$prior_type == 'mu_prior') %>%
    mutate(allele = case_when(tolower(.data$allele) == 'ref' ~ 'ref',
                              tolower(.data$allele) != 'ref' ~ 'alt'),
           draws = map2(.data$alpha_est, .data$beta_est,
                        ~rgamma(n_samp, shape = .x, rate = .y))) %>%
    select(.data$allele, .data$draws) %>%
    unnest(c(.data$draws)) %>%
    filter(!is.na(.data$allele)) %>%
    mutate(iter = rep(1:n_samp, times = 2)) %>%
    spread(.data$allele, .data$draws) %>%
    mutate_if(is.double, log) %>%
    mutate(sim_ts = .data$alt - .data$ref) %>%
    dplyr::select(-.data$iter)

  return(sim_df)
}

#' Summarise prior samples
#'
#' @description Compute summary statistics for input prior simulation draws
#' @param sim_df a data frame of simulation draws
#' @param ts_hdi_prob probability mass to include in the highest density
#'   interval on transcription shift to call MPRA-functional variants
#' @details You can get simulated prior draws using \code{sample_from_prior()}
#' @return a row data frame with columns of summary statistics
#' @examples
#' prior_draws = sample_from_prior(marg_prior_example, n_samp = 1000)
#' summarise_prior_samples(prior_draws)
#' @export
summarise_prior_samples = function(sim_df,
                                   ts_hdi_prob = .95){

  sim_summary = sim_df %>%
    summarise(mean_prior_ts = mean(.data$sim_ts),
              sd_prior_ts = sd(.data$sim_ts),
              prior_ts_hdi = list(coda::HPDinterval(coda::mcmc(.data$sim_ts), prob = ts_hdi_prob)),
              prior_lower_ts = .data$prior_ts_hdi[[1]][1],
              prior_upper_ts = .data$prior_ts_hdi[[1]][2],
              prior_is_func = !between(0, .data$prior_lower_ts, .data$prior_upper_ts))

  return(sim_summary)
}

#' Summarise one RNA prior
#'
#' @description Sample from one RNA prior and compute summary statistics
#' @param prior_df a single RNA prior
#' @param n_samp number of prior draws used in simulation
#' @return a data frame with summary statistics for the input prior
summarise_one_prior = function(prior_df,
                               n_samp){
  prior_samples = sample_from_prior(prior_df, n_samp = n_samp)

  prior_summary = summarise_prior_samples(prior_samples)

  return(prior_summary)
}

#' Summarise a conditional prior
#'
#' @description This function runs Monte Carlo simulations for each variant's
#'   conditional mean priors.
#' @param cond_prior a conditional prior object
#' @param n_cores number of cores for parallelization
#' @param n_samp number of prior simulation draws
#' @return a data frame of RNA priors with prior simulation summary statistics
#' @details This is used to ensure that the prior is not TOO specific. While we
#'   want the prior distribution to accurately reflect prior beliefs on a given
#'   variant's effects, it would be undesirable for the prior to pre-specify
#'   that the variant is functional. If the \code{prior_is_func} column in the
#'   output of this function is all FALSE, this means that no variant is \emph{a
#'   priori} functional.
#'
#'   If one does get a conditional prior that has variants that are \emph{a
#'   priori} functional, this can be addressed by increasing the
#'   \code{min_neighbors} argument in \code{fit_cond_prior}. See the
#'   documentation of that function for details.
#'
#'   Currently this function only simulates the mean parameters as in practice
#'   the dispersion parameters don't vary systematically between alleles and/or
#'   variants.
#' @note simulations for individual variants can be obtained with
#'   \code{summarise_one_prior()} and \code{sample_from_prior()}
#' @examples
#' summarise_cond_prior(cond_prior_example, n_samp = 1000, n_cores = 1)
#' @export
summarise_cond_prior = function(cond_prior,
                                n_cores = 1,
                                n_samp = 5e4){

  cond_prior$rna_priors %>%
    mutate(prior_sim_res = parallel::mclapply(.data$variant_m_prior,
                                              summarise_one_prior,
                                              mc.cores = n_cores,
                                              n_samp = n_samp)) %>%
    unnest(cols = .data$prior_sim_res) %>%
    arrange(desc(abs(.data$mean_prior_ts)))
}

run_dropout_sampler = function(gene_id, gene_data, gene_prior,
                               n_chains = 4,
                               tot_samp,
                               n_warmup,
                               depth_factors,
                               out_dir = NULL) {

  # prepare input ----

  input_counts = gene_data %>%
    dplyr::select(sort(matches('input')))
  output_counts = gene_data %>%
    dplyr::select(sort(matches('output')))
  # ^ We sort the columns so they'll line up with the depth factors below (which
  # are also sorted)

  input_depths = depth_factors %>%
    dplyr::filter(grepl('input', .data$sample_id)) %>%
    dplyr::arrange(.data$sample_id)
  output_depths = depth_factors %>%
    dplyr::filter(grepl('output', .data$sample_id)) %>%
    dplyr::arrange(.data$sample_id)

  prior_mat = gene_prior %>%
    dplyr::select(-matches('prior')) %>%
    tibble::column_to_rownames('param_type') %>%
    as.matrix

  data_list = list(n_grna = nrow(gene_data),
                   n_in = ncol(input_counts),
                   n_out = ncol(output_counts),
                   in_counts = as.matrix(input_counts),
                   out_counts = as.matrix(output_counts),
                   in_depths = as.array(input_depths$depth_factor),
                   out_depths = as.array(output_depths$depth_factor),
                   in_mean_a = prior_mat['input', 'mean_alpha'],
                   in_mean_b = prior_mat['input', 'mean_beta'],
                   in_size_a = prior_mat['input', 'size_alpha'],
                   in_size_b = prior_mat['input', 'size_beta'],
                   out_mean_a = prior_mat['output', 'mean_alpha'],
                   out_mean_b = prior_mat['output', 'mean_beta'],
                   out_size_a = prior_mat['output', 'size_alpha'],
                   out_size_b = prior_mat['output', 'size_beta'])

  n_per_chain = ceiling((tot_samp + n_chains * n_warmup) / n_chains)
  # call sampler ----

  sampler_res = rstan::sampling(object = stanmodels$dropout_model,
                                iter = n_per_chain,
                                data = data_list,
                                chains = n_chains,
                                warmup = n_warmup,
                                cores = 1)

  if(!is.null(out_dir)){
    save(sampler_res,
         file = paste0(out_dir, gene_id, '.RData'))
  }

  # compute summary statistics ----

  summary_df = rstan::summary(sampler_res)$summary %>%
    as.data.frame %>%
    rownames_to_column(var = 'parameter') %>%
    as_tibble

  # return summary ----

  return(summary_df)
}
