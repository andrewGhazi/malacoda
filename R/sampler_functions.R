run_mpra_sampler = function(variant_id, variant_dat, variant_prior,
                            n_chains,
                            tot_samp,
                            n_warmup,
                            n_rna,
                            n_dna,
                            depth_factors,
                            out_dir,
                            save_nonfunctional,
                            ts_hdi_prob,
                            vb_pass = TRUE,
                            vb_prob = .8,
                            ts_rope = NULL) {



  priors = variant_prior
  n_per_chain = ceiling((tot_samp + n_chains * n_warmup) / n_chains)

  data_list = list(n_rna_samples = n_rna,
                   n_dna_samples = n_dna,
                   n_ref = variant_dat %>% filter(tolower(.data$allele) == 'ref') %>% nrow,
                   n_alt = variant_dat %>% filter(tolower(.data$allele) != 'ref') %>% nrow,
                   ref_dna_counts = variant_dat %>% filter(tolower(.data$allele) == 'ref') %>% select(matches('DNA')) %>% as.matrix,
                   alt_dna_counts = variant_dat %>% filter(tolower(.data$allele) != 'ref') %>% select(matches('DNA')) %>% as.matrix,
                   ref_rna_counts = variant_dat %>% filter(tolower(.data$allele) == 'ref') %>% select(matches('RNA')) %>% as.matrix,
                   alt_rna_counts = variant_dat %>% filter(tolower(.data$allele) != 'ref') %>% select(matches('RNA')) %>% as.matrix,
                   rna_depths = depth_factors %>% filter(grepl('RNA', .data$sample_id)) %>% pull(.data$depth_factor),
                   dna_depths = depth_factors %>% filter(grepl('DNA', .data$sample_id)) %>% pull(.data$depth_factor),
                   dna_m_a = priors %>% filter(.data$acid_type == 'DNA', grepl('mu', .data$prior_type)) %>% pull(.data$alpha_est),
                   dna_m_b = priors %>% filter(.data$acid_type == 'DNA', grepl('mu', .data$prior_type)) %>% pull(.data$beta_est),
                   dna_p_a = priors %>% filter(.data$acid_type == 'DNA', !grepl('mu', .data$prior_type)) %>% pull(.data$alpha_est),
                   dna_p_b = priors %>% filter(.data$acid_type == 'DNA', !grepl('mu', .data$prior_type)) %>% pull(.data$beta_est),
                   rna_m_a = priors %>% arrange(desc(.data$allele)) %>% filter(.data$acid_type == 'RNA', grepl('mu', .data$prior_type)) %>% pull(.data$alpha_est),
                   rna_m_b = priors %>% arrange(desc(.data$allele)) %>% filter(.data$acid_type == 'RNA', grepl('mu', .data$prior_type)) %>% pull(.data$beta_est),
                   rna_p_a = priors %>% arrange(desc(.data$allele)) %>% filter(.data$acid_type == 'RNA', !grepl('mu', .data$prior_type)) %>% pull(.data$alpha_est),
                   rna_p_b = priors %>% arrange(desc(.data$allele)) %>% filter(.data$acid_type == 'RNA', !grepl('mu', .data$prior_type)) %>% pull(.data$beta_est))

  if(vb_pass){
    vb_res = rstan::vb(stanmodels$bc_mpra_model,
                       data = data_list)

    vb_hdi = rstan::extract(vb_res,
                            pars = 'transcription_shift')$transcription_shift %>%
      as.matrix %>%
      coda::mcmc() %>%
      coda::HPDinterval(prob = vb_prob)

    if(!between(0, vb_hdi[1], vb_hdi[2])){
      sampler_res = rstan::sampling(stanmodels$bc_mpra_model,
                                    data = data_list,
                                    chains = n_chains,
                                    warmup = n_warmup,
                                    iter = n_per_chain,
                                    cores = 1)
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
                                  cores = 1)
    note = 'mcmc used for posterior evaluation'
  }

  ts_vec = rstan::extract(sampler_res,
                          pars = 'transcription_shift')$transcription_shift
  ref_act = rstan::extract(sampler_res,
                           pars = 'ref_act')$ref_act
  alt_act = rstan::extract(sampler_res,
                           pars = 'alt_act')$alt_act

  ts_hdi_obj = coda::HPDinterval(coda::mcmc(as.matrix(ts_vec)),
                                 prob = ts_hdi_prob)

  is_functional = !between(0, ts_hdi_obj[1], ts_hdi_obj[2])

  dir_ends_in_slash = grepl('/$', out_dir)
  if (!dir_ends_in_slash){
    out_dir = paste0(out_dir, '/')
  }

  if (save_nonfunctional | is_functional){
    save(sampler_res,
         file = paste0(out_dir, variant_id, '.RData'))
  }

  res_df = data_frame(ts_post_mean = mean(ts_vec),
                      ref_act_post_mean = mean(ref_act),
                      alt_act_post_mean = mean(alt_act),
                      is_functional = is_functional,
                      ts_hdi = list(ts_hdi_obj),
                      hdi_lower = ts_hdi_obj[1],
                      hdi_upper = ts_hdi_obj[2],
                      note = note)

  if(!is.null(ts_rope)) {
    ts_rope_mass = sum(ts_vec < ts_rope[2] & ts_vec > ts_rope[1]) / tot_samp
    res_df$ts_rope_mass = ts_rope_mass
  }

  return(res_df)

}

run_crispr_sampler = function() {

}

#' Sample from prior
#'
#' @description Randomly draw activities levels and transcription shifts from
#'   the prior
#' @param prior_df a data frame of the prior for one variant
#' @param n_iter number of prior simulation draws to pull
#' @details This function draws \code{n_iter} draws from the prior for one
#'   variant.
#' @return a data data frame of activity draws and the simulated transcription
#'   shift at each draw
#' @export
sample_from_prior = function(prior_df, n_iter){
  sim_df = prior_df %>%
    mutate(allele = case_when(tolower(.data$allele) == 'ref' ~ 'ref',
                              tolower(.data$allele) != 'ref' ~ 'alt'),
           draws = map2(.data$alpha_est, .data$beta_est,
                        ~rgamma(n_iter, shape = .x, rate = .y))) %>%
    select(.data$allele, .data$draws) %>%
    unnest %>%
    mutate(iter = rep(1:n_iter, times = 2)) %>%
    spread(.data$allele, .data$draws) %>%
    mutate(sim_ts = .data$alt - .data$ref)

  return(sim_df)
}

#' Summarise prior samples
#'
#' @description Compute summary statistics for input prior simulation draws
#' @param sim_df a data frame of simulation draws
#' @return a row data frame with columns of summary statistics
#' @export
summarise_prior_samples = function(sim_df){

  sim_summary = sim_df %>%
    summarise(mean_prior_ts = mean(.data$sim_ts),
              sd_prior_ts = sd(.data$sim_ts),
              prior_ts_hdi = list(coda::HPDinterval(coda::mcmc(.data$sim_ts), prob = .99)),
              prior_lower_ts = .data$prior_ts_hdi[[1]][1],
              prior_upper_ts = .data$prior_ts_hdi[[1]][2],
              prior_is_func = !between(0, .data$prior_lower_ts, .data$prior_upper_ts))

  return(sim_summary)
}

#' Summarise one RNA prior
#'
#' @description Sample from one RNA prior and compute summary statistics
#' @param prior_df a single RNA prior
#' @param n_iter number of prior draws used in simulation
#' @return a data frame with summary statistics for the input prior
#' @export
summarise_one_prior = function(prior_df,
                               n_iter){
  prior_samples = sample_from_prior(prior_df, n_iter = n_iter)

  prior_summary = summarise_prior_samples(prior_samples)

  return(prior_summary)
}

#' Summarise a conditional prior
#'
#' @description This function runs Monte Carlo simulations for each variant's
#'   conditional mean priors.
#' @param cond_prior a conditional prior object
#' @param n_cores number of cores for parallelization
#' @param n_iter number of prior simulation draws
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
#' @export
summarise_cond_prior = function(cond_prior,
                                n_cores = 1,
                                n_iter = 5e4){

  cond_prior$rna_priors %>%
    mutate(prior_sim_res = parallel::mclapply(.data$variant_m_prior,
                                              summarise_one_prior,
                                              mc.cores = n_cores,
                                              n_iter = n_iter)) %>%
    unnest(... = .data$prior_sim_res) %>%
    select(-.data$prior_sim_res) %>%
    arrange(desc(abs(.data$mean_prior_ts)))
}
