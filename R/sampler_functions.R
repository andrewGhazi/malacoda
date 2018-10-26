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
                   n_ref = variant_dat %>% filter(tolower(allele) == 'ref') %>% nrow,
                   n_alt = variant_dat %>% filter(tolower(allele) != 'ref') %>% nrow,
                   ref_dna_counts = variant_dat %>% filter(tolower(allele) == 'ref') %>% select(matches('DNA')) %>% as.matrix,
                   alt_dna_counts = variant_dat %>% filter(tolower(allele) != 'ref') %>% select(matches('DNA')) %>% as.matrix,
                   ref_rna_counts = variant_dat %>% filter(tolower(allele) == 'ref') %>% select(matches('RNA')) %>% as.matrix,
                   alt_rna_counts = variant_dat %>% filter(tolower(allele) != 'ref') %>% select(matches('RNA')) %>% as.matrix,
                   rna_depths = depth_factors %>% filter(grepl('RNA', sample_id)) %>% pull(depth_factor),
                   dna_depths = depth_factors %>% filter(grepl('DNA', sample_id)) %>% pull(depth_factor),
                   dna_m_a = priors %>% filter(acid_type == 'DNA', grepl('mu', prior_type)) %>% pull(alpha_est),
                   dna_m_b = priors %>% filter(acid_type == 'DNA', grepl('mu', prior_type)) %>% pull(beta_est),
                   dna_p_a = priors %>% filter(acid_type == 'DNA', !grepl('mu', prior_type)) %>% pull(alpha_est),
                   dna_p_b = priors %>% filter(acid_type == 'DNA', !grepl('mu', prior_type)) %>% pull(beta_est),
                   rna_m_a = priors %>% arrange(desc(allele)) %>% filter(acid_type == 'RNA', grepl('mu', prior_type)) %>% pull(alpha_est),
                   rna_m_b = priors %>% arrange(desc(allele)) %>% filter(acid_type == 'RNA', grepl('mu', prior_type)) %>% pull(beta_est),
                   rna_p_a = priors %>% arrange(desc(allele)) %>% filter(acid_type == 'RNA', !grepl('mu', prior_type)) %>% pull(alpha_est),
                   rna_p_b = priors %>% arrange(desc(allele)) %>% filter(acid_type == 'RNA', !grepl('mu', prior_type)) %>% pull(beta_est))

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

sample_from_prior = function(prior_df){
  sim_df = prior_df %>%
    mutate(allele = case_when(tolower(allele) == 'ref' ~ 'ref',
                              tolower(allele) != 'ref' ~ 'alt'),
           draws = map2(alpha_est, beta_est, ~rgamma(5e4, shape = .x, rate = .y))) %>%
    select(allele, draws) %>%
    unnest %>%
    mutate(iter = rep(1:5e4, times = 2)) %>%
    spread(allele, draws) %>%
    mutate(sim_ts = alt - ref)

  return(sim_df)
}

summarise_prior_samples = function(sim_df){

  sim_summary = sim_df %>%
    summarise(mean_prior_ts = mean(sim_ts),
              sd_prior_ts = sd(sim_ts),
              prior_ts_hdi = list(coda::HPDinterval(coda::mcmc(sim_ts), prob = .99)),
              prior_lower_ts = prior_ts_hdi[[1]][1],
              prior_upper_ts = prior_ts_hdi[[1]][2],
              prior_is_func = !between(0, prior_lower_ts, prior_upper_ts))

  return(sim_summary)
}

summarise_one_prior = function(prior_df){
  prior_samples = sample_from_prior(prior_df)

  prior_summary = summarise_prior_samples(prior_samples)

  return(prior_summary)
}


summarise_cond_prior = function(cond_prior,
                                n_cores = 1){

  cond_prior$rna_priors %>%
    mutate(prior_sim_res = parallel::mclapply(variant_m_prior,
                                              summarise_one_prior,
                                              mc.cores = n_cores)) %>%
    unnest(... = prior_sim_res) %>%
    arrange(desc(abs(mean_prior_ts)))
}
