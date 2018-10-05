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
                   rna_depths = sample_depths %>% filter(grepl('RNA', sample_id)) %>% pull(depth_factor),
                   dna_depths = sample_depths %>% filter(grepl('DNA', sample_id)) %>% pull(depth_factor),
                   dna_m_a = priors %>% filter(acid_type == 'DNA', grepl('mu', prior_type)) %>% pull(alpha_est),
                   dna_m_b = priors %>% filter(acid_type == 'DNA', grepl('mu', prior_type)) %>% pull(beta_est),
                   dna_p_a = priors %>% filter(acid_type == 'DNA', !grepl('mu', prior_type)) %>% pull(alpha_est),
                   dna_p_b = priors %>% filter(acid_type == 'DNA', !grepl('mu', prior_type)) %>% pull(beta_est),
                   rna_m_a = priors %>% arrange(desc(allele)) %>% filter(acid_type == 'RNA', grepl('mu', prior_type)) %>% pull(alpha_est),
                   rna_m_b = priors %>% arrange(desc(allele)) %>% filter(acid_type == 'RNA', grepl('mu', prior_type)) %>% pull(beta_est),
                   rna_p_a = priors %>% arrange(desc(allele)) %>% filter(acid_type == 'RNA', !grepl('mu', prior_type)) %>% pull(alpha_est),
                   rna_p_b = priors %>% arrange(desc(allele)) %>% filter(acid_type == 'RNA', !grepl('mu', prior_type)) %>% pull(beta_est))

  sampler_res = rstan::sampling(stanmodels$bc_mpra_model,
                                data = data_list,
                                chains = n_chains,
                                warmup = n_warmup,
                                iter = n_per_chain,
                                cores = 1)

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
                      hdi_upper = ts_hdi_obj[2])

  if(!missing(ts_rope)) {
    ts_rope_mass = sum(ts_vec < ts_rope[2] & ts_vec > ts_rope[1]) / tot_samp
    res_df$ts_rope_mass = ts_rope_mass
  }

  return(res_df)

}

run_crispr_sampler = function() {

}
