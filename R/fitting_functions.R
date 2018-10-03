# global imports
#' @importFrom magrittr %>%
#' @importFrom magrittr %<>%
#' @import dplyr
#' @import purrr
#' @import tidyr
empty_fun = function(){

}



#' @title Fit a negative binomial
#'
#' @description Fit a negative binomial distribution to a set of counts
#'
#' @param input_dat a data frame with one column called "counts"
#' @param nb_init optional vector of mean, dispersion initialization point
fit_nb = function(input_dat,
                  nb_init = c(10, 1)) {

  input_vec = input_dat$counts # feeding it a data frame with one column of counts

  fn_to_min = function(param_vec){
    # param_vec[1] nb mean
    # param_vec[2] nb size
    -sum(stats::dnbinom(input_vec,
                        mu = param_vec[1],
                        size = param_vec[2],
                        log = TRUE))
  }

  stats::nlminb(start = nb_init,
                objective = fn_to_min,
                lower = rep(.Machine$double.xmin, 2))

}

fit_gamma = function(input_vec,
                     weights = NULL,
                     gamma_init = c(1, 1)) {

  # if fails to fit, stop() with error message suggesting different inits

  fn_to_min = function(ab_vec){
    -sum(dgamma(input_vec,
                shape = ab_vec[1],
                rate = ab_vec[2],
                log = TRUE))
  }

  stats::nlminb(start = gamma_init,
                objective = fn_to_min,
                lower = rep(.Machine$double.xmin, 2))

}

#' Fit a Bayesian MPRA model
#'
#' @description This function fits a negative-binomial based Bayesian model to
#'   MPRA data. Optional annotations can be included to allow for more
#'   informative conditional priors.
#'
#' @param mpra_data a data frame of MPRA data with 1 column called variant_id,
#'   an allele column, a barcode column, and additional columns per sequencing
#'   sample. Each row is for a single barcode.
#' @param annotations a optional data frame of annotations with identical
#'   variant_ids and an arbitrary number of functional annotations. If omitted,
#'   the prior for a given variant is influenced by all other variants in the
#'   assay equally.
#' @param out_dir path to output directory
#' @param save_nonfunctional logical indicating whether or not to save the
#'   sampler results for variants identified as non-functional
#' @param n_cores number of cores across which to parallelize variant MPRA
#'   samplers
#' @param n_chains number of MCMC chains to run in each sampler
#' @param tot_samp total number of MCMC draws to take, spread evenly across
#'   chains
#' @param n_warmup total number of warmup draws to take from each MCMC chain
#' @param ts_hdi_prob probability mass to include in the highest density
#'   interval on transcription shift to call MPRA-functional variants.
#' @param ts_rope optional length 2 numeric vector describing the boundaries of
#'   the transcription shift region of practical equivalence
#'
#' @details \code{mpra_data} must contain the following groups of columns:
#'   \itemize{ \item{variant_id} \item{allele - either 'ref' or 'alt'}
#'   \item{barcode - a unique index sequence for that row (ideally the same
#'   barcode used in the assay)} \item{at least one column of MPRA counts whose
#'   column name(s) matches 'DNA'} \item{at least one column of MPRA counts
#'   whose column name(s) matches 'RNA'} }
#'
#'   \code{annotations} must contain the same variant_id's used in mpra_data.
#'   Additional columns are used as informative predictors: when estimating the
#'   priors for one variant, other variants with similar annotations will be
#'   upweighted in the prior-fitting process.
#'
#'   Sampler results will be saved to out_dir. By default, only the sampler
#'   results for variants identified as MPRA-functional will be saved. This
#'   behavior can be changed by setting \code{save_nonfunctional} to TRUE.
#'
#'   We've set the sampler parameters (n_chains to n_warmup) to values that work
#'   reasonably well at reasonable speeds for typical MPRA data on typical
#'   hardware. Final analyses and/or models fit to larger MPRA experiments will
#'   likely want to increase n_chains and tot_samp considerably to ensure
#'   precise convergence.
#'
#'   \code{ts_rope} can be used to define a "Region Of Practical Equivalence"
#'   for transcription shift. This is some small-ish region around 0 where
#'   observed posterior samples are "practically equivalent" to 0. Enabling this
#'   option will return the fraction of transcription shift posterior samples
#'   that fall within the defined ROPE along with the usual model outputs. If
#'   this fraction is small, one can say that there is very little posterior
#'   belief that the variant's transcription shift is practically equivalent to
#'   0. If used, the user must be cognizant of defining the region in accordance
#'   with observed noise and effect size levels. Note that the output ROPE
#'   fractions ARE NOT p-values.
#'
#' @return a data frame with a row for each variant_id that specificies the
#'   posterior mean TS, upper and lower HDI bounds, a binary call of functional
#'   or non-functional, and other appropriate outputs
#' @note Sampler results for individaul variants will be saved to the specified
#'   out_dir as they can be several megabytes each
#' @export
fit_mpra_model = function(mpra_data,
                          annotations = NULL,
                          out_dir,
                          save_nonfunctional = FALSE,
                          n_cores = 1,
                          n_chains = 4,
                          tot_samp = 1e4,
                          n_warmup = 500,
                          ts_hdi_prob = .95,
                          ts_rope = NULL) {

  start_time = Sys.time()

  #### Input checks ----
  if (missing(mpra_data)) {
    stop('mpra_data is missing: You must provide MPRA data to fit a MPRA model!')
  }

  if (missing(out_dir)) {
    warning('out_dir is missing: Results will not be saved')
  }

  if (!dir.exists(out_dir)) {
    stop('specified out_dir does not exist')
  }

  if (any(duplicated(mpra_data$variant_id))) {
    stop('variant_id entries in mpra_data are not all unique!')
  }

  if (!missing(annotations)) {
    if (!all(mpra_data$variant_id %in% annotations$variant_id)) {
      stop('Some mpra_data$variant_id\'s missing from annotations')
    }
  }

  if (ts_hdi_prob < 0 | ts_hdi_prob > 1) {
    stop('ts_hdi_prob must be between 0 and 1!')
  }

  # make sure the out_directory ends in a slash, if not, add it
  dir_ends_in_slash = grepl('/$', out_dir)
  if (!dir_ends_in_slash){
    out_dir = paste0(out_dir, '/')
  }

  #### Initial cleanup ----

  mpra_data %<>%
    arrange(variant_id)

  if (!missing(annotations)) {
    annotations %<>%
      arrange(variant_id)
  }

  #### Fit priors ----
  pring('Fitting priors...')
  if (is.null(annotations)) {
    pring('Fitting MARGINAL priors...')
    priors = fit_marg_prior(mpra_data,
                            n_cores = n_cores,
                            rep_cutoff = .15,
                            plot_rep_cutoff = TRUE)
  } else {
    print('Fitting annotation-based conditional priors...')
    priors = fit_cond_prior(mpra_data, annotations, n_cores = n_cores)
  }

  #### Run samplers ----
  print('Running model samplers...')

  n_rna = mpra_data %>% select(matches('RNA')) %>% ncol
  n_dna = mpra_data %>% select(matches('DNA')) %>% ncol
  sample_depths = mpra_data %>%
    gather(sample_id, counts, matches('DNA|RNA')) %>%
    group_by(sample_id) %>%
    summarise(depth_factor = sum(counts) / 1e6)

  analysis_res = mpra_data %>%
    group_by(variant_id) %>%
    nest(.key = variant_dat) %>%
    mutate(sampler_stats = mcmapply(run_mpra_sampler,
                                    variant_id, variant_dat,
                                    MoreArgs = list(priors = priors,
                                                    n_chains = n_chains,
                                                    n_warmup = n_warmup,
                                                    tot_samp = tot_samp,
                                                    n_rna = n_rna,
                                                    n_dna = n_dna,
                                                    depth_factors = sample_depths,
                                                    out_dir = out_dir,
                                                    save_nonfunctional = save_nonfunctional,
                                                    ts_hdi_prob = ts_hdi_prob,
                                                    ts_rope = ts_rope),
                                    mc.cores = n_cores,
                                    SIMPLIFY = FALSE)) %>%
    unnest(... = sampler_stats) %>%
    arrange(desc(abs(ts_post_mean)))

  end_time = Sys.time()
  print(paste0('MPRA data for ', n_distinct(mpra_data$variant_id), ' variants analyzed in ',
               round(digits = 3, end_time - start_time), ' seconds'))
  return(analysis_res)
}

fit_crispr_model = function(crispr_data) {

}
