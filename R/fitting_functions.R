
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

  if(missing(weights)){
    fn_to_min = function(ab_vec){
      -sum(dgamma(input_vec,
                  shape = ab_vec[1],
                  rate = ab_vec[2],
                  log = TRUE))
    }
  } else {
    fn_to_min = function(ab_vec){
      -sum(weights*dgamma(input_vec,
                          shape = ab_vec[1],
                          rate = ab_vec[2],
                          log = TRUE))
    }
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
#' @param annotations an optional data frame of annotations with identical
#'   variant_ids and an arbitrary number of functional annotations. If omitted,
#'   the prior for a given variant is influenced by all other variants in the
#'   assay equally.
#' @param group_df an optional data frame giving group identity by variant_id in mpra_data
#' @param priors optional objects provided by either fit_marg_prior() or
#'   fit_cond_prior.
#' @param out_dir path to output directory
#' @param save_nonfunctional logical indicating whether or not to save the
#'   sampler results for variants identified as non-functional
#' @param n_cores number of cores across which to parallelize variant MPRA
#'   samplers
#' @param n_chains number of MCMC chains to run in each sampler
#' @param tot_samp total number of MCMC draws to take, spread evenly across
#'   chains
#' @param n_warmup total number of warmup draws to take from each MCMC chain
#' @param vb_pass logical indicating whether to use a variational first pass
#' @param vb_prob numeric 0 - 1 indicating probability mass to use as a TS HDI
#'   for identifying "promising" candidates for MCMC followup
#' @param ts_hdi_prob probability mass to include in the highest density
#'   interval on transcription shift to call MPRA-functional variants.
#' @param ts_rope optional length 2 numeric vector describing the boundaries of
#'   the transcription shift region of practical equivalence
#' @param rep_cutoff a representation cutoff quantile (0 to 1)
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
#'   If \code{priors} is provided, any annotations input will be ignored. This
#'   can be useful when you want to fit models again without having to spend
#'   time re-fitting the priors.
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
#'   \code{vb_pass} indicates whether to use a first pass variational check to
#'   see if a given variant is worth running the MCMC sampler. It does this by
#'   checking if a 40% HDI on the variational transcription shift posterior
#'   excludes 0. This speeds up posterior evaluation considerably, but gives
#'   approximate results. If \code{vb_pass} is set to FALSE, all variants get
#'   MCMC.
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
#'   Barcodes below the \code{rep_cutoff} quantile of representation in the DNA
#'   pools are discarded.
#'
#' @return a data frame with a row for each variant_id that specificies the
#'   posterior mean TS, upper and lower HDI bounds, a binary call of functional
#'   or non-functional, and other appropriate outputs
#' @note Sampler results for individaul variants will be saved to the specified
#'   out_dir as they can be several megabytes each
#' @export
fit_mpra_model = function(mpra_data,
                          annotations = NULL,
                          group_df = NULL,
                          out_dir,
                          save_nonfunctional = FALSE,
                          priors = NULL,
                          n_cores = 1,
                          n_chains = 4,
                          tot_samp = 1e4,
                          n_warmup = 500,
                          vb_pass = TRUE,
                          vb_prob = .8,
                          ts_hdi_prob = .95,
                          ts_rope = NULL,
                          rep_cutoff = .15) {

  start_time = Sys.time()

  #### Input checks ----

  if(missing(ts_rope)){
    ts_rope = NULL
  }

  if (missing(mpra_data)) {
    stop('mpra_data is missing: You must provide MPRA data to fit a MPRA model!')
  }

  if (missing(out_dir)) {
    warning('out_dir is missing: Results will not be saved')
  }

  if (!dir.exists(out_dir)) {
    stop('specified out_dir does not exist')
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

  correct_columns = all(grepl('variant_id|allele|barcode|[DR]NA', names(mpra_data)))
  if (!correct_columns){
    stop('mpra_data columns must be: variant_id, allele, barcode, and DNA/RNA columns.\ndplyr::rename(), dplyr::select(), malacoda::count_barcodes(), and the tidyr package might be helpful for preparing your input.')
  }

  # Check that there are 2 alleles for each variant
  variant_allele_counts = mpra_data %>%
    select(.data$variant_id, .data$allele) %>%
    unique %>%
    dplyr::count(.data$variant_id)

  if (!all(variant_allele_counts$n == 2)){
    stop('Non-biallelic variants detected. The variant_id column should be the same for both alleles of a given variant.')
  }

  #### Initial cleanup ----

  mpra_data %<>%
    arrange(.data$variant_id)

  if (!missing(annotations)) {
    annotations %<>%
      arrange(.data$variant_id)
  }

  if (missing(priors)) {
    annotations_given = !is.null(annotations)
    # Fit priors
    message('No annotations provided, fitting marginal priors...')
    if (!annotations_given) {
      print('Fitting MARGINAL priors...')
      priors = fit_marg_prior(mpra_data,
                              n_cores = n_cores,
                              rep_cutoff = .15,
                              plot_rep_cutoff = TRUE)
    } else if (!is.null(group_df)) {
      message('Fitting group-wise priors...')
      priors = fit_grouped_prior(mpra_data,
                                 group_df = group_df,
                                 n_cores = n_cores,
                                 plot_rep_cutoff = TRUE,
                                 rep_cutoff = .15)

    } else {
      message('Fitting annotation-based conditional priors...')
      priors = fit_cond_prior(mpra_data,
                              annotations,
                              n_cores = n_cores,
                              plot_rep_cutoff = TRUE,
                              rep_cutoff = .15,
                              min_neighbors = 30,
                              kernel_fold_increase = 1.3)

      print('Conditional prior fitting done, ')
      save(priors,
           file = paste0(out_dir, 'conditional_prior.RData'))
    }
  } else {
    if (all(class(priors) == 'list')){
      message('Input prior class is list, interpreting as conditional priors.')
      annotations_given = TRUE
    } else if ('group_prior' %in% names(priors)) {
      message('Interpreting input prior as a grouped prior.')
      annotations_given = FALSE
    } else {
      message('Input prior class is not list, interpreting as marginal priors.')
      annotations_given = FALSE
    }
  }


  #### Run samplers ----
  print('Running model samplers...')

  n_rna = mpra_data %>% select(matches('RNA')) %>% ncol
  n_dna = mpra_data %>% select(matches('DNA')) %>% ncol
  sample_depths = get_sample_depths(mpra_data)

  well_represented = get_well_represented(mpra_data,
                                          sample_depths,
                                          rep_cutoff = rep_cutoff,
                                          plot_rep_cutoff = FALSE, # this will have been plotted in the prior fitting already if necessary
                                          verbose = FALSE)


  if (annotations_given) {

    # attach the conditional priors in the form expected by run_mpra_sampler
    sampler_input = mpra_data %>%
      filter(.data$barcode %in% well_represented$barcode) %>%
      group_by(.data$variant_id) %>%
      nest(.key = 'variant_data') %>%
      mutate(variant_prior = map(.data$variant_id,
                                 format_conditional_prior,
                                 cond_priors = priors))

    analysis_res =  sampler_input %>%
      mutate(sampler_stats = parallel::mcmapply(run_mpra_sampler,
                                                .data$variant_id, .data$variant_data, .data$variant_prior,
                                                MoreArgs = list(n_chains = n_chains,
                                                                n_warmup = n_warmup,
                                                                tot_samp = tot_samp,
                                                                n_rna = n_rna,
                                                                n_dna = n_dna,
                                                                depth_factors = sample_depths,
                                                                out_dir = out_dir,
                                                                save_nonfunctional = save_nonfunctional,
                                                                ts_hdi_prob = ts_hdi_prob,
                                                                ts_rope = ts_rope,
                                                                vb_pass = vb_pass,
                                                                vb_prob = vb_prob),
                                                mc.cores = n_cores,
                                                SIMPLIFY = FALSE)) %>%
      unnest(.data$sampler_stats,
             .drop = TRUE,
             .preserve = c(.data$variant_data, .data$variant_prior)) %>%
      arrange(desc(abs(.data$ts_post_mean)))
  } else if ('group_prior' %in% names(priors)) {
    # This block uses a grouped prior.
    analysis_res = mpra_data %>%
      filter(.data$barcode %in% well_represented$barcode) %>%
      group_by(.data$variant_id) %>%
      nest(.key = 'variant_data') %>%
      left_join(group_df, by = 'variant_id') %>%
      left_join(priors, by = 'group_id') %>% # give the grouped_prior by variant_id
      dplyr::rename('variant_prior' = 'group_prior') %>%
      mutate(sampler_stats = parallel::mcmapply(run_mpra_sampler,
                                                .data$variant_id, .data$variant_data, .data$variant_prior,
                                                MoreArgs = list(n_chains = n_chains,
                                                                n_warmup = n_warmup,
                                                                tot_samp = tot_samp,
                                                                n_rna = n_rna,
                                                                n_dna = n_dna,
                                                                depth_factors = sample_depths,
                                                                out_dir = out_dir,
                                                                save_nonfunctional = save_nonfunctional,
                                                                ts_hdi_prob = ts_hdi_prob,
                                                                ts_rope = ts_rope,
                                                                vb_pass = vb_pass,
                                                                vb_prob = vb_prob),
                                                mc.cores = n_cores,
                                                SIMPLIFY = FALSE)) %>%
      unnest(.data$sampler_stats,
             .drop = TRUE,
             .preserve = c(.data$variant_data, .data$variant_prior)) %>%
      arrange(desc(abs(.data$ts_post_mean)))

  } else {
    # This block uses the marg priors

    analysis_res = mpra_data %>%
      filter(.data$barcode %in% well_represented$barcode) %>%
      group_by(.data$variant_id) %>%
      nest(.key = 'variant_data') %>%
      mutate(variant_prior = list(priors)) %>% # give the same marg prior to every variant
      mutate(sampler_stats = parallel::mcmapply(run_mpra_sampler,
                                                .data$variant_id, .data$variant_data, .data$variant_prior,
                                                MoreArgs = list(n_chains = n_chains,
                                                                n_warmup = n_warmup,
                                                                tot_samp = tot_samp,
                                                                n_rna = n_rna,
                                                                n_dna = n_dna,
                                                                depth_factors = sample_depths,
                                                                out_dir = out_dir,
                                                                save_nonfunctional = save_nonfunctional,
                                                                ts_hdi_prob = ts_hdi_prob,
                                                                ts_rope = ts_rope,
                                                                vb_pass = vb_pass,
                                                                vb_prob = vb_prob),
                                                mc.cores = n_cores,
                                                SIMPLIFY = FALSE)) %>%
      unnest(.data$sampler_stats,
             .drop = TRUE,
             .preserve = c(.data$variant_data, .data$variant_prior)) %>%
      arrange(desc(abs(.data$ts_post_mean)))
  }

  save(analysis_res,
       file = paste0(out_dir, 'analysis_res.RData'))

  end_time = Sys.time()
  time_diff = end_time - start_time

  message(paste0('MPRA data for ', n_distinct(mpra_data$variant_id), ' variants analyzed in ',
                 round(digits = 3, end_time - start_time), ' ', attr(time_diff, 'units')))

  return(analysis_res)
}


#' Fit a Bayesian model of dropout CRISPR screen data
#'
#' @description This function fits a Bayesian model of survival/dropout CRISPR
#'   screen data. It uses a negative binomial to model the input and output
#'   counts of sgRNAs, adjusting the results appropriately to account for
#'   sequencing depth.
#' @param dropout_data a data frame of dropout data. See details for column
#'   requirements.
#' @param n_cores number of cores to utilize
#' @inheritParams fit_mpra_model
#' @details \code{dropout_data} requires the following columns:
#'   \itemize{\item{gene_id - gives a unique identifier for each gene}
#'   \item{sgRNA - an identifier for individual sgRNAs (usually the sgRNA
#'   sequence itself)} \item{input count columns - columns of sequencing counts
#'   of the input sgRNA library. Multiple columns for sequencing replicates are
#'   allowed (which require unique identifiers). Column names must contain the
#'   string "input".} \item{output count columns - columns of sequencing counts
#'   of sgRNAs in the output libraries.}}
#' @note Currently this function only supports marginal priors. If you want to
#'   use grouped/conditional priors, contact the malacoda developers.
fit_dropout_model = function(dropout_data,
                             out_dir,
                             n_cores = 1,
                             tot_samp = 1e4,
                             n_warmup = 500,
                             rep_cutoff = .15) {

  #### input checks ----


}
