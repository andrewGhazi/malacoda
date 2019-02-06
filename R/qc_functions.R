
#' Get MPRA sample correlations
#'
#' @description Compute the correlations of MPRA sequencing samples. This can be
#'   a helpful as a QC diagnostic metric. Libraries of the same type (DNA or RNA) should
#'   correlate highly with one another
#' @param mpra_data a data frame of mpra data
#' @return a data frame of pairwise correlations
#' @details the output is ready to be presented to plot_mpra_correlations
#' @export
get_sample_correlations = function(mpra_data){
  mpra_data %>%
    select(matches('[DR]NA')) %>%
    cor() %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = 'sample_1') %>%
    gather("sample_2", 'correlation', -.data$sample_1)
}

#' Cut out sequence lines
#'
#' @description Given a trimmed and quality filtered fastq, cut out just the
#'   sequence lines, write them out, then exit.
#' @inheritParams count_barcodes
#' @inheritParams count_barcodes_in_fastq
cut_out_seqs = function(trimmed_fastq,
                        temp_dir){

  fastq_path_split = stringr::str_split(trimmed_fastq,
                                        pattern = '/',
                                        simplify = TRUE)

  fastq_name = gsub('\\.fastq','', fastq_path_split[1,ncol(fastq_path_split)])

  output_path = paste0(temp_dir, 'seq_only/seq_only_', fastq_name, '.txt')
  seq_only_cmd = paste0('sed -n \'2~4p\' ', trimmed_fastq, ' > ', output_path)
  system(seq_only_cmd)

  return(TRUE)
}

#' Match barcodes for a single fastq
#'
#' @description Count barcode abundance, write out the counts of unmatched
#'   barcodes, return a data frame of counts of the given barcodes
#' @param trimmed_fastq path to a trimmed fastq
#' @inheritParams count_barcodes
count_barcodes_in_fastq = function(trimmed_fastq,
                                   barcode_allele_df,
                                   temp_dir){

  fastq_path_split = stringr::str_split(trimmed_fastq,
                                        pattern = '/',
                                        simplify = TRUE)

  fastq_name = gsub('\\.fastq','', fastq_path_split[1,ncol(fastq_path_split)])

  output_path = paste0(temp_dir, 'seq_only/seq_only_', fastq_name, '.txt')

  # If someone has an absurdly large amount of sequencing this might be slow. Might have to add data.table as a dependency and do it with that, it's > 10x faster

  if (requireNamespace('data.table', quietly = TRUE)){
    barcode = NULL # to avoid "no visible binding..." note in the check
    .N = NULL
    barcode_observations = tibble::as_tibble(data.table::fread(output_path,
                                                               col.names = 'barcode')[,.(n = .N), by = 'barcode'][order(barcode)])
  } else {
    barcode_observations = readr::read_tsv(output_path,
                                           col_names = 'barcode') %>%
      dplyr::count(.data$barcode)

  }

  unmatched = barcode_observations %>%
    filter(!(.data$barcode %in% barcode_allele_df$barcode))

  # These will be available to the user if they set keep_temp = TRUE in count_barcodes()
  readr::write_tsv(unmatched,
            path = paste0(temp_dir, 'seq_only/', fastq_name, '_unmatched_barcode_counts.tsv'))

  sample_name = gsub('trimmed_filtered_', '', fastq_name)

  barcode_counts = barcode_observations %>%
    right_join(barcode_allele_df, by = 'barcode') %>%
    mutate(n = replace(.data$n, is.na(.data$n), 0)) %>%
    dplyr::select('bc_id', 'barcode', 'n') %>%
    magrittr::set_colnames(c('bc_id', 'barcode', sample_name))

  return(barcode_counts)


  # This doesn't work for more than 1021 barcodes without OS-level ulimit changes
  # bc_splitter_cmd = paste0('cat ', trimmed_fastq, ' | ',
  #                          'fastx_barcode_splitter.pl ',
  #                          '--bcfile ', temp_dir, 'bc_file.txt ', # this will have been written out just prior to this function being called
  #                          '--prefix ', temp_dir, 'seq_only/', fastq_name, '/ ',
  #                          '--bol ',
  #                          '--exact ',
  #                          '--suffix .txt')

}

#' Trim and Filter a FASTQ
#'
#' @description Use the FASTX-Toolkit to trim out the barcodes from the FASTQs by position, then quality filter them.
#'
#' @param fastq path to a fastq file
#' @inheritParams count_barcodes
trim_and_filter = function(fastq,
                           temp_dir,
                           quality_cutoff,
                           bc_start,
                           bc_end){

  fastq_path_split = stringr::str_split(fastq,
                                        pattern = '/',
                                        simplify = TRUE)

  fastq_name = fastq_path_split[1,ncol(fastq_path_split)]

  trim_filter_cmd = paste0('fastx_trimmer ',
                           '-f ', bc_start, ' ',
                           '-l ', bc_end, ' ',
                           '-i', fastq,
                           ' | fastq_quality_filter ',
                           '-q ', quality_cutoff, ' ',
                           ' > ', temp_dir, 'trimmed_filtered_fastqs/trimmed_filtered_', fastq_name)

  system(trim_filter_cmd)

  return(TRUE)

}

#' Count barcodes
#'
#' @description Given a mapping between barcodes and alleles and a directory of
#'   MPRA sequencing results, count the abundance of each barcode and return a
#'   data frame of counts.
#' @param barcode_allele_df a data frame giving a mapping between barcodes and
#'   alleles
#' @param fastq_dir a character string giving the path to a directory of MPRA
#'   FASTQ results
#' @param temp_dir a character string giving the path to a directory to use for
#'   temporary intermediate files
#' @param keep_temp a logical indicating whether to avoid discarding the
#'   temporary intermediate files
#' @param bc_start an integer indicating the position in the read where the
#'   barcode starts
#' @param bc_end an integer indicating the position in the read where the
#'   barcode ends
#' @param quality_cutoff an integer indicating the quality cutoff to be used
#' @param n_cores an integer giving the number of cores to use in parallel
#' @details The \code{barcode_allele_df} should have two columns: \code{bc_id}
#'   giving a unique identifier to each barcode and \code{barcode} which gives the
#'   barcode. The bc_id should ideally be descriptive and denote which allele of
#'   which SNP it is associated with. tidyr::unite can be convenient for
#'   preparing this input.
#'
#'   The quality cutoff applies to ALL bases in the barcode. That is, if every
#'   base in the barcode is not at or above the given cutoff, the read is
#'   discarded. The default cutoff, 30, is fairly aggressive, so if you're
#'   losing too many reads you can try lowering it.
#'
#'   The output returns a column of counts for each FASTQ in the input
#'   directory. The column names are taken from the filenames of the fastqs
#'   without the filetype extensions. The usual MPRA practice is to have one
#'   FASTQ per transfection as well as multiple replicates taken from the
#'   plasmid library -- if this is not the case in your experiment you'll need
#'   to do some form of experiment-structure-appropriate aggregation before
#'   proceeding to downstream analysis.
#'
#'   Unmatched reads are also written out to files included in the temporary
#'   intermediates under seq_only/. These can be helpful if you want to regain
#'   reads that would otherwise be discarded by using the error-correctable
#'   barcodes from the freebarcodes package available through
#'   \href{https://github.com/andrewGhazi/mpradesigntools}{mpradesigntools}. If
#'   you want to use these make sure to set \code{keep_temp = TRUE}.
#'
#' @return a data frame of counts by file in the input directory
#' @note This function requires an installation of
#'   \href{http://hannonlab.cshl.edu/fastx_toolkit/index.html}{FASTX-Toolkit} by
#'   the Hannon Lab (it's easy ot install with apt-get) and sed (which comes as
#'   part of any standard Unix-like OS).
#'
#'   The temporary intermediates are comparable in size to the input FASTQ
#'   files, so make sure you have enough disk space available.
#' @export
count_barcodes = function(barcode_allele_df,
                          fastq_dir,
                          temp_dir,
                          keep_temp = FALSE,
                          bc_start,
                          bc_end,
                          quality_cutoff = 30,
                          n_cores = 1,
                          verbose = TRUE){

  #### Input checks ----
  fastx_toolkit_installed = length(system('which fastx_trimmer', intern = TRUE)) == 1 # There's probably a better way to check for this...

  if (!fastx_toolkit_installed){
    stop('It seems the FASTX-Toolkit is not installed. See http://hannonlab.cshl.edu/fastx_toolkit/download.html for installation details.')
  }

  if (!dir.exists(fastq_dir)) {
    stop('The FASTQ directory specified doesn\'t exist!')
  }

  if (!dir.exists(temp_dir)){
    message('The temporary directory specified doesn\'t exist! Creating with dir.create()')
    dir.create(temp_dir)
  }

  if (fastq_dir == temp_dir){
    stop('Please specify a temp_dir other than where the FASTQs are.')
  }

  temp_dir_ends_in_slash = grepl('/$', temp_dir)
  if (!temp_dir_ends_in_slash){
    temp_dir = paste0(temp_dir, '/')
  }

  fastq_dir_ends_in_slash = grepl('/$', fastq_dir)
  if (!fastq_dir_ends_in_slash){
    fastq_dir = paste0(fastq_dir, '/')
  }

  if (any(grepl('[^A-Za-z0-9_]', barcode_allele_df$bc_id))){
    # This might not actually be
    stop('Some barcode identifiers may be ill-formatted. Please restrict IDs to only alphanumeric and underscore characters.')
  }

  if (length(list.files(temp_dir)) != 0){
    stop('Please specify an empty directory as the temp_dir')
  }

  fastqs = list.files(gsub('/$', '', fastq_dir),
                      pattern = '.fastq$',
                      full.names = TRUE)

  if (length(fastqs) == 0) {
    stop('No .fastq files found in fastq_dir! Use gunzip *.fastq.gz if you have gzip-compressed files.')
  }

  if (quality_cutoff != 30 & verbose) {
    message(paste0('Using non-default quality cutoff: Q >= ', quality_cutoff))
  }

  #### Trim and quality filter with the FASTX-Toolkit ----

  dir.create(paste0(temp_dir, 'trimmed_filtered_fastqs/'))

  if (verbose) {message('Trimming & quality filtering fastqs...')}

  tf_cmd = parallel::mclapply(fastqs,
                    trim_and_filter,
                    mc.cores = n_cores,
                    quality_cutoff = quality_cutoff,
                    bc_start = bc_start,
                    bc_end = bc_end,
                    temp_dir = temp_dir)

  dir.create(paste0(temp_dir, 'seq_only/'))

  trimmed_fastqs = list.files(paste0(temp_dir, 'trimmed_filtered_fastqs'),
                              full.names = TRUE)

  if (verbose) {message('Extracting sequence lines...')}

  cut_cmd = parallel::mclapply(trimmed_fastqs,
                     cut_out_seqs,
                     mc.cores = n_cores,
                     temp_dir = temp_dir)

  #### On to the counting ----

  if (verbose) {message('Counting barcodes...')}

  mpra_data = parallel::mclapply(trimmed_fastqs,
                       count_barcodes_in_fastq,
                       mc.cores = n_cores,
                       barcode_allele_df = barcode_allele_df,
                       temp_dir = temp_dir) %>%
    purrr::reduce(.f = dplyr::inner_join, by = c('bc_id', 'barcode'))

  if (verbose) {message('Counting done. If no/few barcodes are found, check if the barcodes in your reads are the reverse complement of how you ordered them. If so, you may need to use DNAString(), reverseComplement(), and toString() from the Biostrings package on your input')}

  #### Cleanup if required ----

  if (!keep_temp){
    unlink(temp_dir, recursive = TRUE)
  }

  #### Return ----

  return(mpra_data)
}
