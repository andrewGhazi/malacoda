#' Example MPRA annotation set
#'
#' DeepSea v0.94 K562 DNase hypersensitivity estimates for the variants in Ulirsch et
#' al., Cell, 2016.
#' @format a data frame containing variant_id and annotation columns
#' @source \url{http://deepsea.princeton.edu/}
"u_deepsea"

#' Example MPRA counts from Ulirsch et al., 2016
#'
#' MPRA counts for 25 randomly chosen variants from Ulirsch et al., Cell, 2016.
#' Original URL from the manuscript seems to be down as of 10/2/2018.
#'
#' @format a data frame giving variant_id, allele identifier, barcode, and
#'   sequencing counts for the two DNA libraries and 6 control RNA libraries
"umpra_example"

#' Example output of fit_mpra_model
#'
#' Example needed for the vignette
#' @format a data frame of mpra modelling results
"example_result"

#' Example activity measurements
#'
#' Example activities for the variants in umpra_example
"activities_example"

#' Example FASTQ data
#'
#' Example FASTQ entries, used for the example in count_barcodes(). Taken from
#' Ulirsch et al., Cell, 2016.
"fastq_examples"

#' Example barcode_allele_df
#'
#' Example barcode_allele_df, used for the example in count_barcodes(). Adapted
#' from Ulirsch et al., Cell, 2016.
"barcode_by_id"

#' Example DeepSea annotations
#'
#' Example DeepSea annotations for the example MPRA data.
"u_deepsea"

#' Example posterior
#'
#' Example posterior for chr6:135426558 1/2 variant in Ulirsch et al., Cell,
#' 2016. This example used too-few MCMC samples for the sake of keeping the file
#' size small.
"example_posterior"

#' Example marginal prior
#'
#' Example marginal prior derived from MPRA data from Ulirsch et al., Cell,
#' 2016.
"marg_prior_example"

#' Example conditional prior
#'
#' Example conditional prior for the example variants taken from from Ulirsch et
#' al., Cell, 2016. This prior is the full, accurate prior using all 7000+
#' variants from that assay, so conditionl priors derived from the subsetted
#' example data included in malacoda will be different.
"cond_prior_example"

#' Example CRISPR dropout screen counts
#'
#' Example dropout screen counts for 50 randomly selected genes from "Optimized
#' libraries for CRISPR-Cas9 genetic screens with multiple modalities", Doench
#' et al., Nature Communications, 2018
"dropout_example"
