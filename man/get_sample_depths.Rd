% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/traditional_analysis_functions.R
\name{get_sample_depths}
\alias{get_sample_depths}
\title{Get Sample Depths}
\usage{
get_sample_depths(mpra_data, depth_multiplier = 1e+06)
}
\arguments{
\item{mpra_data}{a dataframe of MPRA data}

\item{depth_multiplier}{a numeric to divide through the depths to make the
numbers more easily interpretable}
}
\description{
Computes the sum of all barcode counts
}
\note{
The \code{depth_multiplier} input has no effect on the downstream
  analysis as it divides out because both RNA and DNA counts are normalized
  by this same factor. It simply sets the scale of depth factors at an
  easily-readable range.
}
\examples{
get_sample_depths(umpra_example)
# Data from non-subsampled datasets with typical sequencing depths will typically
# show depth factors larger than this example, in the 10-100 range.
}
