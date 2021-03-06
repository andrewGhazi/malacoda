% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting_functions.R
\name{plot_prior_ratios}
\alias{plot_prior_ratios}
\title{Plot prior ratios}
\usage{
plot_prior_ratios(prior_ratios, x_limits = c(-1, 1), n_bins = 100)
}
\arguments{
\item{prior_ratios}{a data frame of prior ratios from \code{get_prior_ratios}}

\item{x_limits}{a length two vector for the x-limits of the histograms}

\item{n_bins}{number of bins in the histogram}
}
\description{
Visualize prior ratios by histogram
}
\examples{
plot_prior_ratios(ratios_example, n_bins = 30)
}
