% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prior_fitting.R
\name{fit_cond_prior}
\alias{fit_cond_prior}
\title{Fit a informative conditional prior}
\usage{
fit_cond_prior(
  mpra_data,
  annotations,
  n_cores = 1,
  plot_rep_cutoff = TRUE,
  rep_cutoff = 0.15,
  min_neighbors = 100,
  kernel_fold_increase = 1.4142,
  verbose = TRUE
)
}
\arguments{
\item{mpra_data}{a data frame of mpra data}

\item{annotations}{a data frame of annotations for the same variants in
mpra_data}

\item{n_cores}{number of cores to parallelize across}

\item{plot_rep_cutoff}{logical indicating whether to plot the representation
cutoff used}

\item{rep_cutoff}{fraction indicating the depth-adjusted DNA count quantile
to use as the cutoff}

\item{min_neighbors}{The minimum number of neighbors in annotation space that
must contribute to prior estimation}

\item{kernel_fold_increase}{The amount to iteratively increase kernel width
by when estimating conditional priors. Smaller values (closer to 1) will
yield more refined priors but take longer.}

\item{verbose}{logical indicating whether to print messages}
}
\value{
A list of two data frames. The first is for the DNA and the second is
  by-variant RNA priors.
}
\description{
Use informative annotations to bias prior estimation towards
  alleles that show similar annotations in the provided annotation space.
}
\details{
The empirical prior returned by this object is "conditional" in the
  sense that the prior estimation weights are conditional on the annotations.

  The DNA prior is still estimated marginally because the annotations should
  not be able to provide any information on the DNA inputs (which are
  presumably only affected by the preparation of the oligonucleotide library
  at the vendor).

  The RNA prior is estimated from the RNA observations of
  other variants in the assay that are nearby in annotation space. A
  multivariate t distribution centered on the variant in question is used to
  weight all other variants in the assay. It is initialized with a very small
  width, and if there are fewer than \code{min_neighbors} that provide
  substantial input to the prior, the width is iteratively increased by a
  factor of \code{kernel_fold_increase} until that condition is satisfied.
  This prevents the prior estimation for variants in sparse regions of
  annotation space from being influenced too heavily by their nearest
  neighbors.
}
\examples{
cond_prior = fit_cond_prior(mpra_data = umpra_example,
                            annotations = u_deepsea,
                            n_cores = 1,
                            rep_cutoff = .15,
                            plot_rep_cutoff = TRUE,
                            min_neighbors = 5)
}
