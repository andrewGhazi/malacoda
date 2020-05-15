# This readme file contains a description of the various .R files in the
# malacoda package. None of these files should be utilized directly unless you
# want to get into messing with the source code. All of the functions / help
# documentation they define can be accessed simply by loading the package with
# library(malacoda). Some functions are only used by the package internally and
# are not exported. Other .R files located elsewhere in the package directory
# are part of the structure laid out by rstantools::rstan_package_skeleton
# (which is used to lay out R packages that use Stan).

# annotation_analysis.R - functionality for analyzing and visualizing the impact
# of annotations used to define conditional priors (over marginal priors)

# data.R - roxygen2 documentation for the package data objects

# fitting_functions.R - NB and Gamma distribution fitting functions, malacoda
# model fitting functions for MPRA and CRISPR datasets

# importing_functions.R - used to define functions to import from dependencies

# malacoda-package.R - used in package documentation

# plotting_functions.R - all visualization functions.

# prior_fitting.R - functions for fitting marginal, conditional, and grouped priors

# qc_functions.R - miscellaneous functions for checking sample correlations, converting FASTQ files to MPRA counts, and decoding sequencing errors using the FreeBarcodes package.

# sampler_functions.R - functions to call Stan

# stanmodels.R - a function from rstantools::rstan_package_skeleton() to make the Stan models accessible from elsewhere in the package

# traditional_analysis_functions.R - functions to run the traditional activity analysis based on t-tests

# zzz.R -  a function from rstantools::rstan_package_skeleton() to make the Stan models accessible from elsewhere in the package
