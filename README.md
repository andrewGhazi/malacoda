
<!-- README.md is generated from README.Rmd. Please edit that file -->
malacoda <img src="man/figures/logo.png" align="right" />
=========================================================

The goal of malacoda is to enable Bayesian analysis of high-throughput genomic assays like massively parallel reporter assays (MPRA) and CRISPR screens.

It uses a negative-binomial-based Bayesian model that offers numerous advantages over traditional null hypothesis significance testing based methods:

-   Models raw data - The model is fit directly to the input counts (MPRA barcodes or gRNAs)
-   The lack of transformations avoids discarding 0 counts as in traditional methods.
-   Prior information - Empirical priors are fit from the observed assay globally, enabling estimate shrinkage that reduces errors due to multiple testing
-   Informative annotations (such as DNase hypersensitivity estimates or gene scores) can be included to further refine the empirical priors by conditional density estimation.
-   The R interface provides clear and interpretable outputs and figures.

malacoda runs a custom Stan model for fast posterior evaluation.

Installation
------------

You can install malacoda from github with:

``` r
# install.packages("devtools")
devtools::install_github("andrewGhazi/malacoda")
```

Example - under construction
----------------------------

This is a basic example which shows you how to solve a common problem:

``` r
## basic example code
```