
<!-- README.md is generated from README.Rmd. Please edit that file -->
malacoda <img src="man/figures/logo.png" align="right" />
=========================================================

The goal of malacoda is to enable Bayesian analysis of high-throughput genomic assays like massively parallel reporter assays (MPRA) and CRISPR screens.

It uses a negative-binomial-based Bayesian model that offers numerous advantages over traditional null hypothesis significance testing based methods:

-   A Bayesian generative model under the hood. This itself confers several advantages:
      -   Raw data modelling - The model is fit directly to the input counts (MPRA barcodes or gRNAs)
      -   No data transformations - The lack of transformations (e.g. a log-ratio) avoids discarding 0 counts as in traditional methods.
-   Prior information - Empirical priors are fit from the observed assay globally, enabling estimate shrinkage that reduces errors due to multiple testing
-   Informative annotations (such as DNase hypersensitivity estimates or gene scores) can be included to further refine the empirical priors by conditional density estimation.

Other features include:

-   custom Stan models for fast posterior evaluation
    -   variational first pass with `rstan::vb` for fast "promising candidate" checks
-   Annotation checking - quantitatively evaluate how much a given genomic annotation source improves empirical prior estimation by prior ratios
-   Configurable plotting functions 
-   Methods for traditional NHST-based MPRA analysis

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

Example output
--------------

![An example activity beeswarm with overlaid activity mean posteriors](man/figures/posterior_beeswarm_example.png)
