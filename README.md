
<!-- README.md is generated from README.Rmd. Please edit that file -->
malacoda <img src="man/figures/logo.png" align="right" title="Evil tail"/>
==========================================================================

The goal of malacoda is to enable Bayesian analysis of high-throughput genomic assays like massively parallel reporter assays (MPRA) and CRISPR screens.

It uses a negative-binomial-based Bayesian model that offers numerous advantages over traditional null hypothesis significance testing based methods:

-   Models raw data - The model is fit directly to the input counts (MPRA barcodes or gRNAs)
    -   The lack of transformations avoids discarding 0 counts as in traditional methods.
-   Prior information - Empirical priors are fit from the observed assay globally, enabling estimate shrinkage that reduces errors due to multiple testing
    -   Informative annotations (such as DNase hypersensitivity estimates or gene scores) can be included to further refine the empirical priors by conditional density estimation.
-   The R interface provides clear and interpretable outputs and figures.

Other features include:

-   custom Stan models for fast posterior evaluation
-   variational Bayes support through `rstan::vb` that allows for quick first pass checks
-   Annotation checking - quantitatively evaluate how much a given genomic annotation source improves empirical prior estimation by prior ratios

Installation
------------

You can install malacoda from github with:

``` r
# install.packages("devtools")
devtools::install_github("andrewGhazi/malacoda")
```

Example
-------

This is a basic example which shows you how to solve a common problem:

``` r
## basic example code
library(malacoda)
marg_prior = fit_marg_prior(umpra_example)
fit_mpra_model(mpra_data = umpra_example,
               out_dir = '/path/to/outputs/',
               priors = marg_prior,
               n_cores = getOption('mc.cores', 2L),
               vb_pass = TRUE,
               save_nonfunctional = TRUE)
```

This will fit the model to each input in the assay (using some example data from [Ulirsch et al., Cell, 2016](https://www.ncbi.nlm.nih.gov/pubmed/27259154)) using a marginal prior, save the outputs for each variant at the specified directory, and return a data frame of summary statistics for each variant, including binary calls of functional/non-functional, posterior means on activity levels & transcription shift.

More sophisticated analyses that use annotations to create informative priors for more sensitive analysis will be described in an upcoming vignette. The required functionality to do so currently exists with `fit_cond_prior`. Other features like annotation checking and traditional NHST analysis will also be explained in the vignette.

Example output
--------------

The sampler outputs for variants are saved in the provided output directory. These are stanfit objects, hence they can be visualized using all the tools provided in packages like [bayesplot](http://mc-stan.org/users/interfaces/bayesplot).

`malacoda` also provides one plotting function of its own (with more planned), `posterior_beeswarm()`. This plots the traditional activity measurements as points in a beeswarm plot along with violins for posteriors on means for each allele. Optional colors can help diagnose unwanted sample-specific bias. That is, all the colors should be mixed within each allele, indicating that activity measurements are not influenced by sample.

![An example activity beeswarm with overlaid activity mean posteriors](man/figures/posterior_beeswarm_example.png)

Contact
-------

Please contact me through Github DM or my BCM email address if you use the package, have feature requests / comments, or want a hexagon sticker!
