## Resubmission
This is a resubmission of malacoda, which was previously released May 10 2019.

* As requested, we have added examples for the packages most important functions, most prominently fit_mpra_model(). Each row of results from fit_mpra_model() involves running an MCMC chain, so the example uses too-few posterior samples for only 3 variants (rows) for the sake of running quickly (~2s on my laptop). This is noted in a comment above the example, with recommendations to the user on how to improve the MCMC. The example_result object included in the package is there to provide a precise example in the vignette; it would take too long to run to be in an example.

## Test environments  
* local Linux Mint 19.1, R 3.6.0
* local macOS Mojave 10.14.1, R 3.6.0
* Ubuntu Linux 16.04 LTS, R-release, GCC via Rhub
* Fedora Linux, R-devel, clang, gfortran via Rhub

Windows is not a target platform for us, so Windows-related errors/warnings/build-failures do not bother us if they are not a problem for the CRAN maintainers.

## R CMD check results 
There were no ERRORs or WARNINGs.

There were 4 NOTEs: 
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Andrew Ghazi <arghazi@bcm.edu>'

New submission

Possibly mis-spelled words in DESCRIPTION:
  CRISPR (9:84, 9:144)
  MPRA (9:75, 9:127)
  barcodes (9:132)
  gRNAs (9:151)
  
This is a resubmission of this package (my first) after the first submission was rejected do to some initial problems.

The words CRISPR, MPRA, barcodes, and gRNAs are all spelled correctly. These seem to be false positives from the spell-checker.
  
* checking for GNU extensions in Makefiles ... NOTE GNU make is a
SystemRequirements.

GNU make was added as a system requirement by rstantools::rstan_package_skeleton(). Removing this requirement would likely interfere with the compilation of the Stan models in this package.

* on local Linux Mint 19.1, R 3.6.0: checking compilation flags used ... NOTE
  Compilation used the following non-portable flag(s):
    ‘-march=native’
    
The -march=native flag is recommended by the RStan team as part of our Stan installation. We have tried leaving off this flag and it doesn't seem to cause any problems, so we don't think its absence will cause issues.

* checking installed package size ... NOTE
  installed size is 41.7Mb
  sub-directories of 1Mb or more:
    doc    1.2Mb
    libs  39.2Mb
    
The large compiled libraries come from the package's Stan models and the associated pre-compiled C++ modules which are fairly large. There doesn't seem to be a way of reducing the size of these compiled Stan models, and they are essential to the function of the package.

## Downstream dependencies

There are no downstream dependencies for this package.

Thank you! :)
