## Test environments  
* local Linux Mint 19, R 3.5.1 
* local Linux Mint 19, R-devel (11/05/2018)
* local macOS High Sierra 10.13.1, R 3.5.1

Windows is not a target platform for us, so Windows-related
errors/warnings/build-failures do not bother us as long as they are not a
problem for the CRAN maintainers.

## R CMD check results 
There were no ERRORs or WARNINGs.

There were 2 NOTEs: 
* checking for GNU extensions in Makefiles ... NOTE GNU make is a
SystemRequirements.
* checking CRAN incoming feasibility ... NOTE
Maintainer: ‘Andrew Ghazi <arghazi@bcm.edu>’

New submission

Found the following (possibly) invalid URLs:
  URL: https://www.bloodgenes.org/RBC_MPRA
    From: man/umpra_example.Rd
    Status: 404
    Message: Not Found

GNU make was added as a system requirement by rstantools::rstan_package_skeleton(). Removing this requirement would likely interfere with the compilation of the Stan models in this package.

This is the first submission of this package.

The mentioned URL is currently down but is the original source of the example data, so we have left it for posterity. The original scientific publication, Ulirsch et al., Cell, 2016, is also linked in the Readme.

## Downstream dependencies

There are no downstream dependencies for this package.

Thank you! :)
