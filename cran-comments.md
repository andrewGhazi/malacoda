## Test environments  
* local Linux Mint 19, R 3.5.1 
* local macOS High Sierra 10.13.1, R 3.5.1

Windows is not a target platform for us, so Windows-related
errors/warnings/build-failures do not bother us as long as they are not a
problem for the CRAN maintainers.

## R CMD check results 
There were no ERRORs or WARNINGs.

There was 1 NOTE: 
* checking for GNU extensions in Makefiles ... NOTE GNU make is a
SystemRequirements.

GNU make was added as a system requirement by rstantools::rstan_package_skeleton(). Removing this requirement would likely interfere with the compilation of the Stan models in this package.

This is the first submission of this package.

## Downstream dependencies

There are no downstream dependencies for this package.

Thank you! :)
