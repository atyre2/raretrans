
<!-- README.md is generated from README.Rmd. Please edit that file -->

# raretrans

[![Build
Status](https://travis-ci.com/atyre2/raretrans.svg?token=fZ6v1nEUhWRrRW3uj9TH&branch=master)](https://travis-ci.com/atyre2/raretrans)

[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/atyre2/raretrans?branch=master&svg=true)](https://ci.appveyor.com/project/atyre2/raretrans)

[![codecov](https://codecov.io/gh/atyre2/raretrans/branch/master/graph/badge.svg?token=NWFl8m4owW)](https://codecov.io/gh/atyre2/raretrans)

[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/raretrans)](https://cran.r-project.org/package=raretrans)

Functions to create matrix population models from a combination of data
on stage/age transitions and Bayesian prior information. This
compensates for structural problems caused by missing observatons of
rare transitions.

## Installation

raretrans is not currently available from CRAN, but you can install the
development version from github with:

``` r
# install.packages("devtools")
devtools::install_github("atyre2/raretrans")
```

If you want the vignettes, you have to add some additional
arguments.

``` r
devtools::install_github("atyre2/raretrans", build = TRUE, build_opts = c("--no-resave-data", "--no-manual"))
```
