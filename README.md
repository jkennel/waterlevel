
waterlevel [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) [![](https://travis-ci.org/jkennel/waterlevel.svg?branch=master)](https://travis-ci.org/jkennel/waterlevel) [![Coverage Status](https://img.shields.io/codecov/c/github/jkennel/waterlevel/master.svg)](https://codecov.io/github/jkennel/waterlevel?branch=master) [![](https://www.r-pkg.org/badges/version/waterlevel?color=green)](https://cran.r-project.org/package=waterlevel) [![](http://cranlogs.r-pkg.org/badges/grand-total/waterlevel?color=green)](https://cran.r-project.org/package=waterlevel)
==========================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================

The **waterlevel** package is designed to simplify the analysis of well water level data in **R**. Currently the focus is on barometric/loading (BE/LE) efficiency, but the package is under active development..

Installation
------------

``` r
remotes::install_github("jkennel/waterlevel")
```

Static Barometric/Loading efficiency
------------------------------------

See the package vignettes for examples the methods.

### Static BE/LE methods

-   be\_clark
-   be\_rahi
-   be\_least\_squares
-   be\_acworth
-   be\_high\_low
-   be\_ratio
-   be\_visual

### Response methods

-   Barometric response functions
    -   Regular lags
    -   Irregular lags
    -   Distributed lags
-   Harmonic analysis
    -   Specified frequencies
    -   Specified Earth tide constituents *earthtide* package
-   Frequency response function

Future work
-----------

-   Decrease dependencies
-   fft methods
-   Add parameter estimation methods
-   Add analysis methods (slug testing, pumping tests)
-   Increase test coverage
-   Increase memory and computational efficiency
-   Gap filling methods (interpolation)
-   Add option to use parallel non-FFT distributed lag method
-   Distribute lag with subset
