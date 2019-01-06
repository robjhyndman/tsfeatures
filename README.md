
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tsfeatures

[![Pending
Pull-Requests](http://githubbadges.herokuapp.com/robjhyndman/tsfeatures/pulls.svg?style=flat)](https://github.com/robjhyndman/tsfeatures/pulls)
[![Travis-CI Build
Status](https://travis-ci.org/robjhyndman/tsfeatures.svg?branch=master)](https://travis-ci.org/robjhyndman/tsfeatures)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/tsfeatures)](https://cran.r-project.org/package=tsfeatures)
[![Downloads](http://cranlogs.r-pkg.org/badges/tsfeatures)](https://cran.r-project.org/package=tsfeatures)
[![Licence](https://img.shields.io/badge/licence-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)

The R package *tsfeatures* provides methods for extracting various
features from time series data.

## Installation

You can install the **stable** version on [R
CRAN](https://cran.r-project.org/package=tsfeatures).

``` r
install.packages('tsfeatures', dependencies = TRUE)
```

You can install the **development** version from
[Github](https://github.com/robjhyndman/tsfeatures) with:

``` r
# install.packages("devtools")
devtools::install_github("robjhyndman/tsfeatures")
```

## Usage

``` r
library(tsfeatures)
mylist <- list(sunspot.year, WWWusage, AirPassengers, USAccDeaths)
myfeatures <- tsfeatures(mylist)
myfeatures
#> # A tibble: 4 x 20
#>   frequency nperiods seasonal_period trend   spike linearity curvature
#>       <dbl>    <dbl>           <dbl> <dbl>   <dbl>     <dbl>     <dbl>
#> 1         1        0               1 0.125 2.10e-5      3.58      1.11
#> 2         1        0               1 0.985 3.01e-8      4.45      1.10
#> 3        12        1              12 0.989 2.12e-8     11.0       1.10
#> 4        12        1              12 0.796 9.67e-7     -2.13      2.85
#> # ... with 13 more variables: e_acf1 <dbl>, e_acf10 <dbl>, entropy <dbl>,
#> #   x_acf1 <dbl>, x_acf10 <dbl>, diff1_acf1 <dbl>, diff1_acf10 <dbl>,
#> #   diff2_acf1 <dbl>, diff2_acf10 <dbl>, seasonal_strength <dbl>,
#> #   peak <dbl>, trough <dbl>, seas_acf1 <dbl>
```

## License

This package is free and open source software, licensed under GPL-3.
