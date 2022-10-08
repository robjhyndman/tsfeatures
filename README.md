
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tsfeatures

<!-- badges: start -->

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/tsfeatures)](https://cran.r-project.org/package=tsfeatures)
[![Downloads](http://cranlogs.r-pkg.org/badges/tsfeatures)](https://cran.r-project.org/package=tsfeatures)
[![Licence](https://img.shields.io/badge/licence-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)
[![R build
status](https://github.com/robjhyndman/tsfeatures/workflows/R-CMD-check/badge.svg)](https://github.com/robjhyndman/tsfeatures/actions)
<!-- badges: end -->

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
#> # A tibble: 4 × 20
#>   frequency nperi…¹ seaso…² trend   spike linea…³ curva…⁴ e_acf1 e_acf10 entropy
#>       <dbl>   <dbl>   <dbl> <dbl>   <dbl>   <dbl>   <dbl>  <dbl>   <dbl>   <dbl>
#> 1         1       0       1 0.125 2.10e-5    3.58    1.11  0.793   2.21    0.702
#> 2         1       0       1 0.985 3.01e-8    4.45    1.10  0.774   0.983   0.461
#> 3        12       1      12 0.991 1.46e-8   11.0     1.09  0.509   0.930   0.296
#> 4        12       1      12 0.802 9.15e-7   -2.12    2.85  0.258   0.341   0.548
#> # … with 10 more variables: x_acf1 <dbl>, x_acf10 <dbl>, diff1_acf1 <dbl>,
#> #   diff1_acf10 <dbl>, diff2_acf1 <dbl>, diff2_acf10 <dbl>,
#> #   seasonal_strength <dbl>, peak <dbl>, trough <dbl>, seas_acf1 <dbl>, and
#> #   abbreviated variable names ¹​nperiods, ²​seasonal_period, ³​linearity,
#> #   ⁴​curvature
```

## License

This package is free and open source software, licensed under GPL-3.
