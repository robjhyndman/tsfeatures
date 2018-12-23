
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tsfeatures

[![Pending
Pull-Requests](http://githubbadges.herokuapp.com/robjhyndman/tsfeatures/pulls.svg?style=flat)](https://github.com/robjhyndman/tsfeatures/pulls)

The R package *tsfeatures* provides methods for extracting various
features from time series data.

## Installation

The **stable** version on R CRAN is coming soon.

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
