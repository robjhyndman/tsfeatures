# Time series feature matrix

`tsfeatures` computes a matrix of time series features from a list of
time series

## Usage

``` r
tsfeatures(
  tslist,
  features = c("frequency", "stl_features", "entropy", "acf_features"),
  scale = TRUE,
  trim = FALSE,
  trim_amount = 0.1,
  parallel = FALSE,
  multiprocess = future::multisession,
  na.action = na.pass,
  ...
)
```

## Arguments

- tslist:

  a list of univariate time series, each of class `ts` or a numeric
  vector. Alternatively, an object of class `mts` may be used.

- features:

  a vector of function names which return numeric vectors of features.
  All features returned by these functions must be named if they return
  more than one feature. Existing functions from installed packages may
  be used, but the package must be loaded first. Functions must return a
  result for all time series, even if it is just NA.

- scale:

  if `TRUE`, time series are scaled to mean 0 and sd 1 before features
  are computed.

- trim:

  if `TRUE`, time series are trimmed by `trim_amount` before features
  are computed. Values larger than `trim_amount` in absolute value are
  set to `NA`.

- trim_amount:

  Default level of trimming if `trim==TRUE`.

- parallel:

  If TRUE, multiple cores (or multiple sessions) will be used. This only
  speeds things up when there are a large number of time series.

- multiprocess:

  The function from the `future` package to use for parallel processing.
  Either
  [`multisession`](https://future.futureverse.org/reference/multisession.html)
  or
  [`multicore`](https://future.futureverse.org/reference/multicore.html).
  The latter is preferred for Linux and MacOS.

- na.action:

  A function to handle missing values. Use `na.interp` to estimate
  missing values.

- ...:

  Other arguments get passed to the feature functions.

## Value

A feature matrix (in the form of a tibble) with each row corresponding
to one time series from tslist, and each column being a feature.

## Author

Rob J Hyndman

## Examples

``` r
mylist <- list(sunspot.year, WWWusage, AirPassengers, USAccDeaths)
tsfeatures(mylist)
#> # A tibble: 4 × 20
#>   frequency nperiods seasonal_period trend      spike linearity curvature e_acf1
#>       <dbl>    <dbl>           <dbl> <dbl>      <dbl>     <dbl>     <dbl>  <dbl>
#> 1         1        0               1 0.125    2.10e-5      3.58      1.11  0.793
#> 2         1        0               1 0.985    3.01e-8      4.45      1.10  0.774
#> 3        12        1              12 0.991    1.46e-8     11.0       1.09  0.509
#> 4        12        1              12 0.802    9.15e-7     -2.12      2.85  0.258
#> # ℹ 12 more variables: e_acf10 <dbl>, entropy <dbl>, x_acf1 <dbl>,
#> #   x_acf10 <dbl>, diff1_acf1 <dbl>, diff1_acf10 <dbl>, diff2_acf1 <dbl>,
#> #   diff2_acf10 <dbl>, seasonal_strength <dbl>, peak <dbl>, trough <dbl>,
#> #   seas_acf1 <dbl>
```
