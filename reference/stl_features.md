# Strength of trend and seasonality of a time series

Computes various measures of trend and seasonality of a time series
based on an STL decomposition. The number of seasonal periods, and the
length of the seasonal periods are returned. Also, the strength of
seasonality corresponding to each period is estimated. The
[`mstl`](https://pkg.robjhyndman.com/forecast/reference/mstl.html)
function is used to do the decomposition.

## Usage

``` r
stl_features(x, ...)
```

## Arguments

- x:

  a univariate time series.

- ...:

  Other arguments are passed to
  [`mstl`](https://pkg.robjhyndman.com/forecast/reference/mstl.html).

## Value

A vector of numeric values.

## Author

Rob J Hyndman
