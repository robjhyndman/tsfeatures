# Time series features based on tiled windows

Computes feature of a time series based on tiled (non-overlapping)
windows. Means or variances are produced for all tiled windows. Then
stability is the variance of the means, while lumpiness is the variance
of the variances.

## Usage

``` r
lumpiness(x, width = ifelse(frequency(x) > 1, frequency(x), 10))

stability(x, width = ifelse(frequency(x) > 1, frequency(x), 10))
```

## Arguments

- x:

  a univariate time series

- width:

  size of sliding window

## Value

A numeric vector of length 2 containing a measure of lumpiness and a
measure of stability.

## Author

Earo Wang and Rob J Hyndman
