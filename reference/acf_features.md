# Autocorrelation-based features

Computes various measures based on autocorrelation coefficients of the
original series, first-differenced series and second-differenced series

## Usage

``` r
acf_features(x)
```

## Arguments

- x:

  a univariate time series

## Value

A vector of 6 values: first autocorrelation coefficient and sum of
squared of first ten autocorrelation coefficients of original series,
first-differenced series, and twice-differenced series. For seasonal
data, the autocorrelation coefficient at the first seasonal lag is also
returned.

## Author

Thiyanga Talagala
