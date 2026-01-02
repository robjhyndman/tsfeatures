# Partial autocorrelation-based features

Computes various measures based on partial autocorrelation coefficients
of the original series, first-differenced series and second-differenced
series

## Usage

``` r
pacf_features(x)
```

## Arguments

- x:

  a univariate time series

## Value

A vector of 3 values: Sum of squared of first 5 partial autocorrelation
coefficients of the original series, first differenced series and
twice-differenced series. For seasonal data, the partial autocorrelation
coefficient at the first seasonal lag is also returned.

## Author

Thiyanga Talagala
