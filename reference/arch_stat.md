# ARCH LM Statistic

Computes a statistic based on the Lagrange Multiplier (LM) test of Engle
(1982) for autoregressive conditional heteroscedasticity (ARCH). The
statistic returned is the \\R^2\\ value of an autoregressive model of
order `lags` applied to \\x^2\\.

## Usage

``` r
arch_stat(x, lags = 12, demean = TRUE)
```

## Arguments

- x:

  a univariate time series

- lags:

  Number of lags to use in the test

- demean:

  Should data have mean removed before test applied?

## Value

A numeric value.

## Author

Yanfei Kang
