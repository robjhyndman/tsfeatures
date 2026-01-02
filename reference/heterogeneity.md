# Heterogeneity coefficients

Computes various measures of heterogeneity of a time series. First the
series is pre-whitened using an AR model to give a new series y. We fit
a GARCH(1,1) model to y and obtain the residuals, e. Then the four
measures of heterogeneity are: (1) the sum of squares of the first 12
autocorrelations of \\y^2\\; (2) the sum of squares of the first 12
autocorrelations of \\e^2\\; (3) the \\R^2\\ value of an AR model
applied to \\y^2\\; (4) the \\R^2\\ value of an AR model applied to
\\e^2\\. The statistics obtained from \\y^2\\ are the ARCH effects,
while those from \\e^2\\ are the GARCH effects.

## Usage

``` r
heterogeneity(x)
```

## Arguments

- x:

  a univariate time series

## Value

A vector of numeric values.

## Author

Yanfei Kang and Rob J Hyndman
