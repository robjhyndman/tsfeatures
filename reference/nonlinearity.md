# Nonlinearity coefficient

Computes a nonlinearity statistic based on Lee, White & Granger's
nonlinearity test of a time series. The statistic is \\10X^2/T\\ where
\\X^2\\ is the Chi-squared statistic from Lee, White and Granger, and T
is the length of the time series. This takes large values when the
series is nonlinear, and values around 0 when the series is linear.

## Usage

``` r
nonlinearity(x)
```

## Arguments

- x:

  a univariate time series

## Value

A numeric value.

## References

Lee, T. H., White, H., & Granger, C. W. (1993). Testing for neglected
nonlinearity in time series models: A comparison of neural network
methods and alternative tests. *Journal of Econometrics*, 56(3),
269-290.

Teräsvirta, T., Lin, C.-F., & Granger, C. W. J. (1993). Power of the
neural network linearity test. *Journal of Time Series Analysis*, 14(2),
209–220.

## Author

Yanfei Kang and Rob J Hyndman

## Examples

``` r
nonlinearity(lynx)
#> nonlinearity 
#>    0.8959046 
```
