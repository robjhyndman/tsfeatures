# Parameter estimates of Holt's linear trend method

Estimate the smoothing parameter for the level-alpha and the smoothing
parameter for the trend-beta. `hw_parameters` considers additive
seasonal trend: ets(A,A,A) model.

## Usage

``` r
holt_parameters(x)

hw_parameters(x)
```

## Arguments

- x:

  a univariate time series

## Value

`holt_parameters` produces a vector of 2 values: alpha, beta.

`hw_parameters` produces a vector of 3 values: alpha, beta and gamma.

## Author

Thiyanga Talagala, Pablo Montero-Manso
