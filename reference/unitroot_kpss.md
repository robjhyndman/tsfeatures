# Unit Root Test Statistics

`unitroot_kpss` computes the statistic for the Kwiatkowski et al. unit
root test using the default settings for the
[`ur.kpss`](https://rdrr.io/pkg/urca/man/ur.kpss.html) function.
`unitroot_pp` computes the statistic for the Phillips-Perron unit root
test using the default settings for the
[`ur.pp`](https://rdrr.io/pkg/urca/man/ur.pp.html) function.

## Usage

``` r
unitroot_kpss(x, ...)

unitroot_pp(x, ...)
```

## Arguments

- x:

  a univariate time series.

- ...:

  Other arguments are passed to the
  [`ur.kpss`](https://rdrr.io/pkg/urca/man/ur.kpss.html) or
  [`ur.kpss`](https://rdrr.io/pkg/urca/man/ur.kpss.html) functions.

## Value

A numeric value

## Author

Pablo Montero-Manso
