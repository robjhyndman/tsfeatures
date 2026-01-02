# Convert mts object to list of time series

An mts object contains a multivariate time series in a matrix, with time
on rows. This is converted into a list of univariate time series.

## Usage

``` r
# S3 method for class 'mts'
as.list(x, ...)
```

## Arguments

- x:

  multivariate time series of class mts.

- ...:

  other arguments are ignored.

## Value

A list of ts objects.

## Author

Rob J Hyndman
