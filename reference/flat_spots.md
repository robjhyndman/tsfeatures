# Longest flat spot

"Flat spots” are computed by dividing the sample space of a time series
into ten equal-sized intervals, and computing the maximum run length
within any single interval.

## Usage

``` r
flat_spots(x)
```

## Arguments

- x:

  a univariate time series

## Value

A numeric value.

## Author

Earo Wang and Rob J Hyndman
