# Time series features based on sliding windows

Computes feature of a time series based on sliding (overlapping)
windows. `max_level_shift` finds the largest mean shift between two
consecutive windows. `max_var_shift` finds the largest var shift between
two consecutive windows. `max_kl_shift` finds the largest shift in
Kulback-Leibler divergence between two consecutive windows.

## Usage

``` r
max_level_shift(x, width = ifelse(frequency(x) > 1, frequency(x), 10))

max_var_shift(x, width = ifelse(frequency(x) > 1, frequency(x), 10))

max_kl_shift(x, width = ifelse(frequency(x) > 1, frequency(x), 10))
```

## Arguments

- x:

  a univariate time series

- width:

  size of sliding window

## Value

A vector of 2 values: the size of the shift, and the time index of the
shift.

## Details

Computes the largest level shift and largest variance shift in sliding
mean calculations

## Author

Earo Wang and Rob J Hyndman
