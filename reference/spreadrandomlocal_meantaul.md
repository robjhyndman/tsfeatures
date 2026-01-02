# Bootstrap-based stationarity measure from software package `hctsa`

100 time-series segments of length `l` are selected at random from the
time series and the mean of the first zero-crossings of the
autocorrelation function in each segment is calculated.

## Usage

``` r
spreadrandomlocal_meantaul(y, l = 50)
```

## Arguments

- y:

  the input time series

- l:

  the length of local time-series segments to analyse as a positive
  integer. Can also be a specified character string: "ac2": twice the
  first zero-crossing of the autocorrelation function

## Value

mean of the first zero-crossings of the autocorrelation function

## References

B.D. Fulcher and N.S. Jones. hctsa: A computational framework for
automated time-series phenotyping using massive feature extraction. Cell
Systems 5, 527 (2017).

B.D. Fulcher, M.A. Little, N.S. Jones Highly comparative time-series
analysis: the empirical structure of time series and their methods. J.
Roy. Soc. Interface 10, 83 (2013).

## Author

Yangzhuoran Yang
