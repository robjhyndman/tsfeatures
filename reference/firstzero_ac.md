# The first zero crossing of the autocorrelation function from software package `hctsa`

Search up to a maximum of the length of the time series

## Usage

``` r
firstzero_ac(y, acfv = stats::acf(y, N - 1, plot = FALSE, na.action = na.pass))
```

## Arguments

- y:

  the input time series

- acfv:

  vector of autocorrelation, if exist, used to avoid repeated
  computation.

## Value

The first zero crossing of the autocorrelation function

## References

B.D. Fulcher and N.S. Jones. hctsa: A computational framework for
automated time-series phenotyping using massive feature extraction. Cell
Systems 5, 527 (2017).

B.D. Fulcher, M.A. Little, N.S. Jones Highly comparative time-series
analysis: the empirical structure of time series and their methods. J.
Roy. Soc. Interface 10, 83 (2013).

## Author

Yangzhuoran Yang
