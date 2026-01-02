# Points inside a given circular boundary in a 2-d embedding space from software package `hctsa`

The time lag is set to the first zero crossing of the autocorrelation
function.

## Usage

``` r
embed2_incircle(
  y,
  boundary = NULL,
  acfv = stats::acf(y, length(y) - 1, plot = FALSE, na.action = na.pass)
)
```

## Arguments

- y:

  the input time series

- boundary:

  the given circular boundary, setting to 1 or 2 in CompEngine. Default
  to 1.

- acfv:

  vector of autocorrelation, if exist, used to avoid repeated
  computation.

## Value

the proportion of points inside a given circular boundary

## References

B.D. Fulcher and N.S. Jones. hctsa: A computational framework for
automated time-series phenotyping using massive feature extraction. Cell
Systems 5, 527 (2017).

B.D. Fulcher, M.A. Little, N.S. Jones Highly comparative time-series
analysis: the empirical structure of time series and their methods. J.
Roy. Soc. Interface 10, 83 (2013).

## Author

Yangzhuoran Yang
