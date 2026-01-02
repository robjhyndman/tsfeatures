# Time of first minimum in the autocorrelation function from software package `hctsa`

Time of first minimum in the autocorrelation function from software
package `hctsa`

## Usage

``` r
firstmin_ac(
  x,
  acfv = stats::acf(x, lag.max = N - 1, plot = FALSE, na.action = na.pass)
)
```

## Arguments

- x:

  the input time series

- acfv:

  vector of autocorrelation, if exist, used to avoid repeated
  computation.

## Value

The lag of the first minimum

## References

B.D. Fulcher and N.S. Jones. hctsa: A computational framework for
automated time-series phenotyping using massive feature extraction. Cell
Systems 5, 527 (2017).

B.D. Fulcher, M.A. Little, N.S. Jones Highly comparative time-series
analysis: the empirical structure of time series and their methods. J.
Roy. Soc. Interface 10, 83 (2013).

## Author

Yangzhuoran Yang

## Examples

``` r
firstmin_ac(WWWusage)
#> [1] 21
```
