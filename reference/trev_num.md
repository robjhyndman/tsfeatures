# Normalized nonlinear autocorrelation, the numerator of the trev function of a time series from software package `hctsa`

Calculates the numerator of the trev function, a normalized nonlinear
autocorrelation, The time lag is set to 1.

## Usage

``` r
trev_num(y)
```

## Arguments

- y:

  the input time series

## Value

the numerator of the trev function of a time series

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
trev_num(WWWusage)
#> [1] 109.1515
```
