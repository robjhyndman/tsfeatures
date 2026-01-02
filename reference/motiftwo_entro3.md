# Local motifs in a binary symbolization of the time series from software package `hctsa`

Coarse-graining is performed. Time-series values above its mean are
given 1, and those below the mean are 0.

## Usage

``` r
motiftwo_entro3(y)
```

## Arguments

- y:

  the input time series

## Value

Entropy of words in the binary alphabet of length 3.

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
motiftwo_entro3(WWWusage)
#> [1] 1.200579
```
