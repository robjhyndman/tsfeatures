# Mode of a data vector from software package `hctsa`

Measures the mode of the data vector using histograms with a given
number of bins as suggestion. The value calculated is different from
`hctsa` and `CompEngine` as the histogram edges are calculated
differently.

## Usage

``` r
histogram_mode(y, numBins = 10)
```

## Arguments

- y:

  the input data vector

- numBins:

  the number of bins to use in the histogram.

## Value

the mode

## References

B.D. Fulcher and N.S. Jones. hctsa: A computational framework for
automated time-series phenotyping using massive feature extraction. Cell
Systems 5, 527 (2017).

B.D. Fulcher, M.A. Little, N.S. Jones Highly comparative time-series
analysis: the empirical structure of time series and their methods. J.
Roy. Soc. Interface 10, 83 (2013).

## Author

Yangzhuoran Yang
