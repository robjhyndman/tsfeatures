# Implements fluctuation analysis from software package `hctsa`

Fits a polynomial of order 1 and then returns the range. The order of
fluctuations is 2, corresponding to root mean square fluctuations.

## Usage

``` r
fluctanal_prop_r1(x)
```

## Arguments

- x:

  the input time series (or any vector)

## References

B.D. Fulcher and N.S. Jones. hctsa: A computational framework for
automated time-series phenotyping using massive feature extraction. Cell
Systems 5, 527 (2017).

B.D. Fulcher, M.A. Little, N.S. Jones Highly comparative time-series
analysis: the empirical structure of time series and their methods. J.
Roy. Soc. Interface 10, 83 (2013).

## Author

Yangzhuoran Yang
