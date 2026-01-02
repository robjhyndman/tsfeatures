# The prediction feature set from software package `hctsa`

Calculate the features that grouped as prediction set, which have been
used in CompEngine database, using method introduced in package `hctsa`.

## Usage

``` r
pred_features(x)
```

## Arguments

- x:

  the input time series

## Value

a vector with prediction features

## Details

Features in this set are `localsimple_mean1`, `localsimple_lfitac`, and
`sampen_first`.

## References

B.D. Fulcher and N.S. Jones. hctsa: A computational framework for
automated time-series phenotyping using massive feature extraction. Cell
Systems 5, 527 (2017).

B.D. Fulcher, M.A. Little, N.S. Jones Highly comparative time-series
analysis: the empirical structure of time series and their methods. J.
Roy. Soc. Interface 10, 83 (2013).

## See also

[`localsimple_taures`](http://pkg.robjhyndman.com/tsfeatures/reference/localsimple_taures.md)

[`sampen_first`](http://pkg.robjhyndman.com/tsfeatures/reference/sampen_first.md)

## Author

Yangzhuoran Yang
