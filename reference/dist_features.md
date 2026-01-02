# The distribution feature set from software package `hctsa`

Calculate the features that grouped as distribution set, which have been
used in CompEngine database, using method introduced in package `hctsa`.

## Usage

``` r
dist_features(x)
```

## Arguments

- x:

  the input time series

## Value

a vector with distribution features

## Details

Features in this set are `histogram_mode_10` and `outlierinclude_mdrmd`.

## References

B.D. Fulcher and N.S. Jones. hctsa: A computational framework for
automated time-series phenotyping using massive feature extraction. Cell
Systems 5, 527 (2017).

B.D. Fulcher, M.A. Little, N.S. Jones Highly comparative time-series
analysis: the empirical structure of time series and their methods. J.
Roy. Soc. Interface 10, 83 (2013).

## See also

[`histogram_mode`](http://pkg.robjhyndman.com/tsfeatures/reference/histogram_mode.md)

[`outlierinclude_mdrmd`](http://pkg.robjhyndman.com/tsfeatures/reference/outlierinclude_mdrmd.md)

## Author

Yangzhuoran Yang
