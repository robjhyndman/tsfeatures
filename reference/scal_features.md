# The scaling feature set from software package `hctsa`

Calculate the features that grouped as scaling set, which have been used
in CompEngine database, using method introduced in package `hctsa`.

## Usage

``` r
scal_features(x)
```

## Arguments

- x:

  the input time series

## Value

a vector with scaling features

## Details

Feature in this set is `fluctanal_prop_r1`.

## References

B.D. Fulcher and N.S. Jones. hctsa: A computational framework for
automated time-series phenotyping using massive feature extraction. Cell
Systems 5, 527 (2017).

B.D. Fulcher, M.A. Little, N.S. Jones Highly comparative time-series
analysis: the empirical structure of time series and their methods. J.
Roy. Soc. Interface 10, 83 (2013).

## See also

[`fluctanal_prop_r1`](http://pkg.robjhyndman.com/tsfeatures/reference/fluctanal_prop_r1.md)

## Author

Yangzhuoran Yang
