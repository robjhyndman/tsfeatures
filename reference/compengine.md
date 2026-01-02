# CompEngine feature set

Calculate the features that have been used in CompEngine database, using
method introduced in package `hctsa`.

## Usage

``` r
compengine(x)
```

## Arguments

- x:

  the input time series

## Value

a vector with CompEngine features

## Details

The features involved can be grouped as `autocorrelation`, `prediction`,
`stationarity`, `distribution`, and `scaling`.

## References

B.D. Fulcher and N.S. Jones. hctsa: A computational framework for
automated time-series phenotyping using massive feature extraction. Cell
Systems 5, 527 (2017).

B.D. Fulcher, M.A. Little, N.S. Jones Highly comparative time-series
analysis: the empirical structure of time series and their methods. J.
Roy. Soc. Interface 10, 83 (2013).

## See also

[`autocorr_features`](http://pkg.robjhyndman.com/tsfeatures/reference/autocorr_features.md)

[`pred_features`](http://pkg.robjhyndman.com/tsfeatures/reference/pred_features.md)

[`station_features`](http://pkg.robjhyndman.com/tsfeatures/reference/station_features.md)

[`dist_features`](http://pkg.robjhyndman.com/tsfeatures/reference/dist_features.md)

[`scal_features`](http://pkg.robjhyndman.com/tsfeatures/reference/scal_features.md)

## Author

Yangzhuoran Yang
