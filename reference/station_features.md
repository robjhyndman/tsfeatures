# The stationarity feature set from software package `hctsa`

Calculate the features that grouped as stationarity set, which have been
used in CompEngine database, using method introduced in package `hctsa`.

## Usage

``` r
station_features(x)
```

## Arguments

- x:

  the input time series

## Value

a vector with stationarity features

## Details

Features in this set are `std1st_der`, `spreadrandomlocal_meantaul_50`,
and `spreadrandomlocal_meantaul_ac2`.

## References

B.D. Fulcher and N.S. Jones. hctsa: A computational framework for
automated time-series phenotyping using massive feature extraction. Cell
Systems 5, 527 (2017).

B.D. Fulcher, M.A. Little, N.S. Jones Highly comparative time-series
analysis: the empirical structure of time series and their methods. J.
Roy. Soc. Interface 10, 83 (2013).

## See also

[`std1st_der`](http://pkg.robjhyndman.com/tsfeatures/reference/std1st_der.md)

[`spreadrandomlocal_meantaul`](http://pkg.robjhyndman.com/tsfeatures/reference/spreadrandomlocal_meantaul.md)

## Author

Yangzhuoran Yang
