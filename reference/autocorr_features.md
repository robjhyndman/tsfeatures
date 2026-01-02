# The autocorrelation feature set from software package `hctsa`

Calculate the features that grouped as autocorrelation set, which have
been used in CompEngine database, using method introduced in package
`hctsa`.

## Usage

``` r
autocorr_features(x)
```

## Arguments

- x:

  the input time series

## Value

a vector with autocorrelation features

## Details

Features in this set are `embed2_incircle_1`, `embed2_incircle_2`,
`ac_9`, `firstmin_ac`, `trev_num`, `motiftwo_entro3`, and
`walker_propcross`.

## References

B.D. Fulcher and N.S. Jones. hctsa: A computational framework for
automated time-series phenotyping using massive feature extraction. Cell
Systems 5, 527 (2017).

B.D. Fulcher, M.A. Little, N.S. Jones Highly comparative time-series
analysis: the empirical structure of time series and their methods. J.
Roy. Soc. Interface 10, 83 (2013).

## See also

[`embed2_incircle`](http://pkg.robjhyndman.com/tsfeatures/reference/embed2_incircle.md)

[`ac_9`](http://pkg.robjhyndman.com/tsfeatures/reference/ac_9.md)

[`firstmin_ac`](http://pkg.robjhyndman.com/tsfeatures/reference/firstmin_ac.md)

[`trev_num`](http://pkg.robjhyndman.com/tsfeatures/reference/trev_num.md)

[`motiftwo_entro3`](http://pkg.robjhyndman.com/tsfeatures/reference/motiftwo_entro3.md)

[`walker_propcross`](http://pkg.robjhyndman.com/tsfeatures/reference/walker_propcross.md)

## Author

Yangzhuoran Yang
