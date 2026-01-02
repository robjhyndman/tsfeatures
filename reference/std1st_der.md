# Standard deviation of the first derivative of the time series from software package `hctsa`

Modified from `SY_StdNthDer` in `hctsa`. Based on an idea by Vladimir
Vassilevsky.

## Usage

``` r
std1st_der(y)
```

## Arguments

- y:

  the input time series. Missing values will be removed.

## Value

Standard deviation of the first derivative of the time series.

## References

cf. http://www.mathworks.de/matlabcentral/newsreader/view_thread/136539

B.D. Fulcher and N.S. Jones. hctsa: A computational framework for
automated time-series phenotyping using massive feature extraction. Cell
Systems 5, 527 (2017).

B.D. Fulcher, M.A. Little, N.S. Jones Highly comparative time-series
analysis: the empirical structure of time series and their methods. J.
Roy. Soc. Interface 10, 83 (2013).

## Author

Yangzhuoran Yang
