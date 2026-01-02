# Simulates a hypothetical walker moving through the time domain from software package `hctsa`

The hypothetical particle (or 'walker') moves in response to values of
the time series at each point. The walker narrows the gap between its
value and that of the time series by 10%.

## Usage

``` r
walker_propcross(y)
```

## Arguments

- y:

  the input time series

## Value

fraction of time series length that walker crosses time series

## References

B.D. Fulcher and N.S. Jones. hctsa: A computational framework for
automated time-series phenotyping using massive feature extraction. Cell
Systems 5, 527 (2017).

B.D. Fulcher, M.A. Little, N.S. Jones Highly comparative time-series
analysis: the empirical structure of time series and their methods. J.
Roy. Soc. Interface 10, 83 (2013).

## Author

Yangzhuoran Yang
