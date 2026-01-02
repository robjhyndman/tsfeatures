# Second Sample Entropy from software package `hctsa`

Modified from the Ben Fulcher version of original code sampenc.m from
http://physionet.org/physiotools/sampen/
http://www.physionet.org/physiotools/sampen/matlab/1.1/sampenc.m Code by
DK Lake (dlake@virginia.edu), JR Moorman and Cao Hanqing.

## Usage

``` r
sampenc(y, M = 6, r = 0.3)
```

## Arguments

- y:

  the input time series

- M:

  embedding dimension

- r:

  threshold

## References

cf. "Physiological time-series analysis using approximate entropy and
sample entropy", J. S. Richman and J. R. Moorman, Am. J. Physiol. Heart
Circ. Physiol., 278(6) H2039 (2000)

B.D. Fulcher and N.S. Jones. hctsa: A computational framework for
automated time-series phenotyping using massive feature extraction. Cell
Systems 5, 527 (2017).

B.D. Fulcher, M.A. Little, N.S. Jones Highly comparative time-series
analysis: the empirical structure of time series and their methods. J.
Roy. Soc. Interface 10, 83 (2013).

## Author

Yangzhuoran Yang
