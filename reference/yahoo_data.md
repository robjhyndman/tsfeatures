# Yahoo server metrics

Downloads and returns aggregated and anonymized datasets from Yahoo
representing server metrics of Yahoo services. If the download fails
(e.g., no internet connection or the server is unavailable), a message
is shown and `NULL` is returned invisibly.

## Usage

``` r
yahoo_data(...)
```

## Arguments

- ...:

  Additional arguments passed to `download.file`

## Value

A matrix of time series with 1437 rows of hourly data, and 1748 columns
representing different servers, or `NULL` if the data could not be
downloaded.

## References

Hyndman, R.J., Wang, E., Laptev, N. (2015) Large-scale unusual time
series detection. In: *Proceedings of the IEEE International Conference
on Data Mining*. Atlantic City, NJ, USA. 14–17 November 2015.
<https://robjhyndman.com/publications/icdm2015/>

## Author

Rob Hyndman, Earo Wang, Nikolay Laptev, Mitchell O'Hara-Wild

## Examples

``` r
yahoo <- yahoo_data()
if (!is.null(yahoo)) {
  plot(yahoo[, 1:10])
  plot(yahoo[, 1:44], plot.type = "single", col = 1:44)
}

```
