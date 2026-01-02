# Yahoo server metrics

Yahoo server metrics

## Usage

``` r
yahoo_data(...)
```

## Arguments

- ...:

  Additional arguments passed to `download.file`

  Downloads and returns aggregated and anonymized datasets from Yahoo
  representing server metrics of Yahoo services.

## Value

A matrix of time series with 1437 rows of hourly data, and 1748 columns
representing different servers.

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
plot(yahoo[,1:10])

plot(yahoo[,1:44], plot.type='single', col=1:44)

```
