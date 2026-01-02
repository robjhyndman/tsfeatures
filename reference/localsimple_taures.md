# The first zero crossing of the autocorrelation function of the residuals from Simple local time-series forecasting from software package `hctsa`

Simple predictors using the past trainLength values of the time series
to predict its next value.

## Usage

``` r
localsimple_taures(y, forecastMeth = c("mean", "lfit"), trainLength = NULL)
```

## Arguments

- y:

  the input time series

- forecastMeth:

  the forecasting method, default to `mean`. `mean`: local mean
  prediction using the past trainLength time-series values. `lfit`:
  local linear prediction using the past trainLength time-series values.

- trainLength:

  the number of time-series values to use to forecast the next value.
  Default to 1 when using method `mean` and 3 when using method `lfit`.

## Value

The first zero crossing of the autocorrelation function of the residuals
