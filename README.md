
<!-- README.md is generated from README.Rmd. Please edit that file -->
tsfeatures
==========

[![Pending Pull-Requests](http://githubbadges.herokuapp.com/robjhyndman/tsfeatures/pulls.svg?style=flat)](https://github.com/robjhyndman/tsfeatures/pulls)

The R package *tsfeatures* provides methods for extracting various features from time series data.

Installation
------------

The **stable** version on R CRAN is coming soon.

You can install the **development** version from [Github](https://github.com/robjhyndman/tsfeatures) with:

``` r
# install.packages("devtools")
devtools::install_github("robjhyndman/tsfeatures")
```

Usage (still in development)
----------------------------

``` r
# Compute features for M3 paper: Kang, Hyndman & Smith-Miles (IJF 2017)

library(tsfeatures)
library(tidyverse)

M3data <- purrr::map(Mcomp::M3, function(x)x$x)
Lambda <- function(x){
  suppressWarnings(forecast::BoxCox.lambda(x, lower=0, upper=1))
}
khs <- tsfeatures(M3data, c("entropy","stl_features","frequency","Lambda")) %>% 
  as_tibble %>%
  select(frequency, entropy, trend, season, acfremainder, Lambda) %>%
  replace_na(list(season=1)) %>%
  rename(
    Period = frequency,
    Entropy = entropy,
    Trend = trend,
    Season = season,
    ACF1 = acfremainder) %>%
  mutate(Period = as.factor(Period))
```

``` r
# Fig 1 of paper
GGally::ggpairs(khs)
```

![](README-ijf2017graphs-1.png)

``` r

# 2-d Feature space (Top of Fig 2)
prcomp(khs[,-1])$x %>%
  as_tibble() %>%
  bind_cols(Period=khs$Period) %>%
  ggplot(aes(x=PC1, y=PC2)) +
    geom_point(aes(col=Period))
```

![](README-ijf2017graphs-2.png)

``` r
library(tsfeatures)
library(tidyverse)
library(forecast)

M3data <- purrr::map(Mcomp::M3, function(x)x$x)

# Compute all features in MARs paper
nsdiffs <- function(x){
  c(nsdiffs=ifelse(frequency(x)==1L, 1L, forecast::nsdiffs(x)))
}
features1 <- c(
  "entropy",
  "stl_features",
  "acf1",
  "max_kl_shift",
  "lumpiness",
  "ndiffs",
  "nonlinearity",
  "frequency",
  "nsdiffs",
  "heterogeneity"
)
features2 <- c(  
  "max_level_shift",
  "max_var_shift"
)

yk <- cbind(
        tsfeatures(M3data, features1),
        tsfeatures(M3data, features2, trim=TRUE)) %>% 
  as_tibble %>%
  select(entropy, trend, linearity, curvature, acf1,
    acfremainder, max_kl_shift, max_level_shift, max_var_shift, 
    ndiffs, Nonlinearity, frequency, season, nsdiffs,
    ARCHtest.p, GARCHtest.p, Boxtest.p, GARCHBoxtest.p, Hetero)
# What happened to var change on remainder?
```

``` r
# Compute all features in Yahoo anomaly paper

yahoo_DT <- c("dat0", "dat1", "dat2", "dat3", "dat4", "dat5")
data(list = yahoo_DT, package = "anomalous")

r <- lapply(yahoo_DT, function(name_DT){
  DT <- get(name_DT)
  return(lapply(seq(ncol(DT)), function(j){
    DT[, j]
  }))
})
list_real <- c(r[[1]], r[[2]], r[[3]], r[[4]])
list_simulated <- c(r[[5]], r[[6]])

# real world yahoo series
myfeatures <- c("entropy", "af1", "stl_features")
mat <- tsfeatures(list_real, myfeatures)

# simulated yahoo series
mat2 <- tsfeatures(list_simulated, myfeatures)
```

License
-------

This package is free and open source software, licensed under GPL-3.
