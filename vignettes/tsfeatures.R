## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = FALSE,
  fig.align = "center"
)
library(tsfeatures)

## ----gh-installation, eval = FALSE---------------------------------------
#  # install.packages("devtools")
#  devtools::install_github("robjhyndman/tsfeatures")

## ------------------------------------------------------------------------
mylist <- list(sunspot.year, WWWusage, AirPassengers, USAccDeaths)
tsfeatures(mylist)

## ------------------------------------------------------------------------
# Function from outside of tsfeatures package being used
is.monthly <- function(x){
  frequency(x) == 12
}
tsfeatures(mylist, features = "is.monthly")

## ------------------------------------------------------------------------
mylist <- list(sunspot.year, WWWusage, AirPassengers, USAccDeaths)
fmat <- tsfeatures(mylist)
z <- embedding(dist(fmat), perplexity=1)
plot(z)

## ------------------------------------------------------------------------
entropy(AirPassengers)

## ------------------------------------------------------------------------
stability(AirPassengers)
lumpiness(AirPassengers)

## ------------------------------------------------------------------------
max_level_shift(AirPassengers)
max_var_shift(AirPassengers)
max_kl_shift(AirPassengers)

## ------------------------------------------------------------------------
crossing_points(AirPassengers)

## ------------------------------------------------------------------------
flat_spots(AirPassengers)

## ------------------------------------------------------------------------
hurst(AirPassengers)

## ------------------------------------------------------------------------
unitroot_kpss(AirPassengers)
unitroot_pp(AirPassengers)

## ------------------------------------------------------------------------
stl_features(AirPassengers)

## ------------------------------------------------------------------------
acf_features(AirPassengers)

## ------------------------------------------------------------------------
pacf_features(AirPassengers)

## ------------------------------------------------------------------------
holt_parameters(AirPassengers)
hw_parameters(AirPassengers)

## ------------------------------------------------------------------------
heterogeneity(AirPassengers)

## ------------------------------------------------------------------------
nonlinearity(AirPassengers)

## ------------------------------------------------------------------------
arch_stat(AirPassengers)

## ------------------------------------------------------------------------
comp <- compengine(AirPassengers)
knitr::kable(comp)

## ------------------------------------------------------------------------
embed2_incircle(AirPassengers, boundary = 1e5)

## ------------------------------------------------------------------------
ac_9(AirPassengers)

## ------------------------------------------------------------------------
firstmin_ac(AirPassengers)

## ------------------------------------------------------------------------
firstzero_ac(AirPassengers)

## ------------------------------------------------------------------------
trev_num(AirPassengers)

## ------------------------------------------------------------------------
motiftwo_entro3(AirPassengers)

## ------------------------------------------------------------------------
str(binarize_mean(AirPassengers))

## ------------------------------------------------------------------------
walker_propcross(AirPassengers)

## ------------------------------------------------------------------------
localsimple_taures(AirPassengers)

## ------------------------------------------------------------------------
sampen_first(AirPassengers)
sampenc(AirPassengers, M = 5, r = 0.3)

## ------------------------------------------------------------------------
std1st_der(AirPassengers)

## ------------------------------------------------------------------------
spreadrandomlocal_meantaul(AirPassengers)

## ------------------------------------------------------------------------
histogram_mode(AirPassengers, numBins = 10)

## ------------------------------------------------------------------------
outlierinclude_mdrmd(AirPassengers)

## ------------------------------------------------------------------------
fluctanal_prop_r1(AirPassengers)

## ----yahoo---------------------------------------------------------------
library(tsfeatures)
library(dplyr)
library(anomalous)

yahoo <- cbind(dat0, dat1, dat2, dat3)
hwl <- bind_cols(
         tsfeatures(yahoo,
           c("acf_features","entropy","lumpiness",
             "flat_spots","crossing_points")),
         tsfeatures(yahoo,"stl_features", s.window='periodic', robust=TRUE),
         tsfeatures(yahoo, "max_kl_shift", width=48),
         tsfeatures(yahoo,
           c("mean","var"), scale=FALSE, na.rm=TRUE),
         tsfeatures(yahoo,
           c("max_level_shift","max_var_shift"), trim=TRUE)) %>%
  select(mean, var, x_acf1, trend, linearity, curvature,
         seasonal_strength, peak, trough,
         entropy, lumpiness, spike, max_level_shift, max_var_shift, flat_spots,
         crossing_points, max_kl_shift, time_kl_shift)

## ----yahoo2--------------------------------------------------------------
# 2-d Feature space
library(ggplot2)
prcomp(na.omit(hwl), scale=TRUE)$x %>%
  as_tibble() %>%
  ggplot(aes(x=PC1, y=PC2)) +
    geom_point()

## ----ijf2017-------------------------------------------------------------
library(tsfeatures)
library(dplyr)
library(tidyr)
library(forecast)

M3data <- purrr::map(Mcomp::M3,
  function(x){
      tspx <- tsp(x$x)
      ts(c(x$x,x$xx), start=tspx[1], frequency=tspx[3])
  })
khs_stl <- function(x,...)
{
  lambda <- BoxCox.lambda(x, lower=0, upper=1, method='loglik')
  y <- BoxCox(x, lambda)
  c(stl_features(y,s.window='periodic', robust=TRUE, ...), lambda=lambda)
}
khs <- bind_cols(
  tsfeatures(M3data, c("frequency", "entropy")),
  tsfeatures(M3data, "khs_stl", scale=FALSE)) %>%
  select(frequency, entropy, trend, seasonal_strength, e_acf1, lambda) %>%
  replace_na(list(seasonal_strength=0)) %>%
  rename(
    Frequency = frequency,
    Entropy = entropy,
    Trend = trend,
    Season = seasonal_strength,
    ACF1 = e_acf1,
    Lambda = lambda) %>%
  mutate(Period = as.factor(Frequency))

## ----ijf2017graphs, dependson="ijf2017"----------------------------------
# Fig 1 of paper
khs %>%
  select(Period, Entropy, Trend, Season, ACF1, Lambda) %>%
  GGally::ggpairs()

# 2-d Feature space (Top of Fig 2)
prcomp(select(khs, -Period), scale=TRUE)$x %>%
  as_tibble() %>%
  bind_cols(Period=khs$Period) %>%
  ggplot(aes(x=PC1, y=PC2)) +
    geom_point(aes(col=Period))

