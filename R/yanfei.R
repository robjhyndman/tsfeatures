#' Heterogeneity coefficients
#'
#' Computes various measures of heterogeneity of a time series. First the series
#' is pre-whitened using an AR model to give a new series y. We fit a GARCH(1,1)
#' model to y and obtain the residuals, e. Then the four measures of heterogeneity
#' are: 
#' (1) the sum of squares of the first 12 autocorrelations of \eqn{y^2}{y^2};
#' (2) the sum of squares of the first 12 autocorrelations of \eqn{e^2}{e^2};
#' (3) the \eqn{R^2}{R^2} value of an AR model applied to \eqn{y^2}{y^2};
#' (4) the \eqn{R^2}{R^2} value of an AR model applied to \eqn{e^2}{e^2}.
#' The statistics obtained from \eqn{y^2}{y^2} are the ARCH effects, while those
#' from \eqn{e^2}{e^2} are the GARCH effects.
#' @param x a univariate time series
#' @return A vector of numeric values.
#' @author Yanfei Kang and Rob J Hyndman
#' @export

heterogeneity <- function(x)
{
  # One possible issue when applied to the ETS/ARIMA comparison is that it will
  # be high for any type of heteroskedasticity, whereas ETS heteroskedasticity
  # is of a particular type, namely that the variation increases with the level
  # of the series. But the GARCH type hetero could be high when the variation
  # changes independently of the level of the series.


  # fit arima model to fit the mean dynamics.
  # x.freq <- frequency(x)
  # if (x.freq != 1){
  #   seasonal <- TRUE
  # }else{
  #   seasonal <- FALSE
  # }
  # x.arma <- auto.arima(x, max.p = 2, max.q = 2,
  #                      max.P = 2, max.Q = 2,
  #                      max.d = 2, max.D = 1,
  #                      seasonal = seasonal)
  #
  # # perform arch test
  # x.arma.archtest <- ArchTest(x.arma$residuals)
  # # x.arma.archtest <- ArchTest(x)
  #
  # # perform box test
  # x.arma.boxtest <- Box.test(x.arma$residuals^2, lag = 12, type = 'Ljung-Box')
  # x.arma.boxtest <- Box.test(x^2, lag = 12, type = 'Ljung-Box')

  # pre-whiten a series before Garch modeling
  x.whitened <- na.contiguous(ar(x)$resid)
  #x.whitened <- prewhiten.ts(x, AR.max = 24, plot = FALSE, verbose = FALSE)$prew_ar

  # perform arch and box test
  x.archtest <- arch_stat(x.whitened)
  LBstat <- sum(acf(x.whitened^2, lag.max=12L, plot=FALSE)$acf[-1L]^2)
  #x.boxtest <- Box.test(x.whitened^2, lag = 12, type = 'Ljung-Box')$p.value

  # fit garch model to capture the variance dynamics.
  garch.fit <- suppressWarnings(tseries::garch(x.whitened, trace=FALSE))
  # garch.fit <- garchFit(~ garch(1,1), data = x, trace = FALSE)

  # For a time series with strong ARCH/GARCH effects, sigma2 will have larger coefficient of variation?
  #sigma2 <- garch.fit$fitted.values[,1]
  #cor.sigma2.x <- suppressWarnings(cor(cbind(sigma2, x),
  #  use="pairwise.complete.obs")[1,2])
  #max.sigma2 <- max(sigma2, na.rm=TRUE)
  #min.sigma2 <- min(sigma2, na.rm=TRUE)
  #range.sigma2 <- max.sigma2 - min.sigma2
  #mean.sigma2 <- mean(sigma2, na.rm=TRUE)
  #std.sigma2 <- sd(sigma2, na.rm=TRUE)
  #cv.sigma2 <- std.sigma2/mean.sigma2

  # compare arch test before and after fitting garch
  garch.fit.std <- residuals(garch.fit)
  x.garch.archtest <- arch_stat(garch.fit.std)
  #archtest.p.diff <- x.garch.archtest - x.archtest

  # compare Box test of squared residuals before and after fitting garch
  #x.garch.boxtest <- Box.test(garch.fit.std^2, lag = 12, type = 'Ljung-Box')$p.value
  LBstat2 <- sum(acf(na.contiguous(garch.fit.std^2), lag.max=12L, plot=FALSE)$acf[-1L]^2)
  #boxtest.p.diff <- x.garch.boxtest - x.boxtest

  output <- c(
    ARCH_ACF = LBstat,
    GARCH_ACF = LBstat2,
    ARCH_R2 = unname(x.archtest),
    GARCH_R2 = unname(x.garch.archtest)
    #Hetero = max(archtest.p.diff, boxtest.p.diff)
    )
  output[is.na(output)] <- 1
  return(output)
}


#' Nonlinearity coefficent
#'
#' Computes a nonlinearity statistic based on Terasvirta's nonlinearity test of a time series.
#' The statistic is \eqn{10X^2/T}{10X^2/T} where \eqn{X^2}{X^2} is the Chi-squared statistic from
#' Terasvirta's test, and T is the length of the time series. This takes large values
#' when the series is nonlinear, and values around 0 when the series is linear.
#' @param x a univariate time series
#' @return A numeric value.
#' @examples
#' nonlinearity(lynx)
#' @author Yanfei Kang and Rob J Hyndman
#' @export

nonlinearity <- function(x)
{
  X2 <- tseries::terasvirta.test(as.ts(x),type = "Chisq")$stat
  c(nonlinearity = 10*unname(X2)/length(x))
}


#' ARCH LM Statistic
#'
#' Computes a statistic based on the Lagrange Multiplier (LM) test of Engle (1982) for
#' autoregressive conditional heteroscedasticity (ARCH). The statistic returned is
#' the \eqn{R^2}{R^2} value of an autoregressive model of order \code{lags} applied
#' to \eqn{x^2}{x^2}.
#' @param x a univariate time series
#' @param lags Number of lags to use in the test
#' @param demean Should data have mean removed before test applied?
#' @return A numeric value.
#' @author Yanfei Kang
#' @export

arch_stat <- function (x, lags = 12, demean=TRUE)
{
  if(length(x) <= 13)
    return(c(ARCH.LM=NA_real_))
  if(demean)
    x <- x - mean(x, na.rm=TRUE)
  mat <- embed(x^2, lags+1)
  fit <- try(lm(mat[, 1] ~ mat[, -1]), silent=TRUE)
  if("try-error" %in% class(fit))
    return(c(ARCH.LM=NA_real_))
  else
  {
    arch.lm <- summary(fit)
    S <- arch.lm$r.squared #* NROW(mat)
#    p2 <- 1 - pchisq(S, df=lags)
    return(c(ARCH.LM=S))
  }
}
