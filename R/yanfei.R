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

heterogeneity <- function(x) {
  # One possible issue when applied to the ETS/ARIMA comparison is that it will
  # be high for any type of heteroskedasticity, whereas ETS heteroskedasticity
  # is of a particular type, namely that the variation increases with the level
  # of the series. But the GARCH type hetero could be high when the variation
  # changes independently of the level of the series.

  # pre-whiten a series before Garch modeling
  x.whitened <- na.contiguous(ar(x)$resid)

  # perform arch and box test
  x.archtest <- arch_stat(x.whitened)
  LBstat <- sum(acf(x.whitened^2, lag.max = 12L, plot = FALSE)$acf[-1L]^2)

  # fit garch model to capture the variance dynamics.
  garch.fit <- suppressWarnings(tseries::garch(x.whitened, trace = FALSE))

  # compare arch test before and after fitting garch
  garch.fit.std <- residuals(garch.fit)
  x.garch.archtest <- arch_stat(garch.fit.std)

  # compare Box test of squared residuals before and after fitting garch
  LBstat2 <- NA
  try(LBstat2 <- sum(acf(na.contiguous(garch.fit.std^2), lag.max = 12L, plot = FALSE)$acf[-1L]^2),
    silent = TRUE
  )
  output <- c(
    arch_acf = LBstat,
    garch_acf = LBstat2,
    arch_r2 = unname(x.archtest),
    garch_r2 = unname(x.garch.archtest)
  )
  # output[is.na(output)] <- 1
  return(output)
}

#' Nonlinearity coefficient
#'
#' Computes a nonlinearity statistic based on Teräsvirta's nonlinearity test of a time series.
#' The statistic is \eqn{10X^2/T}{10X^2/T} where \eqn{X^2}{X^2} is the Chi-squared statistic from
#' Teräsvirta's test, and T is the length of the time series. This takes large values
#' when the series is nonlinear, and values around 0 when the series is linear.
#' @param x a univariate time series
#' @return A numeric value.
#' @examples
#' nonlinearity(lynx)
#' @author Yanfei Kang and Rob J Hyndman
#' @export

nonlinearity <- function(x) {
  X2 <- tryCatch(tseries::terasvirta.test(as.ts(x), type = "Chisq")$stat,
                 error = function(e) NA)
  c(nonlinearity = 10 * unname(X2) / length(x))
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

arch_stat <- function(x, lags = 12, demean = TRUE) {
  if (length(x) <= lags+1) {
    return(c(ARCH.LM = NA_real_))
  }
  if (demean) {
    x <- x - mean(x, na.rm = TRUE)
  }
  mat <- embed(x^2, lags + 1)
  fit <- try(lm(mat[, 1] ~ mat[, -1]), silent = TRUE)
  if ("try-error" %in% class(fit)) {
    return(c(ARCH.LM = NA_real_))
  } else {
    arch.lm <- summary(fit)
    S <- arch.lm$r.squared #* NROW(mat)
    return(c(ARCH.LM = if(is.nan(S)) 1 else S))
  }
}
