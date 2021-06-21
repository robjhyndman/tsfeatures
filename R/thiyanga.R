#' Autocorrelation-based features
#'
#' Computes various measures based on autocorrelation coefficients of the
#' original series, first-differenced series and second-differenced series
#' @param x a univariate time series
#' @return A vector of 6 values: first autocorrelation coefficient and sum of squared of
#' first ten autocorrelation coefficients of original series, first-differenced series,
#' and twice-differenced series.
#' For seasonal data, the autocorrelation coefficient at the first seasonal lag is
#' also returned.
#' @author Thiyanga Talagala
#' @export
acf_features <- function(x) {
  m <- frequency(x)
  if(length(x) > 1) {
    acfx <- acf(x, lag.max = max(10L, m), plot = FALSE, na.action=na.pass)$acf[-1L]
  } else {
    acfx <- NA
  }
  if(length(x) > 10) {
    acfdiff1x <- acf(diff(x, differences = 1), lag.max = 10L, plot = FALSE, na.action = na.pass)$acf[-1L]
  } else {
    acfdiff1x <- NA
  }
  if(length(x) > 11) {
    acfdiff2x <- acf(diff(x, differences = 2), lag.max = 10L, plot = FALSE, na.action = na.pass)$acf[-1L]
  } else {
    acfdiff2x <- NA
  }

  # first autocorrelation coefficient
  acf_1 <- acfx[1L]

  # sum of squares of first 10 autocorrelation coefficients
  sum_of_sq_acf10 <- sum((acfx[seq(10)])^2)

  # first autocorrelation coefficient of differenced series
  diff1_acf1 <- acfdiff1x[1L]

  # Sum of squared of first 10 autocorrelation coefficients of differenced series
  diff1_acf10 <- sum((acfdiff1x[seq(10)])^2)

  # first autocorrelation coefficient of twice-differenced series
  diff2_acf1 <- acfdiff2x[1L]

  # Sum of squared of first 10 autocorrelation coefficients of twice-differenced series
  diff2_acf10 <- sum((acfdiff2x[seq(10)])^2)

  output <- c(
    x_acf1 = unname(acf_1),
    x_acf10 = unname(sum_of_sq_acf10),
    diff1_acf1 = unname(diff1_acf1),
    diff1_acf10 = unname(diff1_acf10),
    diff2_acf1 = unname(diff2_acf1),
    diff2_acf10 = unname(diff2_acf10)
  )

  if (m > 1) {
    output <- c(output, seas_acf1 = unname(acfx[m]))
  }

  return(output)
}

#' Partial autocorrelation-based features
#'
#' Computes various measures based on partial autocorrelation coefficients of the
#' original series, first-differenced series and second-differenced series
#' @param x a univariate time series
#' @return A vector of 3 values: Sum of squared of first 5
#' partial autocorrelation coefficients of the original series, first differenced
#' series and twice-differenced series.
#' For seasonal data, the partial autocorrelation coefficient at the first seasonal
#' lag is also returned.
#' @author Thiyanga Talagala
#' @export
pacf_features <- function(x) {
  m <- frequency(x)
  if(length(x) > 1){
    pacfx <- pacf(x, lag.max = max(5L, m), plot = FALSE)$acf
  } else {
    pacfx <- NA
  }

  # Sum of first 5 PACs squared
  if(length(x) > 5) {
    pacf_5 <- sum((pacfx[seq(5L)])^2)
  } else {
    pacf_5 <- NA
  }

  # Sum of first 5 PACs of difference series squared
  if(length(x) > 6) {
    diff1_pacf_5 <- sum(pacf(diff(x, differences = 1L), lag.max = 5L, plot = FALSE)$acf^2)
  } else {
    diff1_pacf_5 <- NA
  }

  # Sum of first 5 PACs of twice differenced series squared
  if(length(x) > 7) {
    diff2_pacf_5 <- sum(pacf(diff(x, differences = 2L), lag.max = 5L, plot = FALSE)$acf^2)
  } else {
    diff2_pacf_5 <- NA
  }

  output <- c(
    x_pacf5 = unname(pacf_5),
    diff1x_pacf5 = unname(diff1_pacf_5),
    diff2x_pacf5 = unname(diff2_pacf_5)
  )
  if (m > 1) {
    output <- c(output, seas_pacf = pacfx[m])
  }
  return(output)
}

#' Parameter estimates of Holt's linear trend method
#'
#' Estimate the smoothing parameter for the level-alpha and
#' the smoothing parameter for the trend-beta.
#' \code{hw_parameters} considers additive seasonal trend: ets(A,A,A) model.
#' @param x a univariate time series
#' @return \code{holt_parameters} produces a vector of 2 values: alpha, beta.
#'
#' \code{hw_parameters} produces a vector of 3 values: alpha, beta and gamma.
#' @author Thiyanga Talagala, Pablo Montero-Manso
#' @export

holt_parameters <- function(x) {
  # parameter estimates of holt linear trend model
  fit <- forecast::ets(x, model = c("AAN"))
  params <- c(fit$par["alpha"], fit$par["beta"])
  names(params) <- c("alpha", "beta")
  return(params)
}

#' @rdname holt_parameters
#' @export
hw_parameters <- function(x) {
  # parameter estimates of holt winters additive trend seasonal model
  hw_fit <- purrr::possibly(forecast::ets,
    list(par = c(alpha = NA, beta = NA, gamma = NA)))(x, model = c("AAA"))
  return(hw_fit$par[c("alpha", "beta", "gamma")])
}
# #' Autocorrelation coefficient at lag 1 of the residual
# #'
# #' Computes the first order autocorrelation of the residual series of the deterministic trend model
# #' @param x a univariate time series
# #' @return A numeric value.
# #' @author Thiyanga Talagala
# #' @export
# acfresid <- function(x){
#   time <- 1:length(x)
#   linear_mod <- lm(x~time)
#   Res<-resid(linear_mod)
#   return(stats::acf(Res,lag.max=1L,plot=FALSE)$acf[-1])
# }

#' Proportion of zeros
#'
#' Computes proportion of zeros in a time series
#' @param x a univariate time series
#' @param tol tolerance level. Absolute values below this are considered zeros.
#' @return A numeric value.
#' @author Thiyanga Talagala
#' @export
zero_proportion <- function(x, tol = 1e-8) {
  mean(abs(x) < tol, na.rm=TRUE)
}
