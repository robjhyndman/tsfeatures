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
  acfx <- acf(x, lag.max = max(m, 10L), plot = FALSE, na.action = na.pass)
  acfdiff1x <- acf(diff(x, differences = 1), lag.max = 10L, plot = FALSE, na.action = na.pass)
  acfdiff2x <- acf(diff(x, differences = 2), lag.max = 10L, plot = FALSE, na.action = na.pass)

  # first autocorrelation coefficient
  acf_1 <- acfx$acf[2L]

  # sum of squares of first 10 autocorrelation coefficients
  sum_of_sq_acf10 <- sum((acfx$acf[2L:11L])^2)

  # first autocorrelation coefficient of differenced series
  diff1_acf1 <- acfdiff1x$acf[2L]

  # Sum of squared of first 10 autocorrelation coefficients of differenced series
  diff1_acf10 <- sum((acfdiff1x$acf[-1L])^2)

  # first autocorrelation coefficient of twice-differenced series
  diff2_acf1 <- acfdiff2x$acf[2L]

  # Sum of squared of first 10 autocorrelation coefficients of twice-differenced series
  diff2_acf10 <- sum((acfdiff2x$acf[-1L])^2)

  output <- c(
    x_acf1 = unname(acf_1),
    x_acf10 = unname(sum_of_sq_acf10),
    diff1_acf1 = unname(diff1_acf1),
    diff1_acf10 = unname(diff1_acf10),
    diff2_acf1 = unname(diff2_acf1),
    diff2_acf10 = unname(diff2_acf10)
  )

  if (m > 1) {
    output <- c(output, seas_acf1 = unname(acfx$acf[m + 1L]))
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

  pacfx <- pacf(x, lag.max = max(5L, m), plot = FALSE)$acf
  # Sum of squared of first 5 partial autocorrelation coefficients
  pacf_5 <- sum((pacfx[seq(5L)])^2)

  # Sum of squared of first 5 partial autocorrelation coefficients of difference series
  diff1_pacf_5 <- sum((pacf(diff(x, differences = 1), lag.max = 5L, plot = FALSE)$acf)^2)

  # Sum of squared of first 5 partial autocorrelation coefficients of twice differenced series
  diff2_pacf_5 <- sum((pacf(diff(x, differences = 2), lag.max = 5L, plot = FALSE)$acf)^2)

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
  fit <- forecast::holt(x)
  return(c(fit$model$par["alpha"], fit$model$par["beta"]))
}

#' @rdname holt_parameters
#' @export
hw_parameters <- function(x) {
  hw_fit <- NULL
  hw_fit$par <- c(NA, NA, NA)
  try(hw_fit <- forecast::ets(x, model = c("AAA")), silent = TRUE)
  hw_fit$par[1:3]
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
