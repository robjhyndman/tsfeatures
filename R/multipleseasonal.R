
#' Strength of trend and seasonality of a time series
#'
#' Computes various measures of trend and seasonality of a time series based on
#' an STL decomposition. The number of seasonal periods, and the length of the
#' seasonal periods are returned. Also, the strength of seasonality corresponding
#' to each period is estimated. The \code{\link[forecast]{mstl}} function is used
#' to do the decomposition.
#' @param x a univariate time series.
#' @param ... Other arguments are passed to \code{\link[forecast]{mstl}}.
#' @return A vector of numeric values.
#' @author Rob J Hyndman
#' @export

stl_features <- function(x, ...) {
  if ("msts" %in% class(x)) {
    msts <- attributes(x)$msts
    nperiods <- length(msts)
  }
  else if ("ts" %in% class(x)) {
    msts <- frequency(x)
    nperiods <- msts > 1
    season <- 0
  }
  else {
    msts <- 1
    nperiods <- 0L
    season <- 0
  }
  trend <- linearity <- curvature <- season <- spike <- peak <- trough <- acfremainder <- NA

  # STL fits
  stlfit <- forecast::mstl(x, ...)
  trend0 <- stlfit[, "Trend"]
  remainder <- stlfit[, "Remainder"]
  seasonal <- stlfit[, grep("Season", colnames(stlfit)), drop = FALSE]

  # When the maximum frequency is dropped
  tsp(x) <- tsp(trend0)

  # De-trended and de-seasonalized data
  detrend <- x - trend0
  deseason <- forecast::seasadj(stlfit)
  fits <- x - remainder

  # Summary stats
  n <- length(x)
  varx <- var(x, na.rm = TRUE)
  vare <- var(remainder, na.rm = TRUE)
  vardetrend <- var(detrend, na.rm = TRUE)
  vardeseason <- var(deseason, na.rm = TRUE)
  nseas <- NCOL(seasonal)

  # Measure of trend strength
  if (vardeseason / varx < 1e-10) {
    trend <- 0
  } else {
    trend <- max(0, min(1, 1 - vare / vardeseason))
  }

  if (nseas > 0) {
    # Measure of seasonal strength
    season <- numeric(nseas)
    for (i in seq(nseas))
      season[i] <- max(0, min(1, 1 - vare / var(remainder + seasonal[, i], na.rm = TRUE)))

    # Find time of peak and trough for each component
    peak <- trough <- numeric(nseas)
    for (i in seq(nseas))
    {
      startx <- start(x)[2L] - 1L
      pk <- (startx + which.max(seasonal[, i])) %% msts[i]
      th <- (startx + which.min(seasonal[, i])) %% msts[i]
      peak[i] <- ifelse(pk == 0, msts[i], pk)
      trough[i] <- ifelse(th == 0, msts[i], th)
    }
  }

  # Compute measure of spikiness
  d <- (remainder - mean(remainder, na.rm = TRUE))^2
  varloo <- (vare * (n - 1) - d) / (n - 2)
  spike <- var(varloo, na.rm = TRUE)

  # Compute measures of linearity and curvature
  tren.coef <- coef(lm(trend0 ~ poly(seq(n), degree = 2L)))[2L:3L]
  linearity <- tren.coef[1L]
  curvature <- tren.coef[2L]

  # ACF of remainder
  acfremainder <- unname(acf_features(remainder))

  # Assemble results
  output <- c(
    nperiods = nperiods, seasonal_period = msts, trend = trend,
    spike = spike, linearity = unname(linearity), curvature = unname(curvature),
    e_acf1 = acfremainder[1L], e_acf10 = acfremainder[2L]
  )
  if (nseas > 0) {
    output <- c(output, seasonal_strength = season, peak = peak, trough = trough)
  }

  return(output)
}
