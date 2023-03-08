
#' Time series features based on tiled windows
#'
#' Computes feature of a time series based on tiled (non-overlapping) windows.
#' Means or variances are produced for all tiled windows. Then stability is
#' the variance of the means, while lumpiness is the variance of the variances.
#' @param x a univariate time series
#' @param width size of sliding window
#' @return A numeric vector of length 2 containing a measure of lumpiness and
#' a measure of stability.
#' @author Earo Wang and Rob J Hyndman
#' @export

lumpiness <- function(x, width = ifelse(frequency(x) > 1,
                        frequency(x), 10
                      )) {
  x <- scalets(x)
  nr <- length(x)
  lo <- seq(1, nr, by = width)
  up <- seq(width, nr + width, by = width)
  nsegs <- nr / width
  varx <- map_dbl(seq_len(nsegs), function(idx)
    var(x[lo[idx]:up[idx]], na.rm = TRUE))
  if (length(x) < 2 * width) {
    lumpiness <- 0
  } else {
    lumpiness <- var(varx, na.rm = TRUE)
  }
  return(c(lumpiness = lumpiness))
}

#' @rdname lumpiness
#' @export

stability <- function(x, width = ifelse(frequency(x) > 1,
                        frequency(x), 10
                      )) {
  x <- scalets(x)
  nr <- length(x)
  lo <- seq(1, nr, by = width)
  up <- seq(width, nr + width, by = width)
  nsegs <- nr / width
  meanx <- map_dbl(seq_len(nsegs), function(idx)
    mean(x[lo[idx]:up[idx]], na.rm = TRUE))
  if (length(x) < 2 * width) {
    stability <- 0
  } else {
    stability <- var(meanx, na.rm = TRUE)
  }
  return(c(stability = stability))
}

#' Time series features based on sliding windows
#'
#' Computes feature of a time series based on sliding (overlapping) windows.
#' \code{max_level_shift} finds the largest mean shift between two consecutive windows.
#' \code{max_var_shift} finds the largest var shift between two consecutive windows.
#' \code{max_kl_shift} finds the largest shift in Kulback-Leibler divergence between
#' two consecutive windows.
#'
#' Computes the largest level shift and largest variance shift in sliding mean calculations
#' @param x a univariate time series
#' @param width size of sliding window
#' @return A vector of 2 values: the size of the shift, and the time index of the shift.
#' @author Earo Wang and Rob J Hyndman
#' @export

max_level_shift <- function(x, width = ifelse(frequency(x) > 1,
                              frequency(x), 10
                            )) {
  suppressWarnings(rollmean <- try(RcppRoll::roll_mean(x, width, na.rm = TRUE), silent = TRUE))
  if ("try-error" %in% class(rollmean)) {
    maxmeans <- NA_real_
    maxidx <- NA_real_
  } else {
    means <- abs(diff(rollmean, width))
    if (length(means) == 0L) {
      maxmeans <- 0
      maxidx <- NA_real_
    }
    else if (all(is.na(means))) {
      maxmeans <- NA_real_
      maxidx <- NA_real_
    }
    else {
      maxmeans <- max(means, na.rm = TRUE)
      maxidx <- which.max(means) + width - 1L
    }
  }
  return(c(max_level_shift = maxmeans, time_level_shift = maxidx))
}

#' @rdname max_level_shift
#' @export

max_var_shift <- function(x, width = ifelse(frequency(x) > 1,
                            frequency(x), 10
                          )) {
  suppressWarnings(rollvar <- try(RcppRoll::roll_var(x, width, na.rm = TRUE), silent = TRUE))
  if ("try-error" %in% class(rollvar)) {
    maxvar <- NA_real_
    maxidx <- NA_real_
  } else {
    vars <- abs(diff(rollvar, width))

    if (length(vars) == 0L) {
      maxvar <- 0
      maxidx <- NA_real_
    }
    else if (all(is.na(vars))) {
      maxvar <- NA_real_
      maxidx <- NA_real_
    }
    else {
      maxvar <- max(vars, na.rm = TRUE)
      maxidx <- which.max(vars) + width - 1L
    }
  }
  return(c(max_var_shift = maxvar, time_var_shift = maxidx))
}

#' @rdname max_level_shift
#' @export

max_kl_shift <- function(x, width = ifelse(frequency(x) > 1,
                           frequency(x), 10
                         )) {
  gw <- 100 # grid width
  xgrid <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length = gw)
  grid <- xgrid[2L] - xgrid[1L]
  tmpx <- x[!is.na(x)] # Remove NA to calculate bw
  bw <- bw.nrd0(tmpx)
  lenx <- length(x)
  if (lenx <= (2 * width)) {
    return(c(max_kl_shift = NA_real_, time_kl_shift = NA_real_))
  }
  # Using binning algorithm to achieve efficiency but obsecure exact positions.
  # lastrep <- ceiling(lenx/5)
  # group <- rep(1:lastrep, each = 5)[1:lenx]
  # midpoints <- aggregate(x, by = list(group), function(y) y[3L])[, 2]
  # dens.mat <- matrix(, nrow = lastrep, ncol = gw)
  # for (i in 1L:lastrep) {
  #   dens.mat[i, ] <- dnorm(xgrid, mean = midpoints[i], sd = bw)
  # }
  dens.mat <- matrix(, nrow = lenx, ncol = gw)
  for (i in 1L:lenx) {
    dens.mat[i, ] <- dnorm(xgrid, mean = x[i], sd = bw)
  }
  dens.mat <- pmax(dens.mat, dnorm(38))
  rmean <- RcppRoll::roll_mean(dens.mat,
    n = width, na.rm = TRUE, fill = NA,
    align = "right"
  ) # by column
  # lo <- seq(1, lastrep - width + 1)
  # hi <- seq(width + 1, lastrep)
  lo <- seq(1, lenx - width + 1)
  hi <- seq(width + 1, lenx)
  seqidx <- min(length(lo), length(hi))
  kl <- sapply(1:seqidx, function(i) sum(rmean[lo[i], ] *
      (log(rmean[lo[i], ]) - log(rmean[hi[i], ])) *
      grid, na.rm = TRUE))
  diffkl <- diff(kl, na.rm = TRUE)
  if (length(diffkl) == 0L) {
    diffkl <- 0
    maxidx <- NA_real_
  }
  else {
    maxidx <- which.max(diffkl) + width - 1L
  }
  return(c(max_kl_shift = max(diffkl, na.rm = TRUE), time_kl_shift = maxidx))
}

#' Number of crossing points
#'
#' Computes the number of times a time series crosses the median.
#' @param x a univariate time series
#' @return A numeric value.
#' @author Earo Wang and Rob J Hyndman
#' @export
crossing_points <- function(x) {
  midline <- median(x, na.rm = TRUE)
  ab <- x <= midline
  lenx <- length(x)
  p1 <- ab[1:(lenx - 1)]
  p2 <- ab[2:lenx]
  cross <- (p1 & !p2) | (p2 & !p1)
  return(c(crossing_points = sum(cross, na.rm = TRUE)))
}

#' Longest flat spot
#'
#' "Flat spots” are computed by dividing the sample space of a time series into ten equal-sized intervals, and computing the maximum run length within any single interval.
#' @param x a univariate time series
#' @return A numeric value.
#' @author Earo Wang and Rob J Hyndman
#' @export

flat_spots <- function(x) {
  cutx <- try(cut(x, breaks = 10, include.lowest = TRUE, labels = FALSE),
    silent = TRUE
  )
  if ("try-error" %in% class(cutx)) {
    fspots <- NA
  } else {
    rlex <- rle(cutx)
    # Any flat spot
    return(c(flat_spots = max(rlex$lengths)))
    # Low flat spots
    # ones <- (rlex$values == 1)
    # return(max(rlex$lengths[ones]))
  }
}

# shapes <- function(x, width, scale = TRUE, FUN = mean, ...){
#   nr <- length(x)
#   if (nr %% width != 0) {
#       stop("width must be a divisor of the length of the series.")
#     }
#     shapes <- matrix(x, ncol = width, byrow= TRUE)
#     if(scale){
#       dtotal <- apply(shapes, 1, sum)
#       idremove <- which(dtotal == 0)
#       if(length(idremove) > 0){
#         shapes <- shapes[-idremove, ]
#         dtotal <- dtotal[-idremove]
#       }
#       shapes <- t(t(shapes) / dtotal)
#     }
#     xprofile <- apply(shapes, 2, FUN, ...)
#     return(c(shapes=xprofile))
# }

#' Hurst coefficient
#'
#' Computes the Hurst coefficient indicating the level of fractional differencing
#' of a time series.
#' @param x a univariate time series. If missing values are present, the largest
#' contiguous portion of the time series is used.
#' @return A numeric value.
#' @author Rob J Hyndman
#' @export

hurst <- function(x) {
  # Hurst=d+0.5 where d is fractional difference.
  return(c(hurst = suppressWarnings(fracdiff::fracdiff(na.contiguous(x), 0, 0)[["d"]] + 0.5)))
}

#' Unit Root Test Statistics
#'
#' \code{unitroot_kpss} computes the statistic for the Kwiatkowski et al. unit root test
#' using the default settings for the \code{\link[urca]{ur.kpss}} function.
#' \code{unitroot_pp} computes the statistic for the Phillips-Perron unit root test
#' using the default settings for the \code{\link[urca]{ur.pp}} function.
#' @param x a univariate time series.
#' @param ... Other arguments are passed to the \code{\link[urca]{ur.kpss}} or
#' \code{\link[urca]{ur.kpss}} functions.
#' @return A numeric value
#' @author Pablo Montero-Manso
#' @export
unitroot_kpss <- function(x, ...) {
  kpss <- try(urca::ur.kpss(x, ...)@teststat, silent=TRUE)
  if("try-error" %in% class(kpss)) {
    warning("Error in unitroot_kpss")
    kpss <- NA
  }
  return(kpss)
}

#' @rdname unitroot_kpss
#' @export
unitroot_pp <- function(x, ...) {
  pp <- try(urca::ur.pp(x, ...)@teststat, silent = TRUE)
  if("try-error" %in% class(pp)) {
    warning("Error in unitroot_pp")
    pp <- NA
  }
  return(pp)
}
