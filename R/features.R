#' Spectral entropy of a time series
#'
#' Computes the spectral entropy of a time series
#' @param x a univariate time series
#' @return A numeric value.
#' @export

entropy <- function(x) {
  entropy <- try(ForeCA::spectral_entropy(na.contiguous(x))[1L], silent = TRUE)
  if (class(entropy) == "try-error") {
    entropy <- NA
  }
  return(c(entropy=entropy))
}

#' Autocorrelation function at lag 1
#'
#' Computes the first order autocorrelation a time series.
#' @param x a univariate time series
#' @return A numeric value.
#' @export

acf1 <- function(x) {
  return(c(acf1=acf(x, lag.max=1L, plot = FALSE, na.action = na.pass)$acf[2L]))
}

#' Strength of trend and seasonality of a time series
#'
#' Computes the strength of trend and seasonality of a time series using an STL
#' decomposition (or for non-seasonal time series, a penalized regression spline smoother).
#' @param x a univariate time series
#' @param robust A logical variable indicating if robust STL should be applied.
#' @param transform A logical variable indicating if a Box-Cox transform should be applied
#' before the STL decomposition.
#' @param lambda The value of the Box-Cox transformation parameter if \code{transform=TRUE}.
#' If \code{lambda=NULL}, it is automatically selected using \code{\link[forecast]{BoxCox.lambda}}.
#' @param ... Other arguments are passed to \code{\link[forecast]{BoxCox.lambda}}.
#' @return A numeric value.
#' @export

stl_features <- function(x, robust=FALSE, transform=FALSE, lambda=NULL,...) 
{
  x <- as.ts(x)
  tspx <- tsp(x)
  freq <- tspx[3]
  contx <- try(na.contiguous(x), silent = TRUE)
  len.contx <- length(contx)
  varx <- var(contx, na.rm=TRUE)
  tt <- seq_len(len.contx)
  if (length(contx) <= 2 * freq || class(contx) == "try-error") 
  {
    trend <- linearity <- curvature <- season <- spike <- peak <- trough <- acfremainder <- NA
  } 
  else 
  {
    if(transform)
    {
      if(is.null(lambda))
        lambda <- forecast::BoxCox.lambda(contx, ...)
      contx <- forecast::BoxCox(contx, lambda=lambda)
    }
    if (freq > 1L) 
    {
      stlfit <- stl(contx, s.window = "periodic", robust = robust)
      trend0 <- stlfit$time.series[, "trend"]
      seasonal <- stlfit$time.series[, "seasonal"]
      remainder <- stlfit$time.series[, "remainder"]
      # Find time of peak and trough
      starty <- start(contx)[2L]
      pk <- (starty + which.max(seasonal) - 1L) %% freq
      th <- (starty + which.min(seasonal) - 1L) %% freq
      peak <- ifelse(pk == 0, freq, pk) * max(seasonal)
      trough <- ifelse(th == 0, freq, th) * min(seasonal)
      # Compute detrended and deseasonalised data and fitted values
      detrend <- contx - trend0
      deseason <- contx - seasonal
      fits <- trend0 + seasonal
      # Compute measure of seasonal strength
      v.adj <- var(remainder, na.rm = TRUE)
      vardetrend <- var(detrend, na.rm=TRUE)
      if(vardetrend / varx < 1e-10)
        season <- 0
      else
        season <- max(0, min(1, 1 - v.adj/vardetrend))
    }
    else 
    { # No seasonal component, so use spline for trend
      trend0 <- fitted(mgcv::gam(contx ~ s(tt)))
      remainder <- contx - trend0
      deseason <- contx
      v.adj <- var(remainder, na.rm = TRUE)
    }
    # Compute measure of trend strength
    vardeseason <- var(deseason, na.rm=TRUE)
    if(vardeseason / varx < 1e-10)
      trend <- 0
    else
      trend <- max(0, min(1, 1 - v.adj/vardeseason))
    # Compute measure of spikiness
    n <- length(remainder)
    d <- (remainder - mean(remainder, na.rm = TRUE))^2
    varloo <- (v.adj*(n-1) - d) / (n-2)
    spike <- var(varloo, na.rm = TRUE)
    # Compute measures of linearity and curvature
    pl <- poly(tt, degree = 2L)
    tren.coef <- coef(lm(trend0 ~ pl))[2L:3L]
    linearity <- tren.coef[1L]
    curvature <- tren.coef[2L]
    # ACF lag 1 of remainder
    acfremainder <- acf1(remainder)
  }
  output <- c(trend = trend, spike = spike, linearity = unname(linearity),
                curvature = unname(curvature), acfremainder=unname(acfremainder)) 
  if (freq > 1) 
    output <- c(output, season = season, peak = peak, trough = trough)

  return(output)
}


#' Time series features based on tiled windows
#'
#' Computes feature of a time series based on tiled (non-overlapping) windows.
#' Means or variances are produced for all tiled windows. Then stability is
#' the variance of the means, while lumpiness is the variance of the variances.
#' @param x a univariate time series
#' @param width size of sliding window
#' @return A numeric vector of length 2 containing a measure of lumpiness and
#' a measure of stability.
#' @export

lumpiness <- function(x, width=ifelse(frequency(x) > 1,
                       frequency(x), 10))
{
  x <- scalets(x)
  nr <- length(x)
  lo <- seq(1, nr, by = width)
  up <- seq(width, nr + width, by = width)
  nsegs <- nr/width
  varx <- map_dbl(seq_len(nsegs), function(idx)
                 var(x[lo[idx]:up[idx]], na.rm = TRUE))
  return(c(lumpiness = var(varx, na.rm = TRUE)))
}

#' @rdname lumpiness
#' @export

stability <- function(x, width=ifelse(frequency(x) > 1,
                       frequency(x), 10))
{
  x <- scalets(x)
  nr <- length(x)
  lo <- seq(1, nr, by = width)
  up <- seq(width, nr + width, by = width)
  nsegs <- nr/width
  meanx <- map_dbl(seq_len(nsegs), function(idx)
                 mean(x[lo[idx]:up[idx]], na.rm = TRUE))
  return(c(stability = var(meanx, na.rm = TRUE)))
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
#' @export

max_level_shift <- function(x, width=ifelse(frequency(x) > 1,
                       frequency(x), 10))
{
  rollmean <- RcppRoll::roll_mean(x, width, na.rm = TRUE)
  means <- abs(diff(rollmean, width))
  if(all(is.na(means)))
  {
    maxmeans <- NA_real_
    maxidx <- NA_real_
  }
  else
  {
    maxmeans <- max(means, na.rm=TRUE)
    maxidx <- which.max(means)+1L
  }
  return(c(max_level_shift=maxmeans, time_level_shift=maxidx))
}

#' @rdname max_level_shift
#' @export

max_var_shift <- function(x, width=ifelse(frequency(x) > 1,
                       frequency(x), 10))
{
  rollvar <- RcppRoll::roll_var(x, width, na.rm = TRUE)
  vars <- abs(diff(rollvar, width))
  if(all(is.na(vars)))
  {
    maxvar <- NA_real_
    maxidx <- NA_real_
  }
  else
  {
    maxvar <- max(vars, na.rm=TRUE)
    maxidx <- which.max(vars)+1L
  }
  return(c(max_var_shift=maxvar, time_var_shift=maxidx))
}

#' @rdname max_level_shift
#' @export

max_kl_shift <- function(x, width=ifelse(frequency(x) > 1,
                       frequency(x), 10))
{
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
  rmean <- RcppRoll::roll_mean(dens.mat, n = width, na.rm = TRUE, fill = NA,
                               align = "right") # by column
  # lo <- seq(1, lastrep - width + 1)
  # hi <- seq(width + 1, lastrep)
  lo <- seq(1, lenx - width + 1)
  hi <- seq(width + 1, lenx)
  seqidx <- min(length(lo), length(hi))
  kl <- sapply(1:seqidx, function(i) sum(rmean[lo[i], ] *
               (log(rmean[lo[i], ]) - log(rmean[hi[i], ])) *
               grid, na.rm = TRUE))
  diffkl <- diff(kl, na.rm = TRUE)
  maxidx <- which.max(diffkl) + 1L
  return(c(max_kl_shift = max(diffkl), time_kl_shift = maxidx))
}


#' Number of crossing points
#'
#' Computes the number of times a time series crosses the median.
#' @param x a univariate time series
#' @return A numeric value.
#' @export
crossing_points <- function(x)
{
  midline <- median(x, na.rm = TRUE)
  ab <- x <= midline
  lenx <- length(x)
  p1 <- ab[1:(lenx - 1)]
  p2 <- ab[2:lenx]
  cross <- (p1 & !p2) | (p2 & !p1)
  return(c(crossing_points=sum(cross, na.rm=TRUE)))
}

#' Number of flat spots
#'
#' Number of flat spots in a time series
#' @param x a univariate time series
#' @return A numeric value.
#' @export

flat_spots <- function(x) {
  cutx <- try(cut(x, breaks = 10, include.lowest = TRUE, labels = FALSE),
              silent = TRUE)
  if (class(cutx) == "try-error") {
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


#' Hurst coefficent
#'
#' Computes the Hurst coefficient indicating the level of fractional differencing
#' of a time series.
#' @param x a univariate time series. If missing values are present, the largest 
#' contiguous portion of the time series is used.
#' @return A numeric value.
#' @export

hurst <- function(x)
{
  # Hurst=d+0.5 where d is fractional difference.
  return(c(hurst = fracdiff::fracdiff(na.contiguous(x),0,0)[["d"]] + 0.5))
}

#' Lyapunov coefficent
#'
#' Computes the Lyapunov coefficient indicating the level of nonlinearity of
#' of a time series.
#' @param x a univariate time series
#' @return A numeric value.
#' @export

lyapunov <- function(x)
{
  freq <- frequency(x)
  n <- length(x) 
  if(freq >= n)
    stop("Insufficient data")
  Ly <- rep(NA_real_, n-freq)
  for(i in seq_len(n-freq))
  {
    idx <- order(abs(x[i] - x))
    idx <- idx[idx < (n-freq)]
    j <- idx[2]
    Ly[i] <- abs((x[i+freq] - x[j+freq])/(x[i]-x[j]))
    if(Ly[i] > 0)
      Ly[i] <- log(Ly[i])/freq
  }
  return(lyapunov = mean(Ly,na.rm=TRUE))
}

