#' Multiple seasonality features
#'
#' Computes various measures of multiple seasonality of a time series. The number
#' of seasonal periods, and the length of the seasonal periods are return. Also,
#' the strength of seasonality corresponding to each period is estimated using
#' a linear model involving Fourier terms corresponding to the seasonal periods.
#' This is analogous to the \code{stl_features} approach, but using a linear 
#' decomposition rather than an STL decomposition.
#' @param x a univariate time series of class \code{msts} or \code{ts}
#' @param transform A logical variable indicating if a Box-Cox transform should be applied
#' before the STL decomposition. If \code{lambda} is not NULL, then \code{transform} is set to TRUE.
#' @param lambda The value of the Box-Cox transformation parameter.
#' If \code{lambda=NULL} and \code{transform=TRUE}, then lambda is automatically selected
#' using \code{\link[forecast]{BoxCox.lambda}}.
#' @param ... Other arguments are passed to \code{\link[forecast]{BoxCox.lambda}}.
#' @return A vector of numeric values.
#' @export

multiseasonal <- function(x, transform=FALSE, lambda=NULL, ...)
{
  if("msts" %in% class(x))
  {
    msts <- attributes(x)$msts
    nperiods <- length(msts)
  }
  else if("ts" %in% class(x))
  {
    msts <- frequency(x)
    nperiods <- msts > 1
    season <- 0
  }
  else
  {
    msts <- 1
    nperiods <- 0L
    season <- 0
  }
  if(nperiods > 0.1)
  {
    # Now fit a linear model using Fourier terms to estimate seasonal patterns 
    if(transform | !is.null(lambda))
    {
      if(is.null(lambda))
        lambda <- forecast::BoxCox.lambda(x, ...)
      x <- forecast::BoxCox(x, lambda=lambda)
    }
    n <- length(x)
    tt <- seq_len(n)
    order <- pmin(trunc(msts/2), 50L)
    X <- fourier(x, K=order)
    fit <- lm(x ~ X + poly(tt,min(10L, round(n/20))))
    # Compute seasonal components
    seas <- list()
    up <- cumsum(c(1,order))
    for(i in seq_along(msts))
      seas[[i]] <- c(X[,up[i]:(up[i+1]-1)] %*% matrix(coef(fit)[up[i]+seq_len(order[i])]))
    # Compute strength of each seasonal component
    remainder <- residuals(fit)
    vare <- var(remainder, na.rm=TRUE)
    season <- numeric(length(msts))
    for(i in seq_along(msts))
      season[i] <- max(0, min(1, 1 - vare/var(remainder+seas[[i]], na.rm=TRUE)))
  }
  return(c(
    nperiods = nperiods,
    seasonal_period = msts,
    seasonal_strength = season))
}