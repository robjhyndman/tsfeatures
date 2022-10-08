#' Spectral entropy of a time series
#'
#' @description
#' Computes spectral entropy from a univariate normalized
#' spectral density, estimated using an AR model.
#'
#' @details
#' The \emph{spectral entropy} equals the Shannon entropy of the spectral density
#' \eqn{f_x(\lambda)} of a stationary process \eqn{x_t}:
#' \deqn{
#' H_s(x_t) = - \int_{-\pi}^{\pi} f_x(\lambda) \log f_x(\lambda) d \lambda,
#' }
#' where the density is normalized such that
#' \eqn{\int_{-\pi}^{\pi} f_x(\lambda) d \lambda = 1}.
#' An estimate of \eqn{f(\lambda)} can be obtained using \code{\link[stats]{spec.ar}} with
#' the `burg` method.
#'
#' @param x a univariate time series
#' @author Rob J Hyndman
#' @return
#' A non-negative real value for the spectral entropy \eqn{H_s(x_t)}.
#' @seealso \code{\link[stats]{spec.ar}}
#' @references
#' Jerry D. Gibson and Jaewoo Jung (2006). \dQuote{The
#' Interpretation of Spectral Entropy Based Upon Rate Distortion Functions}.
#' IEEE International Symposium on Information Theory, pp. 277-281.
#'
#' Goerg, G. M. (2013). \dQuote{Forecastable Component Analysis}.
#' Proceedings of the 30th International Conference on Machine Learning (PMLR) 28 (2): 64-72, 2013.
#' Available at \url{https://proceedings.mlr.press/v28/goerg13.html}.
#'
#' @examples
#' entropy(rnorm(1000))
#' entropy(lynx)
#' entropy(sin(1:20))
#' @export

entropy <- function(x) {
  #spec <- spectrum(x, plot = FALSE, n.freq = ceiling(length(x)/2 + 1), ...)
  spec <- try(stats::spec.ar(na.contiguous(x), plot=FALSE, method='burg',
                      n.freq = ceiling(length(x)/2 + 1)))
  if ("try-error" %in% class(spec)) {
    entropy <- NA
  } else {
    fx <- c(rev(spec$spec[-1]),spec$spec)/ length(x)
    fx <- fx/sum(fx)
    prior.fx = rep(1 / length(fx), length = length(fx))
    prior.weight = 0.001
    fx <- (1 - prior.weight) * fx + prior.weight * prior.fx
    entropy <- pmin(1, -sum(fx * log(fx, base = length(x))))
  }
  return(c(entropy = entropy))
}
