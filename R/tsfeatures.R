#' Time Series Feature Extraction
#'
#' The tsfeature package provides methods to extract various features from time series data
#'
#' @docType package
#' @name tsfeatures
#' @importFrom stats as.ts bw.nrd0 coef dnorm embed fitted frequency lm spec.ar
#' @importFrom stats median na.contiguous na.pass residuals cor sd tsp "tsp<-" var
#' @importFrom stats quantile acf pacf stl pchisq ar Box.test poly start cmdscale
#' @importFrom purrr map map_dbl
#' @importFrom forecast mstl
NULL
# > NULL

#' Convert mts object to list of time series
#' @method as.list mts
#' @param x multivariate time series of class mts.
#' @param ... other arguments are ignored.
#' @author Rob J Hyndman
#' @export
as.list.mts <- function(x, ...) {
  tspx <- tsp(x)
  listx <- as.list(as.data.frame(x))
  listx <- purrr::map(
    listx,
    function(u) {
      u <- as.ts(u)
      tsp(u) <- tspx
      return(u)
    }
  )
  return(listx)
}
