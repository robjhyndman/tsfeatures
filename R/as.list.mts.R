#' Convert mts object to list of time series
#' 
#' An mts object contains a multivariate time series in a matrix, with time on rows.
#' This is converted into a list of univariate time series.
#' 
#' @method as.list mts
#' @param x multivariate time series of class mts.
#' @param ... other arguments are ignored.
#' @author Rob J Hyndman
#' @return A list of ts objects.
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
