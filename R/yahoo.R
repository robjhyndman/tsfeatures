#' Yahoo server metrics
#'
#' @param ... Additional arguments passed to `download.file`
#'
#' Downloads and returns aggregated and anonymized datasets from Yahoo representing server metrics of Yahoo services.
#'
#' @return A matrix of time series with 1437 rows of hourly data, and 1748 columns representing different servers.
#' @author Rob Hyndman, Earo Wang, Nikolay Laptev, Mitchell O'Hara-Wild
#' @references
#' Hyndman, R.J., Wang, E., Laptev, N. (2015) Large-scale unusual time series detection.
#' In: \emph{Proceedings of the IEEE International Conference on Data Mining}. Atlantic City, NJ, USA. 14–17 November 2015.
#' \url{https://robjhyndman.com/publications/icdm2015/}
#' @examples
#' yahoo <- yahoo_data()
#' plot(yahoo[,1:10])
#' plot(yahoo[,1:44], plot.type='single', col=1:44)
#' 
#' @export
yahoo_data <- function(...){
  yahoo <- 1 # Just to avoid a note about undefined global variables.
  tmp <- tempfile()
  utils::download.file("https://github.com/robjhyndman/tsfeatures/raw/master/extra-data/yahoo.rda", tmp, ...)
  load(tmp)
  yahoo
}