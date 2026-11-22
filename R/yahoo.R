#' Yahoo server metrics
#'
#' Downloads and returns aggregated and anonymized datasets from Yahoo representing server metrics of Yahoo services.
#' If the download fails (e.g., no internet connection or the server is unavailable),
#' a message is shown and `NULL` is returned invisibly.
#'
#' @param ... Additional arguments passed to `download.file`
#'
#' @return A matrix of time series with 1437 rows of hourly data, and 1748 columns representing different servers,
#' or `NULL` if the data could not be downloaded.
#' @author Rob Hyndman, Earo Wang, Nikolay Laptev, Mitchell O'Hara-Wild
#' @references
#' Hyndman, R.J., Wang, E., Laptev, N. (2015) Large-scale unusual time series detection.
#' In: \emph{Proceedings of the IEEE International Conference on Data Mining}. Atlantic City, NJ, USA. 14–17 November 2015.
#' \url{https://robjhyndman.com/publications/icdm2015/}
#' @examples
#' yahoo <- yahoo_data()
#' if (!is.null(yahoo)) {
#'   plot(yahoo[, 1:10])
#'   plot(yahoo[, 1:44], plot.type = "single", col = 1:44)
#' }
#' @export
yahoo_data <- function(...) {
  tmp <- tempfile(fileext = ".rda")
  on.exit(unlink(tmp), add = TRUE)
  ok <- tryCatch(
    {
      utils::download.file(
        "https://github.com/robjhyndman/tsfeatures/raw/master/extra-data/yahoo.rda",
        tmp,
        mode = "wb",
        ...
      )
      TRUE
    },
    error = function(e) {
      message("Unable to download the Yahoo data: ", conditionMessage(e))
      FALSE
    },
    warning = function(w) {
      message("Unable to download the Yahoo data: ", conditionMessage(w))
      FALSE
    }
  )
  if (!ok) {
    return(invisible(NULL))
  }
  load(tmp)
  yahoo
}

# `yahoo` is created by load() inside yahoo_data()
utils::globalVariables("yahoo")
