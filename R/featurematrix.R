#' Time series feature matrix
#'
#' \code{tsfeatures} computes a matrix of time series features from a list of time series
#' @param tslist a list of univariate time series, each of class \code{ts} or a numeric vector.
#' Alternatively, an object of class \code{mts} may be used.
#' @param features a vector of function names which return numeric vectors of features.
#' All features returned by these functions must be named if they return more than one feature.
#' Existing functions from installed packages may be used, but the package must be loaded first.
#' Functions must return a result for all time series, even if it is just NA.
#' @param scale if \code{TRUE}, time series are scaled to mean 0 and sd 1 before features
#' are computed.
#' @param trim if \code{TRUE}, time series are trimmed by \code{trim_amount} before features
#' are computed. Values larger than \code{trim_amount} in absolute value are set to \code{NA}.
#' @param trim_amount Default level of trimming if \code{trim==TRUE}.
#' @return A feature matrix (in the form of a tibble) with each row corresponding to
#' one time series from tslist, and each column being a feature.
#' @examples
#' mylist <- list(sunspot.year, WWWusage, AirPassengers, USAccDeaths)
#' tsfeatures(mylist)
#' @export
tsfeatures <- function(tslist,
                       features = c("frequency","stl_features","entropy","acf1"),
                       scale=TRUE, trim=FALSE, trim_amount=0.1)
{
  if(!is.list(tslist))
    tslist <- as.list(tslist)
  if(scale)
    tslist <- map(tslist, scalets)
  if(trim)
    tslist <- map(tslist, trimts, trim=trim_amount)

  # Compute all features
	flist <- funlist <- list()
  # Assuming that didn't generate an error, we will proceed
	for(i in seq_along(features))
  {
    flist[[i]] <- map(tslist,
      function(x){match.fun(features[i])(x)})
    # Check names
    if(is.null(names(flist[[i]][[1]])))
      flist[[i]] <- map(flist[[i]],
        function(x){
          names(x) <- features[i]
          return(x)})
  }

	# Unpack features into a list of numeric vectors
  featurelist <- list()
	for(i in seq_along(tslist))
    featurelist[[i]] <- unlist(map(flist, function(u)u[[i]]))

  # Find feature names
  featurenames <- map(featurelist, names)
  fnames <- unique(unlist(featurenames))
  if(any(featurenames==""))
    stop("Some unnamed features")

  # Create feature matrix
  fmat <- matrix(NA_real_, nrow=length(tslist), ncol=length(fnames))
  colnames(fmat) <- fnames
  rownames(fmat) <- names(tslist)
  for(i in seq_along(tslist))
    fmat[i,featurenames[[i]]] <- featurelist[[i]][featurenames[[i]]]

  return(tibble::as_tibble(fmat))
}

# Scale time series
scalets <- function(x)
{
  y <- as.ts(as.numeric(scale(x, center=TRUE, scale=TRUE)))
  tsp(y) <- tsp(x)
  return(y)
}

# Trim time series
trimts <- function(x, trim = 0.1)
{
  qtl <- quantile(x, c(trim, 1 - trim), na.rm = TRUE)
  x[x < qtl[1L] | x > qtl[2L]] <- NA
  return(x)
}
