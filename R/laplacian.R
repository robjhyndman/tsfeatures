#' Finding embedding based on distances
#'
#' Various methods for embedding a distance matrix into a low-dimensional matrix
#'
#' @param distances an object of class "dist" (essentially the lower triangle of distances matrix).
#' @param k  embedding dimension
#' @param h  bandwidth for computing the similarity matrix. Only used for Laplacian methods,
#' apart from LaplacianMDS where hs is set to a large h.
#' @param ... any other arguments are ignored.
#'
#' @examples
#' mylist <- list(sunspot.year, WWWusage, AirPassengers, USAccDeaths)
#' fmat <- tsfeatures(mylist)
#' z <- embedding(dist(fmat))
#' plot(z)
#' @export

embedding <- function(
  distances,
  k = 2,
  method = c("Laplacian","Lrw","Lsym","Lsym2",
            "MDS","MDSiso","monoMDS","DPM","Rtsne"),
  h = median(distances), ...)
{
  method <- match.arg(method)

  if(!("dist" %in% class(distances)))
    stop("distances must be a matrix of distances")

  if(method == "Laplacian")
  {
    #Laplacian eigenmap with regular h
    w <- similarity(distances, h=h, ...)
    n <- NROW(w)
    D <- diag(rowSums(w))
    ei <- geigen::geigen(D-w,D,symmetric=TRUE)
    ei <- ei$vectors[,2:(k+1),drop=FALSE]
  }
  else if(method=="Lrw")
  {
    #Laplacian eigenmap. normalized Lrw
    w <- similarity(distances, h=h, ...)
    n <- NROW(w)
    D <- diag(1/rowSums(w))
    ei <- geigen::geigen(diag(n) - D %*% w,D,symmetric=TRUE)
    ei <- ei$vectors[,2:(k+1),drop=FALSE]
  }
  else if(method=="Lsym")
  {
    #Laplacian eigenmap. normalized Lsym
    w <- similarity(distances, h=h, ...)
    n <- NROW(w)
    wden <- rowSums(w)
    D <- diag(1/wden)
    Dhalf <- diag(1/sqrt(wden))
    ei <- geigen::geigen(diag(n) - Dhalf %*% w %*% Dhalf,D,symmetric=TRUE)
    ei <- ei$vectors[,2:(k+1),drop=FALSE]
  }
  else if(method=="Lsym2")
  {
    #Laplacian eigenmap. normalized Lsym
    w <- similarity(distances, h=h, ...)
    n <- NROW(w)
    wden <- rowSums(w)
    D <- diag(1/wden)
    Dhalf <- diag(1/sqrt(wden))
    ei <- eigen(Dhalf %*% w %*% Dhalf,symmetric=TRUE)
    ei <- ei$vectors[,1:k,drop=FALSE]
  }
  else if(method=="Rtsne")
  {
    ei <- Rtsne::Rtsne(distances, dims=k, perplexity=9)$Y
  }
  else if(method=="MDS")
    ei <- cmdscale(distances, k=k)
  else if(method=="MDSiso")
  {
    # Multidimensional scaling
    mds <- MASS::isoMDS(distances, k=k)
    ei <- mds$points
  }
  else if(method=="monoMDS")
  {
    # Multidimensional scaling
    mds <- vegan::monoMDS(distances, k=k, model="local")
    ei <- mds$points
  }
  else if(method=="DPM")
  {
    # Density preserving map
    ei <- dpm(distances, h=h, m=k, ...)
    colnames(ei) <- paste("Comp",1:k, sep="")
    rownames(ei) <- attr(distances, "Labels")
    return(structure(scale(ei),class="embedding"))
  }

  else
    stop("Not implemented")

  colnames(ei) <- paste("Comp",1:k, sep="")
  rownames(ei) <- attr(distances,"Labels")

  # Scale embedding
  ei <- scale(ei)
  # Then take signs so medians are positive
  # Only purpose of this is to avoid arbitrary changes in sign for a component
  med_ei <- apply(ei, 2, median)
  j <- med_ei < 0
  ei[,j] <- -ei[,j]

  return(tibble::as_tibble(ei))
}

# Compute similarity matrix based on pairwise distances
# 3 different kernels are possible
# The h parameter is selected to be very large if the argument is omitted

similarity <- function(distances, h,
              kernel=c("Gaussian","Epanechnikov","Triangle"))
{
  if(class(distances) != "dist")
    stop("distances should be of class 'dist'")
  kernel <- match.arg(kernel)

  if(missing(h))
    h <- 1000*max(distances)

  distances <- as.matrix(distances)
  if(kernel=="Gaussian")
    sims <- exp(-(distances/h)^2)
  else if(kernel=="Epanechnikov")
    sims <- pmax(1-(distances/h)^2, 0)
  else
    sims <- pmax(1-distances/h, 0)

  return(sims)
}



# Compute similarity matrix based on pairwise distances
# Only implements Epanechnikov kernel
# distances must be a dist class.
# returns sparse object
sparsesimilarity <- function(distances, h)
{
  if(!is.element("dist",class(distances)))
    stop("distances must be of class 'dist'")

  if(missing(h))
    h <- 1000*max(distances)

  n <- attr(distances,"Size")
  k <- distances < h
  col <- rep(1:(n-1),(n-1):1)[k]
  row <- matrix(rep(1:n,n), n,n)[lower.tri(matrix(0,n,n))]
  row <- row[k]
  v <- 1-(distances[k]/h)^2
  sims <- Rcsdp::simple_triplet_sym_matrix(row,col,v,n=n)

  return(sims)
}

# Compute row means of sparse symmetric matrix
rowMeansSparse <- function(x)
{
  result <- numeric(x$n)
  for(i in seq(x$n))
  {
    k <- (x$i==i | x$j==i)
    result[i] <- sum(x$v[k])
  }
  return(result/x$n)
}

# Return number of non-zeros in each row of sparse symmetric matrix
rowNonZero <- function(x)
{
  result <- numeric(x$n)
  for(i in seq(x$n))
  {
    k <- (x$i==i | x$j==i)
    result[i] <- sum(k)
  }
  return(result)
}

# Add two sparse symmetric matrices of equal size
addSparse <- function(x,y)
{
  if(x$n != y$n)
    stop("Matrices not the same size")

  # Combine rows, colums and non-zero elements
  i <- c(x$i,y$i)
  j <- c(x$j,y$j)
  v <- c(x$v,y$v)
  # Find duplicates
  z <- duplicated(cbind(i,j))
  if(any(z))
  {
    #Split duplicates into separate vectors
    i2 <- i[z]
    j2 <- j[z]
    v2 <- v[z]
    i <- i[!z]
    j <- j[!z]
    v <- v[!z]
    # Add together any duplicate values
    for(k in seq_along(i2))
    {
      l <- which(i==i2[k] & j==j2[k])
      v[l] <- v[l] + v2[k]
    }
  }
  return(Rcsdp::simple_triplet_sym_matrix(i,j,v,n=x$n))
}


dij <- function(i,j,n)
{
  Rcsdp::simple_triplet_sym_matrix(
    i=c(i,j,i),
    j=c(i,j,j),
    v=c(1,1,-1),n=n)
}

dpm <- function(d, h, m)
{
  n <- attr(d,"Size")
  w <- sparsesimilarity(d, h=h)
  f <- rowMeansSparse(w) - rowNonZero(w)/n
  b <- c(f[f<0],0)

  A <- list()
  nA <- 0
  for(i in seq(n))
  {
    k <- (w$i==i | w$j==i)
    if(any(k))
    {
      idx <- which(k)
      i0 <- w$i[idx]
      j0 <- w$j[idx]
      tmp <- Rcsdp::.simple_triplet_zero_sym_matrix(w$n)
      if(length(idx)>0)
      {
        for(j in seq_along(idx))
          tmp <- addSparse(tmp, dij(i0[j],j0[j],n))
      }
      tmp$v <- -tmp$v/h^2/n
      nA <- nA+1
      A[[nA]] <- list(tmp)
    }
  }
  A[[nA+1]] <- list(matrix(1,n,n))
  C <- list(Rcsdp::.simple_triplet_diag_sym_matrix(1,n))
  tmp <- Rcsdp::csdp(C, A, b, list(type="s", size=n))

  if(tmp$status!=0)
    warning("Not converged")
  return(tmp$X[[1]][,1:m,drop=FALSE])
}

