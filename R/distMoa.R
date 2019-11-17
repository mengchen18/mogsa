#' Calculate the distance matrix from an object of class
#' \code{\link{moa-class}}.
#' 
#' A convenient function to calculate the distance matrix from an object of
#' class \code{\link{moa-class}}.
#' 
#' 
#' @param x An object of class \code{\link{moa-class}}.
#' @param nf Integer; the number of component used to calculate the distance.
#' Default setting (NA) will keep all the axes.
#' @param tol Numerical; the tolerance of component with low variance.
#' @param method passed to function \code{dist}
#' @param diag passed to function \code{dist}
#' @param upper passed to function \code{dist}
#' @param p passed to function \code{dist}
#' @return An object of class \code{dist}, see function "dist".
#' @author Chen Meng
#' @examples
#' 
#' # see examples in \code{\link{mbpca}}
#' 
#' data("NCI60_4arrays")
#' moa <- mbpca(NCI60_4arrays, ncomp = 10, k = "all", method = "globalScore", option = "lambda1", 
#'              center=TRUE, scale=FALSE)
#' 
#' dst <- distMoa(moa)
#' 
#' 
#' 
distMoa <-
function(x, nf=NA, tol=1e-5, method = "euclidean", diag = FALSE, upper = FALSE, p = 2) {
  if (is.na(nf))
    nf <- Inf
  if (inherits(x, "moa")) {
    nfi <- x@eig > tol
    x <- moaScore(x)
  }
  if (nf > ncol(x) | nf > sum(nfi)) {
    nf <- min(ncol(x), sum(nfi))
    cat(paste("nf set to ", nf, ".\n", sep = ""))
  }
  nfi[-(1:nf)] <- FALSE
  x <- x[, nfi]
  dist(x, method = method, diag = diag, upper = upper, p = p)
}
