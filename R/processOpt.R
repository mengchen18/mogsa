#' preprocessing of input data in \code{\link{mbpca}}.
#' 
#' An internal function called by \code{\link{mbpca}}.
#' 
#' 
#' @param x A list of matrices, rows are observations and columns are variables
#' @param center A logical variable indicates whether columns should be
#' centered
#' @param scale A logical variable indicates whether columns should be scaled
#' @param option A charater string could be one of c("lambda1", "inertia",
#' "uniform") to indicate how the different matrices should be normalized. If
#' "lambda1", the matrix is divided by its the first singular value, if
#' "inertia", the matrix is divided by its total inertia (sum of square), if
#' "uniform", none of them would be done.
#' @return A \code{list} of normalized matrix.
#' @export
#' @author Chen Meng
processOpt <-
function(x, center=TRUE, scale=FALSE, option = c("lambda1", "inertia", "uniform")) {
  
  opt <- match.arg(option)  
  
  if (is.null(names(x)))
    names(x) <- paste("data", 1:length(x), sep = "_")

  x <- lapply(x, scale, center, scale)
  if (opt == "lambda1") {
    w <- sapply(x, function(xx) 1/svd.solver(xx, nf = 1)$d)
  } else if (opt == "inertia") {
    w <- sapply(x, function(xx) 1/sqrt(sum(xx^2)))
  } else if (opt == "uniform") {
    w <- rep(1, length(x))
  }
  mapply(SIMPLIFY = FALSE, function(xx, ww) xx*ww, xx=x, ww=w)
}
