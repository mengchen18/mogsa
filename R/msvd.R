#' SVD based algorithm to calculate block Score and global scores for
#' \code{\link{mbpca}}.
#' 
#' An internal function called by \code{\link{mbpca}}. It returns the result
#' comparable with nipalsSoftK, but way faster since it uses the SVD algorithm.
#' No sparse opertors in this function.
#' 
#' 
#' @param x The input matrix, rows are observations, columns are variables
#' @param svd.sol A function object to specify the preferred SVD solver,
#' default is \code{svd}.
#' @return an \code{list} object contains the following elements:
#' 
#' \code{tb} - the block scores
#' 
#' \code{pb} - the block loadings
#' 
#' \code{t} - the global scores
#' 
#' \code{w} - the wegihts of block scores to construct the global scor
#' @author Chen Meng
#' @seealso \code{\link{nipalsSoftK}}
msvd <-
function(x, svd.sol=svd) {
  nd <- length(x)
  nVar <- sapply(x, ncol)
  if (!is.null(names(x)))
    idx <- rep(names(x), times = nVar) else
      idx <- rep(1:length(x), times = nVar)
  nm <- lapply(x, colnames)
  
  tab <- do.call("cbind", x)
  
  dc <- svd.sol(tab)
  
  res <- list()
  res$t <- dc$u[, 1, drop=FALSE] * dc$d[1]
  pb <- split(dc$v[, 1, drop=FALSE], idx)
  for (i in names(x)) names(pb[[i]]) <- colnames(x[[i]]) 
  res$pb <- lapply(pb, function(x) as.matrix(x/sqrt(sum(x^2))))
  
  res$tb <- mapply(SIMPLIFY = FALSE, function(m, v) {
    m %*% v
  }, m=x, v=res$pb[names(x)])
  
  tm <- do.call("cbind", res$tb)
  res$w <- t(tm) %*% res$t / c(t(res$t) %*% res$t)
  rownames(res$w) <- names(x)
  
  res <- res[c("tb", "pb", "t", "w")]
  res$tb <- res$tb[names(x)]
  res$pb <- res$pb[names(x)]
  return(res)
}
