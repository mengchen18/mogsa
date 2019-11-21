#' NIPALS algorithm with soft thresholding operator
#' 
#' An internal function called by \code{\link{mbpca}}.
#' 
#' @param x The input matrix, rows are observations, columns are variables
#' @param maxiter # of maximum interation the algorithm can run
#' @param k The number (>=1) or proportion (<1) of variables want to keep.  It
#' could be a single value or a vector has the same length as x so the sparsity
#' of individual matrix could be different.
#' @return an \code{list} object contains the following elements:
#' 
#' \code{tb} - the block scores
#' 
#' \code{pb} - the block loadings
#' 
#' \code{t} - the global scores
#' 
#' \code{w} - the wegihts of block scores to construct the global score.
#' @author Chen Meng
#' @export
#' @seealso \code{\link{msvd}}
nipalsSoftK <-
function(x, maxiter, k) {
  
  if (length(k) < length(x))
    k <- rep(k, length.out = length(x))
  
  # t <- svd(do.call("cbind", x))$u[, 1]
  t <- svd.solver(do.call("cbind", x), nf = 1)$u
  
  regproj <- function(xb, t, k) { 
    pb <- t(xb) %*% t / c(t(t) %*% t)
    pb <- pb/sqrt(sum(pb^2))
    pb <- softK(pb, k)  # soft-thresholding
    tb <- xb %*% pb
    list(tb=tb, pb=pb) # t-score, p-loading
  }
  
  for (i in 1:maxiter) {
    told <- t
    rp <- mapply(SIMPLIFY = FALSE, function(x, y) regproj(x, t, y), x=x, y=k)
    tm <- sapply(rp, "[[", "tb")
    if (is.list(tm))
      tm <- do.call(cbind, tm)
    w <- t(tm) %*% t / c(t(t) %*% t)
    w <- w/sqrt(sum(w^2))
    #  w <- w/sum(w)
    t <- tm %*% w
    if (all.equal(c(t), c(told))[1] == TRUE)
      break
    if (i == maxiter)
      cat("  Note: maximum number of iterations was reached, algrithm may not converge.\n")
  }
  
  res <- list(tb = lapply(rp, "[[", "tb"),
              pb = lapply(rp, "[[", "pb"),
              t = t,
              w = w)
  
  return(res)
}


