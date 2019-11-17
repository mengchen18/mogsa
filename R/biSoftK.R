#' NIPALS algorithm with soft thresholding operator on rows and columns
#' 
#' An internal function called by \code{\link{mbpca}}.
#' 
#' This function also use the NIPALS algorithm, but it generalized nipalsSoftK
#' from several aspects: 1. Allowing sparsity on both columns and rows of
#' matrices 2. Allowing weights for columns and rows 3. Allowing loading and/or
#' score vectors of blocks to be unit length 4. Allowing only positive number
#' in loading and score vectors
#' 
#' @param x The input matrix, rows are observations, columns are variables
#' @param maxiter Number of maximum interation the algorithm can run
#' @param kp The number (>=1) or proportion (<1) of variables want to keep.  It
#' could be a single value or a vector has the same length as x so the sparsity
#' of individual matrix could be different.
#' @param kt The number (>=1) or proportion (<1) of non-zero scores for
#' obvservations.
#' @param weight.p The weight of variables. It could be 1) a vector has the
#' same length as x, one value for each table/block; 2) one number, all
#' variables share the same weight or 3) a list of vectors, the length of each
#' vector should be the same with the columns numbers of the corresponding
#' table/block, so every variables has a unique weight.
#' @param weight.t The weight for observation. For accepted values or formats,
#' see weight.p.
#' @param pos Logical value, if only non-negaitve values in the loading and
#' score vectors.
#' @param unit.pb Logical value, whether the length of table/block loading
#' should be unit length.
#' @param unit.tb Logical value, whether the length of table/block score should
#' be unit length.
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
#' # examples
#' data("NCI60_4arrays")
#' d <- lapply(NCI60_4arrays[2:4], function(x) scale(t(x)))
#' # 
#' x1 <- biSoftK(d, maxiter = 1000, kp = Inf, kt = Inf,
#'               weight.p = rep(1, length(d)),
#'               weight.t = rep(1, length(d)),
#'               pos = FALSE,
#'               unit.pb = TRUE, unit.tb = TRUE)
#' 
#' x1s <- biSoftK(d, maxiter = 1000, kp = Inf, kt = Inf,
#'               weight.p = 1,
#'               weight.t = 1,
#'               pos = FALSE,
#'               unit.pb = TRUE, unit.tb = TRUE)
#' 
#' identical(x1, x1s)
#' barplot(c(x1$t), col= as.factor(substr(colnames(NCI60_4arrays$agilent), 1, 2)))
#' 
#' x2 <- biSoftK(d, maxiter = 1000, kp = 30, kt = Inf,
#'               weight.p = rep(1, length(d)),
#'               weight.t = rep(1, length(d)),
#'               pos = FALSE,
#'               unit.pb = TRUE, unit.tb = TRUE)
#' barplot(c(x2$t), col= as.factor(substr(colnames(NCI60_4arrays$agilent), 1, 2)))
#' 
#' x3 <- biSoftK(d, maxiter = 1000, kp = 40, kt = Inf,
#'               weight.p = rep(1, length(d)),
#'               weight.t = rep(1, length(d)),
#'               pos = TRUE,
#'               unit.pb = TRUE, unit.tb = TRUE)
#' barplot(c(x3$t), col= as.factor(substr(colnames(NCI60_4arrays$agilent), 1, 2)))
#' 
#' x4 <- biSoftK(d, maxiter = 1000, kp = 30, kt = 6,
#'               weight.p = rep(1, length(d)),
#'               weight.t = rep(1, length(d)),
#'               pos = FALSE,
#'               unit.pb = TRUE, unit.tb = TRUE)
#' barplot(c(x4$t), col= as.factor(substr(colnames(NCI60_4arrays$agilent), 1, 2)),
#'         names.arg = colnames(NCI60_4arrays$agilent), las = 2)
#' 
#' plot(x4$t, x4$tb$hgu133)
#' plot(x4$t, x4$tb$hgu133p2)
#' plot(x4$t, x4$tb$hgu95)
#' heatmap(t(d$hgu133[, which(x4$pb$hgu133 != 0) ]))
#' barplot(t(d$hgu133[, which.max(x4$pb$hgu133) ]), las = 2)
#' barplot(t(d$hgu133[, which.max(x4$pb$hgu133) ]), las = 2)

biSoftK <- function (x, maxiter, kp, kt, weight.p, weight.t, pos = FALSE, unit.pb = TRUE, unit.tb = FALSE) {
  
  if (length(kp) < length(x)) 
    kp <- rep(kp, length.out = length(x))
  if (length(kt) < length(x))
    kt <- rep(kt, length.out = length(x))
  
  regproj <- function(xb, t, kp, kt, wp, wt, unit.pb, unit.tb, pos) {
    pb <- t(xb) %*% t/c(t(t) %*% t)
    if (!unit.pb)
      pb <- normvec(pb)
    pb <- softK(pb, kp, w = wp, pos = pos)
    if (unit.pb)
      pb <- normvec(pb)
    tb <- xb %*% pb
    tb <- softK(tb, kt, w = wt, pos = pos)
    if (unit.tb)
      tb <- normvec(tb)
    list(tb = tb, pb = pb)
  }
  
  # t <- fast.svd(do.call("cbind", x))$u[, 1, drop = FALSE]
  t <- svd.solver(do.call("cbind", x), nf = 1)$u
  for (i in 1:maxiter) { 
    told <- t
    rp <- mapply(SIMPLIFY = FALSE, function(x, kp, kt, wp, wt, upb, utb, pos) {
      regproj(x, t, kp, kt, unit.pb = upb, unit.tb = unit.tb, wp = wp, wt = wt, pos = pos) ## stopped here
    }, x = x, kp = kp, kt = kt, wp = weight.p, wt = weight.t, pos = pos, upb = unit.pb, utb = unit.tb)
    
    tm <- sapply(rp, "[[", "tb")
    w <- t(tm) %*% t/c(t(t) %*% t)
    w <- w/sqrt(sum(w^2))
    t <- tm %*% w
    if (isTRUE(all.equal(c(t), c(told)))) 
      break
    if (i == maxiter) 
      cat("  Note: maximum number of iterations was reached, algrithm may not converge.\n")
  }
  res <- list(tb = lapply(rp, "[[", "tb"), pb = lapply(rp, "[[", "pb"), t = t, w = w)
  return(res)
}


