#' Weighted singular value decomposition (SVD)
#' 
#' The weighted version of singular value decomposition.
#' 
#' The weighted version of generalized singular value decomposition (SVD) of
#' matrix A = UDV' with the constraints U'D1U = I and V'D2V = I D1 and D2 are
#' two matrices express constraints imposed on the rows and the columns of
#' matrix A.
#' 
#' @param X A numeric matrix whose wSVD decomposition is to be computed.
#' @param D1 A square matrix or vector. The left constraint/weight matrix
#' (symmetric and positive in diagonal). The dimension of D1 should be the same
#' with the number of rows in X. A vector input will be converted to a diagnal
#' matrix.
#' @param D2 A square matrix or vector. The right constraint/weight matrix
#' (symmetric, positive in diagonal). The dimension of D1 should be the same
#' with the number of columns in X. A vector input will be converted to a
#' diagnal matrix.
#' @return d - singular values
#' 
#' u - left singular vectors
#' 
#' v - right singular vectors
#' 
#' D1 - the left weight matrix (directly from input)
#' 
#' D2 - the right weight matrix (directly from input)
#' @author Chen Meng
#' @seealso svd
#' @references Herve Abdi. Singular Value Decomposition (SVD) and Generalized
#' Singular Value Decomposition (GSVD)
#' http://www.utdallas.edu/~herve/Abdi-SVD2007-pretty.pdf
#' @keywords SVD generalized SVD weighted SVD
#' @examples
#' 
#'     set.seed(56)
#'     m <- matrix(rnorm(15), 5, 3)
#'     wl <- rnorm(5)
#'     wr <- runif(3)
#'     s <- wsvd(X=m, D1=wl, D2=wr)
#'     # t(s$u) %*% diag(wl) %*% s$u
#'     # t(s$v) %*% diag(wr) %*% s$v
#'     # all.equal(m, as.matrix(s$u) %*% diag(s$d) %*% t(s$v))
#' 
wsvd <- function(X, D1=diag(1, nrow(X)), D2=diag(1, ncol(X))) {
  
  if (is.vector(D1))
    D1 <- diag(D1)
  if (is.vector(D2))
    D2 <- diag(D2)
  i1 <- identical(D1, t(D1))
  i2 <- identical(D2, t(D2))
  if (!(i1 & i2))
    warning("non-symetric distance matrix")
  
  r <- list(D1=D1, D2=D2)
  X <- matpower(D1, 1/2) %*% X
  X <- X %*% matpower(D2, 1/2)
  result <- svd(X)
  r$d <- result$d
  r$u <- matpower(D1, -1/2) %*% result$u
  r$v <- matpower(D2, -1/2) %*% result$v
  
  if (all(r$u[, 1] < 0)) {
    r$u <- r$u * (-1)
    r$v <- r$v * (-1)
  }
  return(r[c("d", "u", "v", "D1", "D2")])
}
