#' compute the power of a matrix
#' 
#' the power of a matrix
#' 
#' The power of a matrix is calculated in two steps: decompostion step: x=UDV'
#' and the reconstruction step: x^n=U*D^n*V' In the reconstruction, the
#' singular vectors with a singular value more than tol are kept.
#' 
#' @param x a numerical matrix object that the power of which should be
#' calculated
#' @param n The matrix to the power of
#' @param nf The number of axes kept in the calculation of SVD and
#' reconstruction
#' @param tol The tolerance of the axis, singular vectors with singular value
#' lower than tol will be ignored in the reconstruction.
#' @return A matrix x^n
#' @note Called by the wsvd function.
#' @author Chen Meng
#' @seealso See Also \code{\link{wsvd}}
#' @export
#' @examples
#' 
#'   set.seed(56)
#'   m <- matrix(rnorm(15), 5, 3)
#'   s <- matpower(m, 2)
#'   s <- matpower(m, -2)
#' 
matpower <- function(x, n, nf=min(dim(x)), tol=1e-7) {
  
  if (inherits(x, "matrix")) {
    s <- svd(x, nu = nf, nv = nf)
    m <- sum(s$d > tol)
    nf <- min(nf, m)
    x <- s$u[, 1:nf] %*% diag(s$d[1:nf]^n) %*% t(s$v[, 1:nf])
  } else 
    stop("x need to be a matrix object!")
  return(x)
}

