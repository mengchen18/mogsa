# singular valude decomposition using NIPALS algorithm
#' @description The function may be faster when x is a sparseMatrix object and 
#'   a few components need to be calculated. In addition, NA values could be 
#'   presented in x.
#' @param x a numeric matrix whose SVD decomposition is to be computed  
#' @param nf the number of components want to computed
#' @param ... other parameters passed to \code{\link[nipals]{nipals}}, except
#'   the following argument:
#'   center (=FALSE), scale (=FALSE), gramschmidt (=FALSE)#' 
#' @import nipals
#' @seealso \code{\link{svd}}
#' @return the same as svd, list of three components, d, u, v
#' @examples 
#' # library(corpcor)
#' m <- matrix(sample(c(rnorm(10000), rep(0, 9990000))), 10000, 1000)
#' ms <- Matrix(m, sparse=TRUE)
#' system.time(
#'   s1 <- svd(m)
#' )
#' system.time(
#'   s2 <- svd.niplas(ms, nf = 1)
#' ) 


svd.niplas <- function(x, nf, ... ) {
  decomp <- nipals(x, ncomp = nf, scale = FALSE, 
                   center = FALSE, gramschmidt = FALSE, ...)
  r <- list()
  r$d <- decomp$eig
  r$u <- decomp$scores
  r$v <- decomp$loadings
  r
}
