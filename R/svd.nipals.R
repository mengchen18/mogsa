#' singular valude decomposition using NIPALS algorithm
#' @description The function may be faster when x is a sparseMatrix object and 
#'   a few components need to be calculated. In addition, NA values could be 
#'   presented in x.
#' @param x a numeric matrix whose SVD decomposition is to be computed  
#' @param nf the number of singular vectors want to computed
#' @param ... other parameters passed to \code{\link[nipals]{nipals}}, except
#'   the following argument:
#'   center (=FALSE), scale (=FALSE), gramschmidt (=FALSE)#' 
#' @import Matrix
#' @import nipals
#' @seealso \code{\link{svd}}
#' @return the same as svd, list of three components, d, u, v
#' @export
#' @examples 
#' # small matrix
#' m <- matrix(rnorm(200), 20, 10)
#' s <- svd.nipals(m)
#' 
#' # large sparse matrix
#' \dontrun{
#' library(Matrix)
#' m <- matrix(sample(c(rnorm(10000), rep(0, 9990000))), 10000, 1000)
#' ms <- Matrix(m, sparse=TRUE)
#' system.time(
#'   s1 <- svd(m)
#' )
#' system.time(
#'   s2 <- svd.nipals(ms, nf = 1)
#' ) 
#' }
#' 
svd.nipals <- function(x, nf, ... ) {
  if (missing(nf))
    nf <- min(dim(x))
  decomp <- nipals::nipals(x, ncomp = nf, scale = FALSE, 
                   center = FALSE, gramschmidt = FALSE, ...)
  r <- list()
  r$d <- decomp$eig
  r$u <- decomp$scores
  r$v <- decomp$loadings
  r
}


#' Singular Value Decomposition of a Matrix
#' @description different methods are automatically selected to calculate svd
#'   according to matrix properties, such as missing values, object class, size. 
#' @param x a numeric or complex matrix whose SVD decomposition is to be computed.
#' @param nf the number of singular vectors want to computed. This will be ignored 
#'   if fast.svd or svd function is used. 
#' @param opts.svd.nipals A list of named parameters passed to \code{\link{svd.nipals}}
#' @param opts.svds A list of named parameters passed to \code{\link[RSpectra]{svds}}
#' @param opts.fast.svd A list of named parameters passed to \code{\link[corpcor]{fast.svd}}
#' @param opts.svd A list of named parameters passed to \code{\link[base]{svd}}
#' @details There are 4 different options could be used:
#'   1. \code{\link{svd.nipals}} When there are missing values in the matrix
#'   2. \code{\link[RSpectra]{svds}} when x is an object of class \code{Matrix} or only a small 
#'     number of singular vectors to be calculated (<3)
#'   3. \code{\link[corpcor]{fast.svd}} when x is big fat or thin matrix
#'   4. \code{\link[base]{svd}} if not any other cases
#' @author Chen Meng
#' @import corpcor
#' @import Matrix
#' @import RSpectra
#' @export
#' @examples 
#' m <- matrix(sample(c(rnorm(1000), rep(0, 999000))), 10000, 100)
#' decomp <- svd.solver(m)
#' attr(decomp, "solver")
#' 
#' decomp <- svd.solver(m[1:200, ])
#' attr(decomp, "solver")
#' \dontrun{
#' library(Matrix) 
#' ms <- Matrix(m, sparse=TRUE)
#' decomp <- svd.solver(ms, nf = 2)
#' attr(decomp, "solver")
#' 
#' mm <- Matrix(m)
#' decomp <- svd.solver(mm, nf = 2)
#' attr(decomp, "solver")
#' 
#' mna <- m
#' mna[sample(1:length(mna), size = 100)] <- NA
#' decomp <- svd.solver(mna, nf = 2)
#' attr(decomp, "solver")
#' 
#' mnas <- Matrix(mna, sparse = TRUE)
#' decomp <- svd.solver(mnas, nf = 2)
#' attr(decomp, "solver")
#' }

svd.solver <- function(
  x, nf, opts.svd.nipals=list(), opts.svds=list(), opts.fast.svd = list(), opts.svd=list()
  ) {
  if (miss.nf <- missing(nf))
    nf <- min(dim(x))
  
  smallnf <- nf < 3
  
  rt <- nrow(x)/ncol(x)
  if (any(is.na(x))) {
    r <- do.call(svd.nipals, c(list(x=x, nf = nf), opts.svd.nipals))
    solver <- "svd.nipals"
    # cat("svd.nipals used. \n" )
  } else if (inherits(x, "Matrix") || smallnf) {
    r <- do.call(svds, c(list(A=x, k= nf), opts.svds))
    solver <- "svds"
    # cat("svds used. \n" )
  } else if (length(x) >= 1e6 & abs(log10(rt)) >= 1) {
    if (!miss.nf)
      message("function 'fast.svd' used, all singular vectors will be computed")
    r <- do.call(fast.svd, c(list(m=x), opts.fast.svd))
    solver <- "fast.svd"
    # cat("fast.svd used. \n" )
  } else {
    if (!miss.nf)
      message("Function 'svd' used, all singular vectors will be computed")
    r <- do.call(svd, c(list(x=x), opts.svd))
    solver <- "svd"
    # cat("svd used. \n" )
  }
  if (inherits(x, "Matrix")) {
    r$u <- Matrix(r$u)
    r$v <- Matrix(r$v)
  }  
  attr(r, "solver") <- solver
  r
}

c <- function(x, ...) {
  if (inherits(x, "Matrix"))
    as.numeric(x) else
      base::c(x, ...)
}