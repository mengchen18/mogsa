#' Extract global scores from an object of class \code{\link{moa-class}}.
#' 
#' Extract global scores from an object of class \code{\link{moa-class}}.
#' 
#' 
#' @param moa An object of class \code{\link{moa-class}}
#' @return A matrix of global score
#' @author Chen Meng
#' @seealso \code{\link{moaCoef}}
#' @examples
#' 
#' # see examples in \code{\link{mbpca}}
#' 
#' data("NCI60_4arrays")
#' moa <- mbpca(NCI60_4arrays, ncomp = 10, k = "all", method = "globalScore", option = "lambda1", 
#'              center=TRUE, scale=FALSE)
#' 
#' genes <- moaCoef(moa)
#' scr <- moaScore(moa)
#' 
#' 
moaScore <- function(moa) moa@fac.scr



#' Extract the loadings/coefficients from an object of class
#' \code{\link{moa-class}}.
#' 
#' Extract the loadings/coefficients from an object of class
#' \code{\link{moa-class}}.
#' 
#' 
#' @param moa An object of class \code{\link{moa-class}}.
#' @return It returns a list consist of two components:
#' 
#' \code{coefMat} - the loading matrix
#' 
#' \code{nonZeroCoef} - it is a \code{list} of \code{data.frame} to list the
#' non-zero coefficient variable in each of loading vectors and data sets. The
#' element names are in a format as
#' 
#' "xxxx.yy.zzz"
#' 
#' xxxx - are the data names, tells the data set where a varirable is from
#' 
#' yy - the number of Axes, for example, "V1" indicate the variable has a
#' non-zero coefficient in the first loading vector.
#' 
#' zzz - could be either "pos" (coefficient >0) or "neg" (coefficient < 0)
#' 
#' The \code{data.frame} has two columns, the first column is the ID of a
#' variable the second column is the coefficient/loading.
#' @author Chen Meng
#' @seealso \code{\link{moaScore}}
#' @examples
#' 
#' # see examples in \code{\link{mbpca}}
#' 
#' data("NCI60_4arrays")
#' moa <- mbpca(NCI60_4arrays, ncomp = 10, k = "all", method = "globalScore", option = "lambda1", 
#'              center=TRUE, scale=FALSE)
#' 
#' genes <- moaCoef(moa)
#' scr <- moaScore(moa)
#' 
#' 
moaCoef <- function(moa) {
  mm <- moa@loading
  m <- as.data.frame(mm)
  nd <- length(moa@data)
  fac <- rep(1:nd, moa@tab.dim[1, ])
  ms <- split(m, fac)
  
  r <- lapply(ms, function(x) {
    var <- rownames(x)
    lapply(x, function(y) {
      names(y) <- var
      neg <- sort(y[y<0], decreasing = FALSE)
      pos <- sort(y[y>0], decreasing = TRUE)
      
      neg <- data.frame(id=names(neg), coef=neg)
      pos <- data.frame(id=names(pos), coef=pos)
      list(neg=neg, pos=pos)
      })
  })
  names(r) <- names(moa@data)
  list(coefMat = mm, nonZeroCoef = unlist(lapply(r, unlist, recursive=FALSE), recursive = FALSE))
}

