#' Bootstrap mbpca to estimate the coherence of different data sets
#' 
#' Bootstrap mbpca to estimate the coherence of different data sets and
#' estimate the number of components should be included in an analysis.
#' 
#' Bootstrap method were used to determine the components that are presenting
#' significant concordant structure between datasets.
#' 
#' @param moa An object of \code{\link{moa}} returned by \code{\link{mbpca}}.
#' @param mc.cores Integer; number of cores used in bootstrap. This value is
#' passed to function \code{mclapply}
#' @param B Integer; number of bootstrap
#' @param replace Logical; sampling with or without replacement
#' @param resample Could be one of "sample", "gene" or "total". "sample" and
#' "gene" means sample-wise and variable-wise resampling, repectively. "total"
#' means total resampling.
#' @param log Could be "x", "y" or "xy" for plot log axis
#' @param ncomp Passed to function \code{\link{mbpca}}. In most of cases, user
#' do not need to specify this argument because it could be inferred from
#' \code{moa}.
#' @param method Passed to function \code{\link{mbpca}}.In most of cases, user
#' do not need to specify this argument because it could be inferred from
#' \code{moa}.
#' @param maxiter Passed to function \code{\link{mbpca}}.In most of cases, user
#' do not need to specify this argument because it could be inferred from
#' \code{moa}.
#' @param svd.solver Passed to function \code{\link{mbpca}}.In most of cases,
#' user do not need to specify this argument because it could be inferred from
#' \code{moa}.
#' @param plot Logical; whether the result should be plotted.
#' @return It returns a matrix, columns are eigenvalues for different
#' components. Each rows is a bootstramp sample.
#' @author Chen Meng
#' @export
#' @examples
#' 
#' # see examples in \code{\link{mbpca}}
#' 
bootMbpca <-
function(moa, mc.cores=1, B = 100, replace=TRUE, resample=c("sample", "gene", "total"), log="y",
                      ncomp = NULL, method = NULL, maxiter=1000, svd.solver=c("svd", "fast.svd", "propack"), plot=TRUE) {
  data <- moa@data
  call <- moa@call
  
  if (is.null(ncomp)) 
    ncomp <- call$ncomp
  if (is.null(method)) {
    method <- call$method
    cat(paste("method is set to '", method, "'.\n", sep = ""))
  }
  if (!is.null(call$maxiter) && call$maxiter > 0 && is.integer(call$maxiter)) {
    maxiter <- call$maxiter
    cat(paste("maxiter is set to ", maxiter, ".\n", sep = ""))
  }  
  if (is.null(call$k))
    call$k <- "all"
  if (call$k != "all") {
    call$k <- "all"
    call$verbose <- FALSE
    moa <- eval(call)
  }
  
  ncomp <- min(c(ncomp), length(moa@eig))
  svd.solver <- match.arg(svd.solver)
  resample <- match.arg(resample)
  btk <- bootMbpcaK(data, B = B, mc.cores=mc.cores, replace = replace, resample=resample,
                    option = "uniform", center = FALSE, scale = FALSE,
                    ncomp = ncomp, k = "all", method = method, maxiter=maxiter, svd.solver=svd.solver)
  
  if (plot) {
    sc <- min(ncol(btk), length(moa@eig))
    isc <- 1:sc
    
    boxplot(rbind(moa@eig[isc], btk[, isc]), col=NA, border = NA, log=log)
    
    sds <- apply(btk, 2, sd)
    means <- apply(btk, 2, mean)
    points(isc,  means, pch=15, col="red")
    points(isc,  means+1.96*sds, col="red", pch="_")
    points(isc,  means-1.96*sds, col="red", pch="_")
    segments(isc, means+1.96*sds, isc, means-1.96*sds, col="red")
    
    lines(isc, moa@eig[isc], pch=20)
    points(isc, moa@eig[isc], pch=20)
  }
  
  colnames(btk) <- paste("PC", 1:ncol(btk), sep="")
  rownames(btk) <- paste("sample", 1:nrow(btk), sep="")
  btk
}
