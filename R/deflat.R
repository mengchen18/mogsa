#' deflat function used by \code{\link{mbpca}}
#' 
#' An internal function called by \code{\link{mbpca}}.
#' 
#' 
#' @param x A \code{list} of \code{matrix} want to deflat
#' @param t The global scores returned by \code{\link{msvd}} or
#' \code{\link{nipalsSoftK}}
#' @param tb The block scores returned by \code{\link{msvd}} or
#' \code{\link{nipalsSoftK}}
#' @param pb The block loadings returned by \code{\link{msvd}} or
#' \code{\link{nipalsSoftK}}
#' @param method A charater to specify the deflation strateg, could be one of
#' c("globalScore", "blockLoading", "blockScore").
#' @return A \code{list} of deflated \code{matrix}
#' @author Chen Meng
deflat <-
function(x, t, tb, pb, method="globalScore") {
  # globalScore, blockScore, blockLoading
  switch(method,
         "globalScore" = lapply(x, function(xb) { xb - t %*% t(t) %*% xb / c(t(t) %*% t) }),
         "blockLoading" = mapply(SIMPLIFY = FALSE, function(xb, pb) {xb - xb %*% pb %*% t(pb)}, x, pb),
         "blockScore" = mapply(SIMPLIFY = FALSE, function(xb, tb) {xb - tb %*% t(tb) %*% xb / c(t(tb) %*% tb)}, x, tb))
}
