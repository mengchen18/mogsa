#' multiple omics data integration and gene set analysis
#' 
#' The main function called by users, omics data analysis and gene set
#' annotation.  A wrapper function of \code{\link{moa}} and
#' \code{\link{sup.moa}}.
#' 
#' A wrapper function of \code{\link{moa}} and \code{\link{sup.moa}}.
#' 
#' @param x An object of class \code{list} or \code{\link{moa-class}}. A list
#' would be a list of data frame.
#' @param sup An object of class \code{list} or \code{\link{moa.sup-class}}. A
#' list would be a list of supplementary data.
#' @param nf The number of principal components used to reconstruct, only used
#' when x is a an object of \code{list}.
#' @param factors The index of principal components used in the projection,
#' used when non-consecutive PC to be included in the analysis.
#' @param proc.row Preprocessing of rows. If x is a object of \code{list}, it
#' is passed \code{moa}
#' @param w.data Weights of datasets. If x is a object of \code{list}, it is
#' passed \code{moa}
#' @param w.row Weight of row. If x is a object of \code{list}, it is passed
#' \code{moa}
#' @param statis A logical indicates if statis algrithm should be used. If x is
#' a object of \code{list}, it is passed \code{moa}
#' @param ks.stat The logical indicates if the p-value should be calculated
#' using K-S statistic (the method used in "ssgsea" in GSVA package).  Default
#' is FALSE, which means using the z-score method. See \code{\link{sup.moa}}.
#' @param ks.B An integer to indicate the number of bootstrapping samples to
#' calculated the p-value of KS statistic.
#' @param ks.cores An integer indicate the number of cores to be used in
#' bootstrapping. It is passed to function \code{mclapply} in the
#' \code{parallel} package.
#' @param p.adjust.method The method of p value adjustment, passed to
#' \code{p.adjust} function.
#' @return An object of class \code{\link{mgsa-class}}.
#' @note This function will be changed to a generic function for "S4-style"
#' programming.
#' @author Chen Meng
#' @seealso \code{\link{moa}} and \code{\link{sup.moa}}
#' @references Preprint: Meng, C., Kuster, B., Peters, B., Culhane, AC.,
#' Moghaddas Gholami, A., moGSA: integrative single sample gene-set analysis of
#' multiple omics data. doi: http://dx.doi.org/10.1101/046904 Haenzelmann, S.,
#' Castelo, R. and Guinney, J. GSVA: Gene set variation analysis for microarray
#' and RNA-Seq data. BMC Bioinformatics, 14:7, 2013.  Barbie, D.A. et al.
#' Systematic RNA interference reveals that oncogenic KRAS-driven cancers
#' require TBK1. Nature, 462(5):108-112, 2009.
#' @keywords MVA supplementary data projection
#' @export
#' @examples
#' 
#'   # library(mogsa)
#'   # loading gene expression data and supplementary data
#'   data(NCI60_4array_supdata)
#'   data(NCI60_4arrays)
#' 
#'   # using a list of data.frame as input
#'   mgsa1 <- mogsa(x = NCI60_4arrays, sup=NCI60_4array_supdata, nf=9,
#'                  proc.row = "center_ssq1", w.data = "inertia", statis = TRUE)
#'   mgsa1x <- mogsa(x = NCI60_4arrays, sup=NCI60_4array_supdata, factors = c(1,3,6),
#'                  proc.row = "center_ssq1", w.data = "inertia", statis = TRUE)
#'   # using moa as input
#'   ana <- moa(NCI60_4arrays, proc.row = "center_ssq1", w.data = "inertia", statis = TRUE)
#'   smoa <- sup.moa(ana, sup=NCI60_4array_supdata, nf=3)
#'   mgsa2 <- mogsa(x = ana, sup=NCI60_4array_supdata, nf=9)
#'   mgsa3 <- mogsa(x = ana, sup=smoa)
#' 
mogsa <- function(x, sup, nf=NULL, factors = NULL, proc.row=NULL, w.data=NULL, 
  w.row=NULL, statis=FALSE, ks.stat=FALSE, ks.B = 1000, ks.cores = NULL, p.adjust.method = "fdr") {
  
  # if sup is NULL .....
  # extract data and moa
  if (inherits(x, "list")) {
    if (is.null(nf) & is.null(factors))
      stop("x is an object of \"list\", nf or factors need to be set.")
    if (is.null(proc.row))
      stop("x is an object of \"list\", proc.row need to be set.")
    if (is.null(w.data))
      stop("x is an object of \"list\", w.data need to be set.")
    if (is.null(statis))
      stop("x is an object of \"list\", statis need to be set.")
    data <- x
    r <- moa(data=data, proc.row=proc.row, 
      w.data=w.data, w.row=w.row, statis=statis)
    
    # sup data
    if (inherits(sup, "list")) {
      supr <- sup.moa (X=r, sup=sup, nf=nf, factors = factors, 
        ks.stat=ks.stat, ks.B = ks.B, ks.cores = ks.cores, p.adjust.method = p.adjust.method)
    } else if (inherits(sup, "moa.sup")) {
      stop("sup cannot be an object of class moa.sup if x is an object of class list.")
    } 
  } else if (inherits(x, "moa")) {
    if (!is.null(proc.row))
      cat("x is an object of \"moa\", proc.row is not used")
    if (!is.null(w.data))
      cat("x is an object of \"moa\", w.data is not used")
    if (!is.null(w.row))
      cat("x is an object of \"moa\", w.row is not used")
    if (!is.null(statis))
      cat("x is an object of \"moa\", statis is not used")
    data <- x@data
    r <- x

    # sup data
    if (inherits(sup, "list")) {
      supr <- sup.moa (X=r, sup=sup, nf=nf, factors = factors, 
        ks.stat=ks.stat, ks.B = ks.B, ks.cores = ks.cores, p.adjust.method = p.adjust.method)
    } else if (inherits(sup, "moa.sup")) {
      if (!is.null(nf))
        cat("x is an object of \"moa\" and sup is an object of class \"sup.moa\", nf is not used")
      supr <- sup
    }
  }
  
  new("mgsa",
      call=match.call(),
      moa=r,
      sup=supr)
}
