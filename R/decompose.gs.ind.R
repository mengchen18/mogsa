#' Data-wise or PC-wise decomposition of gene set scores for a single
#' observation.
#' 
#' Barplot of decomposed gene set scores, either with respect to datasets or
#' axes.
#' 
#' type=1 (the data-pc mode), the axes/PCs are represented as the narrow bars
#' with different colors and the background wide bars behind narrow bars are
#' gene set scores for datasets, which is calculated from the sum of all
#' underlying individual axes/PC scores.  When type=2 (the pc-data mode) the
#' interpreation of narrow and wide bars are in the other way around. If
#' type=3, both are shown.
#' 
#' This function could only be used to check the decomposition of gene set
#' scores of a single observation. So the function is not efficent when the
#' number of observation is large. Another function
#' \code{\link{decompose.gs.group}}, could be used in this case, particularly
#' when the cluster information of the observation panel is available.
#' 
#' @param x An object of class \code{\link{mgsa-class}} or
#' \code{\link{moa.sup-class}}
#' @param gs The gene set want to exam.
#' @param obs The observations want to exam.
#' @param type Which type of plot.  type=1 - the data-pc mode; type=2 - the
#' pc-data mode; type=3 - both. See detail.
#' @param nf The number of axes/PCs to be calculated and plotted.
#' @param plot A logical indicates if a plot should be drawn
#' @param col.data The bar color of datasets
#' @param col.pc The bar color of PCs
#' @param legend A logical if legend should be shown
#' @return Return nothing or a matrix depends on how argument \code{plot} is
#' set.
#' @author Chen Meng
#' @seealso See Also as \code{\link{decompose.gs.group}}
#' @export
#' @references TBA
#' @examples
#' 
#'   # library(mogsa)
#'   # loading gene expression data and supplementary data
#'   data(NCI60_4array_supdata)
#'   data(NCI60_4arrays)
#'   mgsa <- mogsa(x = NCI60_4arrays, sup=NCI60_4array_supdata, nf=9,
#'                 proc.row = "center_ssq1", w.data = "inertia", statis = TRUE)
#' 
#'   allgs <- colnames(NCI60_4array_supdata[[1]])
#'   # plot
#'   decompose.gs.ind(x=mgsa, gs=allgs[5], obs="BR.MDA_MB_231", type=2, nf=5)
#'   # or
#'   decompose.gs.ind(x=getmgsa(mgsa, "sup"), gs=allgs[5], obs="BR.MDA_MB_231", type=3, nf=5)
#' 
decompose.gs.ind <- function(x, gs, obs, type=3, nf=2, plot=TRUE,
                   col.data=NULL, col.pc=NULL, legend=TRUE) {
  
  # type = 1 - the data - pc mode
  # type = 2 - the PC -data mode
  # type = 3 - both
  
  if (inherits(x, "mgsa"))
    x <- x@sup
  scl <- lapply(x@score.sep, function (x) x[1:nf])
  m <- sapply(scl, function(x) 
    sapply(x, function(y) mean(y[gs, obs]))) # sum changed here 
  nc <- ncol(m)
  nr <- nrow(m)
  
  if (plot) {
    lgt1 <- lgt2 <- NULL
    if (legend) {
      lgt1 <- names(x@score.pc)[1:nf]
      lgt2 <- names(x@score.data)
    }
    
    if (type == 3)
      layout(matrix(1:2, 1, 2))
    if (type %in%  c(1, 3)) { # data-pc
      ma.d <- max(c(m, colSums(m)))
      mi.d <- min(c(m, colSums(m)))
      barplot(colSums(m), width=nr, space=1/nr, border=NA, ylim=c(mi.d, ma.d), col=col.data)
      barplot(m, beside=TRUE, add=TRUE, border=FALSE, col=col.pc, legend.text=lgt1,
              axes=FALSE)
    } 
    if (type %in% c(2, 3)) { # pc - data
      ma.p <- max(c(m, rowSums(m)))
      mi.p <- min(c(m, rowSums(m)))
      barplot(rowSums(m), width=nc, space=1/nc, border=NA, ylim=c(mi.p, ma.p), col=col.pc)
      barplot(t(m), beside=TRUE, add=TRUE, border=FALSE, col=col.data, legend.text=lgt2,
              axes=FALSE)
    } 
  } 
  return(invisible(m))
}
# ==============================================================================
# ==                                                                          ==
# ==                    the box.gs plot                                       ==
# ==                                                                          ==
# ==============================================================================


#' boxplot of gene set variables across all samples.
#' 
#' boxplot to show the variables (e.g. gene expression) of a gene set across
#' all samples.
#' 
#' This is a convenient function used to explore the expression of a set of
#' features/genes
#' 
#' @param x An object of calss \code{\link{mgsa-class}} or
#' \code{\link{moa.sup-class}}
#' @param gs Gene set want to be explored
#' @param moa An obejct of class \code{\link{moa}}. It is required if x is an
#' object of class \code{\link{moa.sup-class}}
#' @param col The coler code for samples
#' @param layout The layout control, see examples.
#' @param plot A logical indicates whether the result should be ploted. If
#' FALSE, a list of expression matrix of the gene set genes is returned.
#' Otherwise nothing returned.
#' @param obs.order Can be used to reorder the martrix, could be used when
#' clustering result is available.
#' @param \dots The arguments passed to \code{\link{boxplot}}
#' @return Do not return anything (plot=TRUE) or return a list of matrix
#' (plot=FALSE) depends on plot arugment.
#' @author Chen meng
#' @export
#' @examples
#' 
#'   # library(mogsa)
#'   # loading gene expression data and supplementary data
#'   data(NCI60_4array_supdata)
#'   data(NCI60_4arrays)
#'   mgsa <- mogsa(x = NCI60_4arrays, sup=NCI60_4array_supdata, nf=9,
#'                 proc.row = "center_ssq1", w.data = "inertia", statis = TRUE)
#' 
#'   allgs <- colnames(NCI60_4array_supdata[[1]])
#'   colcode <- as.factor(sapply(strsplit(colnames(NCI60_4arrays$agilent), split="\\."), "[", 1))
#'   a <- box.gs.feature(x=mgsa, gs=allgs[5], type=3, col=colcode, plot=FALSE)
#'   box.gs.feature(x=mgsa, gs=allgs[5], type=3, col=colcode, plot=TRUE, layout=matrix(1:4, 2, 2))
#' 
box.gs.feature <- function(x, gs, moa=NULL, col=1, layout=NULL, plot=TRUE, obs.order=NULL, ...) {
  # x - either mgsa or moa.sup
  if (inherits(x, "mgsa")) {
    if (!is.null(moa))
      warning("x is an object of class mgsa, moa argument is ignored!")
    data <- x@moa@data
    x <- x@sup
  } else if (inherits(x, "moa.sup")) {
    if (!inherits(moa, "moa"))
      stop("x is an object of class moa.sup, moa argument is required!")
    data <- moa@data
  }

  if (is.null(obs.order))
    obs.order <- 1:ncol(data[[1]])
  col <- rep(col, length(obs.order))[obs.order]

  gsidx <- lapply(x@sup, function(x) as.logical(x[, gs]))
  mats <- mapply(SIMPLIFY=FALSE, function(x, i) x[i, obs.order], x=data, i=gsidx)
# sapply(mats, dim)
  n <- length(mats)
  if (is.null(layout))
    layout <- matrix(1:n, 1, n)
  if (plot) {
    layout(layout)
    for (i in 1:n) {
      if(length(mats[[i]]) == 0)
        plot(0, pch=NA) else
          boxplot(mats[[i]], col=col , ...)
    }
  } else {
    return (mats)
  }
}
