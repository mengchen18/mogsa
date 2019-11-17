# moa - multiple omics data analysis 
# 
# input arguments
# data 
#   a list of data.frame or matrix that rows represent variables (genes)
#   and coloumns represent measurements/observations (samples, cell lines)
# method
#   either "statis" or "mfa"
# full.analysis
#   a logical indicate wether the full analysis should be performed
#   set as TRUE in default. In the bootstrapping preocedure, it set as
#   FALSE
# 
# detail
#   if full analysis is FALSE, returns
#     eig - eigen values
#     eig.vec - eigen vector (columns, in omics data, samples/observations)
#     loading - eigen vector for rows (in omics data, genes/molecules)
#     fac.scr - factor score
#     tab.dim - the dimension of each table
#     method - which method, "statis" or "mfa"
#     full.analysis - whether full analysis is performed
#   if full analysis is TRUE, returns
#     eig, eig.vec, loading, tab.dim, method, full.analysis, fac.scr -see above
#     tau - explained variance (percentage) by each eigenvector
#     fac.scr - factor.scores, calculated see abdi. 2012 statis ...
#     partial.fs - partial factor scores
#     ctr.obs - contribution of observations/samples
#     ctr.var - contribution of variables
#     ctr.tab - contribution of tables
#     partial.eig - partial eigen value



#' Multiple omics data analysis using MFA or STATIS
#' 
#' Analysis multiple omics data using MFA or STATIS. The input multiple tables
#' are in a form that columns are samples and rows are variables/features.
#' 
#' Different methods employs different precessing of row and datasets.  For
#' multipple factorial analysis (MFA), the rows of each dataset are first
#' centered and scaled, then each dataset is weighted by the reverse of its
#' first eigenvalue (proc.row=center_ssq1, w.data="lambda1"). This algorithm
#' does not have a well defined criterion to be optimized (see reference).
#' 
#' If statis=TRUE, the statis algorithm will be used, that is, each dataset
#' will be further weighted so that datasets closer to the overall structure
#' will receive a higher weight.
#' 
#' @param data A list of \code{data.frame} or \code{matrix} that contains the
#' input datas, the columns in all datasets should be samples/observations
#' (which need to be matched) and rows should be variables.
#' @param proc.row Preprocessing of rows of datasets, should be one of
#' \code{none} - no preprocessing, \code{center} - center only,
#' \code{center_ssq1} - center and scale (sum of squred values equals 1),
#' \code{center_ssqN} - center and scale (sum of squred values equals the
#' number of columns), \code{center_ssqNm1} - center and scale (sum of squred
#' values equals the number of columns - 1) MFA corresponds to
#' "proc.row=center_ssq1" and 'w.data="lambda1"'
#' @param w.data The weights of each separate dataset, should be one of
#' \code{uniform} - no weighting,
#' 
#' \code{lambda1} - weighted by the reverse of the first eigenvalue of each
#' individual dataset
#' 
#' or \code{inertia} - weighted by the reverse of the total inertia.  See
#' detail.
#' @param w.row If it is not null, it should be a list of positive numerical
#' vectors, the length of which should be the same with the number of rows of
#' each dataset to indicated the weight of rows of datasets.
#' @param statis A logical indicates whether STATIS method should be used. See
#' details.
#' @param moa Logical; whether the output should be converted to an object of
#' class \code{\link{moa-class}}
#' @return An object of class \code{\link{moa-class}}.
#' @author Chen Meng
#' @seealso \code{\link{sup.moa}}, \code{\link{mogsa}}. More about plot see
#' \code{\link{moa-class}}.
#' @references Herve Abdi, Lynne J. Williams, Domininique Valentin and Mohammed
#' Bennani-Dosse. STATIS and DISTATIS: optimum multitable principal component
#' analysis and three way metric multidimensional scaling. WIREs Comput Stat
#' 2012. Volume 4, Issue 2, pages 124-167 Herve Abdi, Lynne J. Williams,
#' Domininique Valentin. Multiple factor analysis: principal component analysis
#' for multitable and multiblock data sets. WIREs Comput Stat 2013
#' @keywords PCA MVA MFA STATIS
#' @export
#' @examples
#' 
#' 
#'   # library(mogsa)
#'   # loading data
#'   data(NCI60_4arrays)
#'   # run analysis
#'   ana <- moa(NCI60_4arrays, proc.row = "center_ssq1", w.data = "inertia", statis = TRUE)
#'   # plot
#'   # plot eigen value
#'   plot(ana, value = "eig", type = 2)
#'   # plot the normalized (percentage) eigen value
#'   plot(ana, value = "tau", type = 2)
#'   # ploting the observations
#'   colcode <- as.factor(sapply(strsplit(colnames(NCI60_4arrays$agilent), split="\\."), "[", 1))
#'   plot(ana, type = 1, value = "obs", col=colcode)
#'   plot(ana, type = 2, value = "obs", col=colcode, data.pch=1:4)
#'   # plot variables/features in each data sets
#'   plot(ana, value = "var", layout=matrix(1:4, 2, 2))
#'   # plot the RV coefficients for the data sets
#'   plot(ana, value = "RV")
#' 
#'   # to extract the components representing significant concordance structures between datasets
#'   bt <- bootMoa(moa = ana, proc.row = "center_ssq1", w.data = "inertia", statis = TRUE, B = 20)
#' 
moa <- function(data, proc.row="center_ssq1", w.data="inertia", w.row=NULL, statis=FALSE, moa=TRUE) {
  
  kd <- data  
  data <- lapply(data, as.matrix)
  nRows <- sapply(data, nrow)

  # check 
  if (!is.null(w.row)) {
    wr <- do.call("c", w.row)
    if (any(wr <= 0))
      stop ("postive row weights are required.")
    if (! all.equal(sapply(w.row, length), nRows))
      stop("w.row should be a list of vector that contain the weight for each row in data! The length of 
           the element of w.row should equals the # of rows in each data.")
  }
  
  checkRowName <- sapply(data, function(x) is.null(rownames(x)))
  checkColName <- sapply(data, function(x) is.null(colnames(x)))
  checkTabName <- is.null(names(data))
  if (any(c(checkRowName, checkColName, checkTabName)))
    stop ("Table name or Colnames or Rownames are missing in all/one data!")
  nCols <- sapply(data, ncol)
  if (length(unique(nCols)) != 1) 
    stop ("Unequal number of columns! Columns need to be matched.")
  
  # preprocessing of row 
  nObs <- unique(nCols)
  nData <- length(data)
  data <- lapply(data, .proc.row, method = proc.row)
  
  # create weight
  wObs <- rep(1/nObs, nObs)
  wD <- .w.data(data, method = w.data, statis=statis)
  wData <- rep(sqrt(wD$w), nRows)
  if (!is.null(w.row))
    wData <- wData * sqrt(wr)

  # data concatenation, weighting and decomposition
  d1 <- .concateTabs(wD$x)
  Xt <- d1$tab * wData / sqrt(nObs)
  sing <- svd(Xt)
  
  if (!moa)
    return(sing)

  decom <- .read.svd(x=sing, data=wD$x, 
                      M=wObs, A=wData^2,
                      design=rep(names(data), nRows))

  decom$data <- kd
  decom$tab.dim <- data.frame(sapply(data, dim), row.names=c("row", "col"))
  decom$call <- match.call()
  
  decom$proc.row <- proc.row
  decom$w.data <- w.data
  decom$w.row <- w.row
  
  colnames(decom$partial.eig) <- paste("PC", 1:ncol(decom$partial.eig), sep="")
  result <- new("moa",
                eig = decom$eig,
                tau = decom$tau,
                partial.eig = decom$partial.eig,
                eig.vec = decom$eig.vec,
                fac.scr = decom$fac.scr,
                loading = decom$loading,
                partial.fs = decom$partial.fs,
                ctr.obs = decom$ctr.obs,
                ctr.var = decom$ctr.var,
                ctr.tab = decom$ctr.tab,
                RV = decom$RV,
                w.row = decom$alpha,
                w.data = wData^2,
                data = decom$data,
                tab.dim = decom$tab.dim,
                call = decom$call)
  return(result)
}

