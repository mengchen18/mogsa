utils::globalVariables(c("convertIdentifiers", "nodes"))

#' Combine two objects of class \code{mgsa} into one.
#' 
#' This function could only be used to combine two "mgsa" objects at present;
#' using "Reduce" function to combine more.
#' 
#' 
#' @name combine-methods
#' @aliases combine combine-methods combine,mgsa,mgsa-method
#' @docType methods
#' @param x one mgsa object
#' @param y another mgsa object
#' @param ...  ignored. Only two mgsa objects could be combined, using "Reduce"
#' to combine more than two sets.
#' @return A combined object of class \code{mgsa} will be returned.
#' @section Methods: \describe{ \item{list("signature(x = \"mgsa\", y =
#' \"mgsa\")")}{ To combine two objects of \code{mgsa}.  } This function could
#' only be used to combine two "mgsa" objects; using "Reduce" function to
#' combine more.  }
#' @keywords combine mogsa mgsa-class
#' @examples
#' 
#'   # library(mogsa)
#'   # loading gene expression data and supplementary data
#'   data(NCI60_4array_supdata)
#'   data(NCI60_4arrays)
#'   # split gene set annotation into two sets.
#'   sup1 <- lapply(NCI60_4array_supdata, function(x) x[, 1:10])
#'   sup2 <- lapply(NCI60_4array_supdata, function(x) x[, -(1:10)])
#'   # project two sets of annotation
#'   mgsa1 <- mogsa(x = NCI60_4arrays, sup=sup1, nf=9,
#'                 proc.row = "center_ssq1", w.data = "inertia", statis = TRUE)
#'   mgsa2 <- mogsa(x = NCI60_4arrays, sup=sup2, nf=9,
#'                  proc.row = "center_ssq1", w.data = "inertia", statis = TRUE)
#'   # combine two indenpendent mgsa sets
#'   mgsa_comb <- combine(mgsa1, mgsa2)
#'   dim(getmgsa(mgsa1, "score"))
#'   dim(getmgsa(mgsa2, "score"))
#'   dim(getmgsa(mgsa_comb, "score"))
#' 
NULL





#' Class \code{"mgsa"}
#' 
#' mgsa class here.
#' 
#' 
#' @name mgsa-class
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new("mgsa", ...)}.
#' @author Chen Meng
#' @seealso \code{\linkS4class{moa}} and \code{\linkS4class{moa.sup}}
#' @keywords classes
#' @examples
#' 
#'   showClass("mgsa")
#'   # library(mogsa)
#'   # loading gene expression data and supplementary data
#'   data(NCI60_4array_supdata)
#'   data(NCI60_4arrays)
#'   # split gene set annotation into two sets.
#'   sup1 <- lapply(NCI60_4array_supdata, function(x) x[, 1:10])
#'   sup2 <- lapply(NCI60_4array_supdata, function(x) x[, -(1:10)])
#'   # project two sets of annotation
#'   mgsa1 <- mogsa(x = NCI60_4arrays, sup=sup1, nf=9,
#'                 proc.row = "center_ssq1", w.data = "inertia", statis = TRUE)
#'   mgsa2 <- mogsa(x = NCI60_4arrays, sup=sup2, nf=9,
#'                  proc.row = "center_ssq1", w.data = "inertia", statis = TRUE)
#'   # combine two indenpendent mgsa sets
#'   mgsa_comb <- combine(mgsa1, mgsa2)
#'   dim(getmgsa(mgsa1, "fac.scr"))
#'   dim(getmgsa(mgsa2, "fac.scr"))
#'   dim(getmgsa(mgsa_comb, "fac.scr"))
#' 
NULL





#' Class \code{"moa"}
#' 
#' moa class object
#' 
#' 
#' @name moa-class
#' @aliases moa-class plot,moa,missing-method
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new("moa", ...)}. %% ~~ describe objects here ~~
#' @author Chen Meng
#' @references Herve Abdi, Lynne J. Williams, Domininique Valentin and Mohammed
#' Bennani-Dosse. STATIS and DISTATIS: optimum multitable principal component
#' analysis and three way metric multidimensional scaling. WIREs Comput Stat
#' 2012. Volume 4, Issue 2, pages 124-167
#' 
#' Herve Abdi, Lynne J. Williams, Domininique Valentin. Multiple factor
#' analysis: principal component analysis for multitable and multiblock data
#' sets. WIREs Comput Stat 2013
#' @keywords classes
#' @examples
#' 
#'     showClass("moa")
#'     # load("R/mogsa/data/NCI60_4arrays.rda")
#'     data(NCI60_4arrays)
#'     ana <- moa(NCI60_4arrays, proc.row = "center_ssq1", w.data = "inertia", statis = TRUE)
#' 
#'     plot(ana, value="eig")
#'     plot(ana, value="tau", type=2)
#' 
NULL





#' Class \code{"moa.sup"}
#' 
#' moa.sup class desc.
#' 
#' 
#' @name moa.sup-class
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new("moa.sup", ...)}.
#' @author Chen Meng
#' @seealso objects to See Also as \code{\link{decompose.gs.ind}},
#' \code{\link{box.gs.feature}}, \code{\link{plotGS}},
#' \code{\link{decompose.gs.group}}.
#' @keywords classes
#' @examples
#' 
#' 	showClass("moa.sup")
#' 	data(NCI60_4array_supdata)
#' 	data(NCI60_4arrays)
#' 
#' 	sapply(NCI60_4array_supdata, dim)
#' 	ana <- moa(NCI60_4arrays, proc.row = "center_ssq1", w.data = "inertia", statis = TRUE)
#' 	plot(ana, value="eig")
#' 	smoa <- sup.moa(ana, sup=NCI60_4array_supdata, nf=5)
#' 
NULL





#' Multiple omics clustering and gene set analysis
#' 
#' Modern "omics" technologies enable quantitative monitoring of the abundance
#' of various biological molecules in a high-throughput manner, accumulating an
#' unprecedented amount of quantitative information on a genomic scale. Gene
#' set analysis is a particularly useful method in high throughput data
#' analysis since it can summarize single gene level information into the
#' biological informative gene set levels. This package provide a method do the
#' gene set analysis based on multiple omics data that describing the same set
#' of observations/samples.
#' 
#' \tabular{ll}{ Package: \tab mogsa\cr Type: \tab Package\cr Version: \tab
#' 1.3.1\cr Date: \tab 2016-01-19\cr License: \tab GPL-2\cr Depends: \tab
#' methods\cr } The main function in the package is "mogsa", see the function
#' help manu for more details.
#' 
#' @name mogsa-package
#' @docType package
#' @author Chen Meng Maintainer: Chen Meng <chen.meng@@tum.de>
#' @references Chen Meng, Dominic Helm, Martin Frejno, and Bernhard Kuster.
#' moCluster: Identifying Joint Patterns Across Multiple Omics Data Sets.
#' Journal of Proteome Research 2016.
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
#'   # using moa as input
#'   ana <- moa(NCI60_4arrays, proc.row = "center_ssq1", w.data = "inertia", statis = TRUE)
#'   smoa <- sup.moa(ana, sup=NCI60_4array_supdata, nf=3)
#'   mgsa2 <- mogsa(x = ana, sup=NCI60_4array_supdata, nf=9)
#'   mgsa3 <- mogsa(x = ana, sup=smoa)
#' 
NULL





#' supp data for Microarray gene expression profiles of the NCI 60 cell lines
#' from 4 different platforms
#' 
#' Supplmentary to NCI60_4arrays.
#' 
#' 
#' @name NCI60_4array_supdata
#' @docType data
#' @format The format is: List of 4 \code{matrix} \itemize{
#' \item\$agilent:\code{matrix} containing 300 rows and 60 columns.  300 gene
#' expression log ratio measurements of the NCI60 cell lines, by Agilent
#' platform.
#' 
#' \item\$hgu133:\code{matrix} containing 298 rows and 60 columns.  298 gene
#' expression log ratio measurements of the NCI60 cell lines, by H-GU133
#' platform.
#' 
#' \item\$hgu133p2:\code{matrix} containing 268 rows and 60 columns.  268 gene
#' expression log ratio measurements of the NCI60 cell lines, by H-GU133 plus
#' 2.0 platform.
#' 
#' \item\$hgu95:\code{matrix} containing 288 rows and 60 columns.  288 gene
#' expression log ratio measurements of the NCI60 cell lines, by H-GU95
#' platform.  }
#' @return \code{NCI60_4array_supdata} will be loaded in your working space.
#' @keywords datasets NCI-60 Microarray
NULL





#' Microarray gene expression profiles of the NCI 60 cell lines from 4
#' different platforms
#' 
#' The 60 human tumour cell lines are derived from patients with leukaemia,
#' melanoma, lung, colon, central nervous system, ovarian, renal, breast and
#' prostate cancers. The cell line panel is widely used in anti-cancer drug
#' screen. In this dataset, a subset of microarray gene expression of the NCI
#' 60 cell lines from four different platforms are combined in a list, which
#' could be used as input to \code{mcia} directly.
#' 
#' 
#' @name NCI60_4arrays
#' @docType data
#' @format The format is: List of 4 \code{data.frame}s \itemize{
#' \item\$agilent:\code{data.frame} containing 300 rows and 60 columns.  300
#' gene expression log ratio measurements of the NCI60 cell lines, by Agilent
#' platform.
#' 
#' \item\$hgu133:\code{data.frame} containing 298 rows and 60 columns.  298
#' gene expression log ratio measurements of the NCI60 cell lines, by H-GU133
#' platform.
#' 
#' \item\$hgu133p2:\code{data.frame} containing 268 rows and 60 columns.  268
#' gene expression log ratio measurements of the NCI60 cell lines, by H-GU133
#' plus 2.0 platform.
#' 
#' \item\$hgu95:\code{data.frame} containing 288 rows and 60 columns.  288 gene
#' expression log ratio measurements of the NCI60 cell lines, by H-GU95
#' platform.  }
#' @return \code{NCI60_4arrays} will be loaded in your working space.
#' @references Reinhold WC, Sunshine M, Liu H, Varma S, Kohn KW, Morris J,
#' Doroshow J, Pommier Y CellMiner: A Web-Based Suite of Genomic and
#' Pharmacologic Tools to Explore Transcript and Drug Patterns in the NCI-60
#' Cell Line Set.  Cancer Research. 2012 Jul, 15;72(14):3499-511
#' @source Cell Miner \url{http://discover.nci.nih.gov/cellminer/}
#' @keywords datasets NCI-60 Microarray
NULL





#' Methods for function \code{plot}
#' 
#' Methods for function \code{plot}
#' 
#' 
#' @name plot-methods
#' @docType methods
#' @section Methods: \describe{
#' 
#' \item{list("signature(x = \"moa\", y = \"missing\")")}{ plot "moa" object }
#' Argument "value" sould be one of "eig", "tau", "obs", "var" and "RV"\
#' 
#' if value = "eig", the eigenvalue would be plotted as scree plot. The
#' following arguments could be set:\
#' 
#' type=1 - The type of plot to show eigenvalues.  (type=1: the eigenvalue are
#' plotted; type=2: partial eigenvalue shown as concatenated bars; type=3:
#' partial eigenvalue shown as bars side by side; type=4: matplot view of
#' eigenvales, lty need to be set; type=5; the two dimensional plot of partial
#' eigenvalues, axes and pch need to be set in this case.) \ axes=NULL - The
#' axes selected to plot \ n=NULL - Top n eigenvalues to be drawn \ tol=1e-5 -
#' The tolerance of eigenvalue, eigenvalues lower than this value will not be
#' shown. \ legend=NULL - legend to put, a character string as calling legend
#' function \ col=NULL - The color of partial eigenvalues from each data set \
#' lty=1 - The line type used in the matplot, used when type =4 \ pch=NULL -
#' the pch to draw 2D partial eigen plot, when type = 5 used \ lg.x="topright"
#' - The position of legend \ lg.y=NULL - Poistion argument passed to function
#' "legend" \ ... - other arguemnts passed to functions \ \
#' 
#' if value = "tau", the same with eig, but in the eigenvalues are scaled to 1
#' \
#' 
#' if value = "obs", the observation space will be shown, the following
#' argument could be set:\
#' 
#' axes=1:2 - Which axes should be draw\ type=1 - Which type, see below (for
#' type=1: the center points draw; type=2: the separate factor scores linked by
#' lines; ... will be passed to function "points")\ data.pch=20 - the pch of
#' dataset, if type=1, the first one is used\ col=1 - the color of
#' observations, recycled used by data.frame\ label=FALSE - A logical indicates
#' if labels should be shown\ lg.x="topright" - Position of legend \ lg.y=NULL
#' - Position of legend \ xlim=NULL - The x limit \ ylim=NULL - The y limit \
#' label.cex=1 - the cex of text \ ... \
#' 
#' var - the separate gene view, layout can be specified \
#' 
#' RV - the heatmap of RV coefficients }
#' @keywords methods moa-class
NULL





#' Methods for function \code{print}
#' 
#' Methods for function \code{print}
#' 
#' 
#' @name print-methods
#' @aliases print-methods print print,moa-method print,moa.sup-method
#' print,mgsa-method
#' @docType methods
#' @section Methods: \describe{
#' 
#' \item{list("signature(object = \"moa\")")}{ print "moa" class }
#' 
#' \item{list("signature(object = \"moa.sup\")")}{ print "sup.moa" class }
#' 
#' \item{list("signature(object = \"mgsa\")")}{ print "mgsa" class } }
#' @keywords methods moa-class sup.moa-class mgsa-class
NULL





#' Methods for function \code{show}
#' 
#' Methods for function \code{show}
#' 
#' 
#' @name show-methods
#' @aliases show-methods show show,moa-method show,moa.sup-method
#' show,mgsa-method
#' @docType methods
#' @section Methods: \describe{
#' 
#' \item{list("signature(object = \"moa\")")}{ show "moa" class }
#' 
#' \item{list("signature(object = \"moa.sup\")")}{ show "sup.moa" class }
#' 
#' \item{list("signature(object = \"mgsa\")")}{ show "mgsa" class } }
#' @keywords methods moa-class sup.moa-class mgsa-class
NULL





#' Methods for function \code{summary}
#' 
#' Methods for function \code{summary}
#' 
#' 
#' @name summary-methods
#' @aliases summary-methods summary summary,moa-method summary,moa.sup-method
#' summary,mgsa-method
#' @docType methods
#' @section Methods: \describe{
#' 
#' \item{list("signature(object = \"moa\")")}{ summary "moa" class }
#' 
#' \item{list("signature(object = \"moa.sup\")")}{ summary "sup.moa" class }
#' 
#' \item{list("signature(object = \"mgsa\")")}{ summary "mgsa" class } }
#' @keywords methods moa-class sup.moa-class mgsa-class
NULL



