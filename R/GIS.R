#' calculate gene influential scores of genes in a gene set.
#' 
#' Calculate the gene influential score of individual feature to the overall
#' variance of GS score.  Using a leave-one-out procedure (See detail).
#' 
#' The evaluation of the importance of a single feature is calculated in the
#' supervised or unsupervised manner.
#' 
#' In the unsupervise manner, the value is calculated by:
#' 
#' log2(var(GS_-i)/var(GS))
#' 
#' where GS is the gene set score, and the GS_-i is a recalculate of gene set
#' score without i'th feature. var() is the variance.
#' 
#' In the supervised manner, the value is caluclated as the F-ratio over a
#' class vector:
#' 
#' log2(F(GS_-i)/F(GS))
#' 
#' Where F() is the calculation of F-ratio. The unsupervised GIS is encouraged
#' since it works better for most of the cases in practice.
#' 
#' @param x An object of class \code{\link{mgsa-class}}.
#' @param geneSet A charater string or number to indicated the gene sets under
#' conserderation.
#' @param nf The number of PCs used in the caluculation of gene set scores.
#' The default is NA, which means using all the PCs in the mogsa. This should
#' work for most of the cases.
#' @param barcol The color of the bars, which is used to distinguish
#' features/genes from different datasets, so its length should be the same as
#' the number of data sets.
#' @param topN An positive integer specify the number of top influencers that
#' should to returned.
#' @param plot A logical indicate if the result should be plotted.
#' @param Fvalue A logical indicate if the GIS should be calculated in a
#' supervised manner.
#' @param ff The vector indicates the group of columns for calculating the
#' F-ratio when Fvalue=TRUE.
#' @param cor A logical indicates whether use correlation between reconstructed
#' expression with GSS.  This is faster than the standard GIS.
#' @return An object of class \code{data.frame} contains three columns. The
#' first column is the feature name, the second columns is the gene influential
#' score. The third columns indicates from where the feature/gene is selected.
#' @author Chen Meng
#' @seealso see \code{\link{annotate.gs}}
#' @references TBA
#' @export
#' @examples
#' 
#' 	# library(mogsa)
#' 	# loading gene expression data and supplementary data
#' 	data(NCI60_4array_supdata)
#' 	data(NCI60_4arrays)
#' 	mgsa <- mogsa(x = NCI60_4arrays, sup=NCI60_4array_supdata, nf=9,
#' 	              proc.row = "center_ssq1", w.data = "inertia", statis = TRUE)
#' 	allgs <- colnames(NCI60_4array_supdata[[1]])
#' 
#' 	# unsupervised measurement
#' 	GIS(mgsa, allgs[1], topN = 5)
#' 
#' 	# supervised measurement
#' 	tissueType <- as.factor(sapply(strsplit(colnames(NCI60_4arrays$agilent), split="\\."), "[", 1))
#' 	GIS(mgsa, allgs[1], topN = 5, Fvalue = TRUE, ff = tissueType)
#' 	# more PCs to calcualte
#' 	GIS(mgsa, allgs[1], nf = 20, topN = 5, Fvalue = TRUE, ff = tissueType)
#' 
GIS <- function(x, geneSet, nf=NA, barcol=NA, topN=NA, plot=TRUE, Fvalue=FALSE, ff=NA, cor=FALSE) {
  
  if (!inherits(x, "mgsa"))
    stop("x should be an object of class mgsa!")
  
  fvalue <- function(a, f) {
    n <- length(f)
    ng <- length(unique(f))
    ma <- mean(a)
    ssb <- sum(tapply(a, f, function(x) length(x)*(mean(x)-ma)^2))
    ssw <- sum(tapply(a, f, function(x) sum((x-mean(x))^2)))
    ssb*(n-ng)/ssw/(ng-1)
  }
  
  exprs <- x@moa@data
  rn <- lapply(exprs, rownames)
  rn <- unlist(rn)
  
  if (is.na(nf))
    nf <- length(x@sup@score.pc)
  if (is.na(barcol[1]))
    barcol <- 1:length(exprs)
  
  scor <- x@sup@score[geneSet, , drop=FALSE]
  gsi <- lapply(x@sup@sup, function(x) x[, geneSet])
  coln <- sapply(gsi, function(x) sum(x!=0))
  col_code <- rep(barcol, coln)
  gsi <- unlist(gsi)
  genes_idx <- which(gsi != 0)
  
  if (cor) {
    if (Fvalue)
      warning("cor = TRUE, Fvalue is ignored.")
    xmat <- as.matrix(x@moa@fac.scr[, 1:nf]) %*% t(as.matrix(x@moa@loading)[, 1:nf])
    ef <- cor(t(scor), xmat[, genes_idx], )
    cn <- colnames(ef)
    ef <- c(ef)
    names(ef) <- cn
  } else {
    if (!Fvalue) {
      sdscor <- sd(scor)
    } else if (Fvalue & length(ff)==length(scor)) {
      sdscor <- fvalue(scor, as.factor(ff))
    } else 
      stop("if Fvalue is TRUE, ff need to be defined!")
    gsimat <- sapply(genes_idx, function(i, gsi) {
      gsi[i] <- 0
      return(gsi)}, gsi=gsi)
    coor <- t(gsimat) %*% as.matrix(x@moa@loading)[, 1:nf]
    scor_rm <- as.matrix(x@moa@fac.scr[, 1:nf]) %*% t(coor)
    
    if (!Fvalue) {
      sdscor_rm <- apply(scor_rm, 2, sd)
    } else if (Fvalue & length(ff)==length(scor)) {
      sdscor_rm <- rowFtests(t(scor_rm), fac=as.factor(ff))[, "statistic"]
    }
    names(sdscor_rm) <- rn[genes_idx]
    ef <- -log2(sdscor_rm/sdscor)
  }  
  
  ef <- ef/max(ef)
  oef <- order(ef, decreasing=FALSE)
  ef <- ef[oef]
  ef <- rev(ef)
  
  if (plot) {
    barplot(rev(ef), horiz=TRUE, col=col_code[oef], border=col_code[oef], names.arg = NA)
    legend(x="bottomright", col=barcol, legend=names(coln), pch=15)
  }
  
  dn <- names(x@moa@data)
  fd <- as.factor(rev(col_code[oef]))
  if (!is.null(dn)) levels(fd) <- dn
  r <- data.frame(feature = names(ef), GIS = ef, data=fd)
  r$feature <- as.character(r$feature)
  
  if (!is.na(topN))
    r <- r[1:min(topN, nrow(r)), ]
  return(r)
}




#' Summary annotation information of a gene set
#' 
#' Retrive variables/features (genes) mapped to the annotated data sets in a
#' gene set. Also returns the the information about presence and absence of a
#' feature for a specific data set.
#' 
#' 
#' @param mgsa An object of class \code{\link{mgsa-class}} or
#' \code{\link{moa.sup-class}}.
#' @param gs The name of a geneset
#' @return This function returns a data.frame.  The first column shows the name
#' of features.  The last column is for the count of how many data sets has the
#' corresponding features. Columns in the middle contains logical value
#' indicating whether a feature is presented in a particular data set.
#' @author Chen Meng
#' @seealso see \code{\link{GIS}}
#' @examples
#' 
#'   # library(mogsa)
#'   # loading gene expression data and supplementary data
#'   data(NCI60_4array_supdata)
#'   data(NCI60_4arrays)
#'   mgsa <- mogsa(x = NCI60_4arrays, sup=NCI60_4array_supdata, nf=9,
#'                 proc.row = "center_ssq1", w.data = "inertia", statis = TRUE)
#'   allgs <- colnames(NCI60_4array_supdata[[1]])
#'   annotate.gs(mgsa, allgs[1])
#' 
annotate.gs <- function(mgsa, gs) {
  if (length(gs) > 1)
    stop("gs has to be an one element character string or integer")
  rn <- sapply(mgsa@moa@data, rownames)
  r <- lapply(mgsa@sup@sup, function(x) x[, gs] != 0)
  gs <- mapply(SIMPLIFY=FALSE, function(x, g) x[g], x=rn, g=r)
  var <- unique(unlist(gs))
  data <- sapply(gs, function(x) var %in% x) 
  stat <- rowSums(data)
  df <- data.frame(var=var, data=data, stat=stat)
  df <- df[order(df$stat, decreasing=TRUE), ]
  return(df)
}

