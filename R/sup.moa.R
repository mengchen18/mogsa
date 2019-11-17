#' Projecting supplementary tables on object of class \code{moa-class}.
#' 
#' Projecting supplementary tables on \code{\link{moa-class}}
#' 
#' Projecting supplementary tables on \code{\link{moa-class}}, for details see
#' reference.
#' 
#' @param X An object of class \code{\link{moa-class}}
#' @param sup A list of data.frames contains supplementary data.
#' @param nf The number of principal components used in the projection.
#' @param factors The index of principal components used in the projection,
#' used when non-consecutive PC to be included in the analysis.
#' @param ks.stat The logical indicates if the p-value should be calculated
#' using K-S statistic (the method used in "ssgsea" in GSVA package).  Default
#' is FALSE, which means using the z-score method.
#' @param ks.B An integer to indicate the number of bootstrapping samples to
#' calculated the p-value of KS statistic.
#' @param ks.cores An integer indicate the number of cores to be used in
#' bootstrapping. It is passed to function \code{mclapply} in the
#' \code{parallel} package.
#' @param p.adjust.method The method of p value adjustment, passed to
#' \code{p.adjust} function.
#' @return An object of class \code{\link{moa.sup-class}}.
#' @author Chen Meng
#' @references Herve Abdi, Lynne J. Williams, Domininique Valentin and Mohammed
#' Bennani-Dosse. STATIS and DISTATIS: optimum multitable principal component
#' analysis and three way metric multidimensional scaling. WIREs Comput Stat
#' 2012. Volume 4, Issue 2, pages 124-167 Haenzelmann, S., Castelo, R. and
#' Guinney, J. GSVA: Gene set variation analysis for microarray and RNA-Seq
#' data. BMC Bioinformatics, 14:7, 2013.  Barbie, D.A. et al. Systematic RNA
#' interference reveals that oncogenic KRAS-driven cancers require TBK1.
#' Nature, 462(5):108-112, 2009.
#' @keywords data projection supplementary data
#' @export
#' @examples
#' 
#'     # library(mogsa)
#'     # loading gene expression data and supplementary data
#'     data(NCI60_4array_supdata)
#'     data(NCI60_4arrays)
#'     # check the dimension of each supplementary data to see how many gene set annotated the data
#'     sapply(NCI60_4array_supdata, dim)
#'     # run analysis
#'     ana <- moa(NCI60_4arrays, proc.row = "center_ssq1", w.data = "inertia", statis = TRUE)
#'     plot(ana, value="eig")
#'     # projectin supplementary data
#'     smoa <- sup.moa(ana, sup=NCI60_4array_supdata, nf=3)
#'     # heatmap visualize the gene set scores
#'     heatmap(slot(smoa, "score"))
#' 
sup.moa <- function(X, sup, nf = 2, factors = NULL, 
  ks.stat=FALSE, ks.B = 1000, ks.cores = NULL, p.adjust.method = "fdr") {

  if (is.null(nf) & is.null(factors))
    stop("nf or factors need to be specified.")

  if (!is.null(factors)) {
    if (!is.null(nf))
      message("factors is given, nf will be ignored.")
    pcomp <- factors
  } else if (!is.null(nf)) 
    pcomp <- 1:nf

  # sup, nf
  N <- length(sup)
  fs <- X@fac.scr 
  
  # sup data and moa data have the same order
  repr <- sapply(X@data, dim)[1, ]
  nn <- names(X@data)
  load <- split(X@loading, f=rep(nn, repr))
  load <- load[nn]
  
  cols <- sapply(sup, colSums)
  if (!is.matrix(cols))
    cols <- matrix(cols, nrow = length(sup))
  w <- rowSums(cols)

  if (any(w == 0)) {
    message("unrelated gene sets are detected and will be removed")
    sup <- lapply(sup, function(x) x[, w > 0])
  }
    
  # normsup <- lapply(sup, function(x, w) {
  #   sweep(x, 2, w, "/")
  # }, w=w)
  normsup <- sup # change here

  GSCoordinate_sep <- mapply(SIMPLIFY=FALSE, function(load, sup, A) {
    a <- t(sup * A) %*% as.matrix(load[, pcomp, drop=FALSE])
    colnames(a) <- paste("PC", pcomp, sep="")
    return(a)
  }, load=load, sup=normsup, A = split(X@w.data, names(X@w.data))[nn])

  GSCoordinate_comb <- Reduce("+", GSCoordinate_sep)
  
  contribution <- lapply(GSCoordinate_sep, function(supcor, score) {
    a <- lapply(1:length(pcomp), function(i) {
      r <- outer(supcor[, i], score[, pcomp[i]])
      colnames(r) <- rownames(score)
      return(r)
    })
    a[is.na(a)] <- 0
    names(a) <- paste("PC", pcomp, sep="")
    return(a)
  }, score=fs)
  
  contribution_dataset <- lapply(contribution, function(x) {
    Reduce("+", x)
  })
  
  contribution_pc <- lapply(1:length(pcomp), function(i, cont) {
    a <- lapply(cont, function(x) x[[i]])
    Reduce("+", a)
  }, cont=contribution)
  names(contribution_pc) <- paste("PC", pcomp, sep="")
  
  contribution_total <- Reduce("+", contribution_dataset) 

  csup <- do.call("rbind", sup)
  if (!ks.stat) {
    pmat <- .signifGS(X=X, sup=csup, A = X@w.data, 
      score=contribution_total, factors=pcomp)
    attr(pmat, "method") <- "zscore"
    } else {
      if (is.null(ks.cores)) 
        ks.cores <- getOption("mc.cores", 2L) 
      cat("running bootstrapping for p values of KS.stat ...\n")
      pmat <- .ks.pval(X, sup, ks.B=ks.B, A = X@w.data, factors=pcomp, mc.cores = ks.cores)
      attr(pmat, "method") <- "KS.stat"
    }
  
  pmatadj <- matrix(p.adjust(pmat, method = p.adjust.method), nrow(pmat), ncol(pmat))
  attr(pmatadj, "method") <- p.adjust.method

  res <- new("moa.sup", 
    sup = sup,
    coord.comb = GSCoordinate_comb,
    coord.sep = GSCoordinate_sep,
    score = contribution_total,
    score.data = contribution_dataset,
    score.pc = contribution_pc,
    score.sep = contribution,
    p.val = pmat,
    p.val.corrected = pmatadj
    )
  return(res)
}


.signifGS <- function(X, sup, A, score, factors) {
  
  # define function 
  ff <- function(x, n, score, infinite=FALSE) {
    lx <- length(x) 
    if (infinite)
      sf <- 1 else
        sf <- sqrt((lx-n)/(lx-1))
    sum_sd <- sf * sd(x)/sqrt(n) * n
    sum_mean <- mean(x) * n
    pp <- abs(pnorm(score, mean = sum_mean, sd = sum_sd))
    2 * rowMin(cbind(pp, 1-pp))
  }
  # reconstuct matrix using nf PCs
  U <- as.matrix(X@loading[, factors, drop=FALSE]) 
  D <- diag(sqrt(X@eig[factors]), nrow = length(factors))
  V <- as.matrix(X@eig.vec[, factors, drop=FALSE])
  rec <- (U %*% D %*% t(V)) * A
  # the number of feature in each GS
  supn <- colSums(sup != 0)
  # calculate the P value
  pmat <- sapply(1:ncol(rec), function(i) ff(rec[, i], supn, score = score[, i]))
  if (!is.matrix(pmat))
    pmat <- matrix(pmat, ncol = ncol(rec))
  
  colnames(pmat) <- colnames(score)
  rownames(pmat) <- rownames(score)
  return(pmat)
}
