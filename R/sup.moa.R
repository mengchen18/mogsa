# Input arguments:
#   X - the output of mpSTATIS
#   sup - the supplementary table (nVar * nIndividual)
#   axes - which axes should be used

sup.moa <- function(X, sup, nf=2, 
  ks.stat=FALSE, ks.B = 1000, ks.cores = NULL) {

  # sup, nf
  N <- length(sup)
  fs <- X@fac.scr 
  
  # sup data and moa data have the same order
  repr <- sapply(X@data, dim)[1, ]
  nn <- names(X@data)
  load <- split(X@loading, f=rep(nn, repr))
  load <- load[nn]
  
  w <- rowSums(sapply(sup, colSums))
  if (any(w == 0)) 
    stop("unrelated pathways involved")
  # normsup <- lapply(sup, function(x, w) {
  #   sweep(x, 2, w, "/")
  # }, w=w)
  normsup <- sup # change here

  GSCoordinate_sep <- mapply(SIMPLIFY=FALSE, function(load, sup, A) {
    a <- t(sup * A) %*% as.matrix(load[, 1:nf, drop=FALSE])
    colnames(a) <- paste("PC", 1:nf, sep="")
    return(a)
  }, load=load, sup=normsup, A = split(X@w.data, names(X@w.data)))

  GSCoordinate_comb <- Reduce("+", GSCoordinate_sep)
  
  contribution <- lapply(GSCoordinate_sep, function(supcor, score) {
    a <- lapply(1:nf, function(i) {
      r <- outer(supcor[, i], score[, i])
      colnames(r) <- rownames(score)
      return(r)
    })
    a[is.na(a)] <- 0
    names(a) <- paste("PC", 1:nf, sep="")
    return(a)
  }, score=fs)
  
  contribution_dataset <- lapply(contribution, function(x) {
    Reduce("+", x)
  })
  
  contribution_pc <- lapply(1:nf, function(i, cont) {
    a <- lapply(cont, function(x) x[[i]])
    Reduce("+", a)
  }, cont=contribution)
  names(contribution_pc) <- paste("PC", 1:nf, sep="")
  
  contribution_total <- Reduce("+", contribution_dataset) 

  csup <- do.call("rbind", sup)
  if (!ks.stat) {
    pmat <- .signifGS(X=X, sup=csup, A = X@w.data, 
      score=contribution_total, nf=nf)
    attr(pmat, "method") <- "zscore"
    } else {
      if (is.null(ks.cores)) 
        ks.cores <- getOption("mc.cores", 2L) 
      cat("running bootstrapping for p values of KS.stat ...\n")
      pmat <- .ks.pval(X, sup, ks.B=ks.B, A = X@w.data, nf=nf, mc.cores = ks.cores)
      attr(pmat, "method") <- "KS.stat"
    }
  
  res <- new("moa.sup", 
    sup = sup,
    coord.comb = GSCoordinate_comb,
    coord.sep = GSCoordinate_sep,
    score = contribution_total,
    score.data = contribution_dataset,
    score.pc = contribution_pc,
    score.sep = contribution,
    p.val = pmat
    )
  return(res)
}


.signifGS <- function(X, sup, A, score, nf) {
  
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
  U <- as.matrix(X@loading[, 1:nf, drop=FALSE]) 
  D <- diag(sqrt(X@eig[1:nf]), nrow = nf)
  V <- as.matrix(X@eig.vec[, 1:nf, drop=FALSE])
  rec <- (U %*% D %*% t(V)) * A
  # the number of feature in each GS
  supn <- colSums(sup != 0)
  # calculate the P value
  pmat <- sapply(1:ncol(rec), function(i) ff(rec[, i], supn, score = score[, i]))
  colnames(pmat) <- colnames(score)
  rownames(pmat) <- rownames(score)
  return(pmat)
}