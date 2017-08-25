# both rows and columns could be sparse

mbpca2 <- function (
  
  x, ncomp, method, 
  kp = "all", kt = "all", 
  unit.p = FALSE, unit.t = FALSE, 
  center = TRUE, scale = FALSE, 
  option = "uniform", maxiter = 1000, moa = TRUE, verbose = TRUE, 
  svd.solver = c("svd", "fast.svd", "propack")) {
  
    x <- lapply(x, t)
    call <- match.call()
    x <- mogsa:::processOpt(x, center = center, scale = scale, option = option)
    nc <- sapply(x, ncol)
    keepAllb <- kp[1] == "all"
    keepAllt <- kt[1] == "all"
    
    prddata <- lapply(x, t)
    ssl <- match.arg(svd.solver)
    svdf <- switch(ssl, svd = svd, fast.svd = fast.svd, 
                   propack = function(X) propack.svd(X, neig = 1, opts = list(kmax = 20)))
    for (i in 1:ncomp) {
      if (verbose) 
        cat(paste("calculating component ", i, " ...\n", 
                  sep = ""))
      if (keepAllb & keepAllt) 
        r <- mogsa:::msvd(x, svd.sol = svdf)
      else {
        if (keepAllb == "all")
          keepAllb <- Inf
        if (keepAllt == "all")
          keepAllt <- Inf
        if (length(kp) < length(x))
          kp <- rep(kp, length.out = length(x))
        if (length(kt) < length(x))
          kt <- rep(kt, length.out = length(x))
        r <- biSoftK(x, maxiter = maxiter, kp = kp, kt = kt, unit.pb = unit.p, unit.tb = unit.t)
      }
      x <- mogsa:::deflat(x, r$t, r$tb, r$pb, method)
      if (i == 1) 
        res <- r
      else {
        res$t <- cbind(res$t, r$t)
        res$w <- cbind(res$w, r$w)
        res$tb <- mapply(cbind, res$tb, r$tb, SIMPLIFY = FALSE)
        res$pb <- mapply(cbind, res$pb, r$pb, SIMPLIFY = FALSE)
      }
    }
    if (moa) 
      res <- mogsa:::toMoa(prddata, res, call = call)
    return(res)
  }
