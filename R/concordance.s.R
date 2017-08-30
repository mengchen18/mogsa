sconcord <- function(x, y, ncomp = 1, kp = "all", kt = "all", dmod = 1, unit.p = TRUE, unit.t = TRUE, 
                     center = TRUE, scale = FALSE, option = "uniform") {
  
  ndata <- length(x)
  pb <- sapply(x, nrow)
  pbt <- sum(pb)
  nsample <- ncol(y)
  py <- nrow(y)
  
  f.pb <- rep(names(x), pb)
  f.tb <- rep(names(x), each = nsample)
  
  loading.x <- matrix(NA, pbt, ncomp)
  loading.y <- matrix(NA, py, ncomp)
  score.x <- matrix(NA, nsample*ndata, ncomp)
  score.y <- matrix(NA, nsample, ncomp)
  eig <- c()
  
  # preprocessing of x and y
  y <- t(scale(t(y), center = center, scale = scale))
  x <- mogsa:::processOpt(x, center = center, scale = scale, option = option)
  
  for (i in 1:ncomp) {
    cat(paste("calculating component", i, "...\n"))
    
    ## construct covariance/correlation matrix
    # ever loop if deflat x and/or y
    cmats <- lapply(x, tcrossprod, y)
    
    ## MBPCA
    res <- mbpca2(cmats, ncomp = 1, method = "globalScore", kp = kp, kt = kt, moa = FALSE, verbose = FALSE, 
                  unit.p = unit.p, unit.t = unit.t, center = FALSE, scale = FALSE)
    
    ## output
    t1 <- normvec(res$t)
    fac.y <- crossprod(y, t1)
    fac.x <- mapply(SIMPLIFY = FALSE, crossprod, x, res$pb)
    loading.x[, i] <- unlist(res$pb)
    score.x[, i] <- unlist(fac.x)
    loading.y[, i] <- t1
    score.y[, i] <- fac.y
    eig <- c(eig, res$eig)
    
    ## deflation
    delta.y <- tcrossprod(t1, fac.y)
    y <- y - delta.y
    
    # if (dmod == 1 | dmod == 2) {
    #   # only deflat Y, dmod = 1
    #   delta.y <- tcrossprod(t1, fac.y)
    #   y <- y - delta.y
    # }
    # if (dmod == 2) {
    #   # deflat original x and y
    #   delta.x <- mapply(SIMPLIFY = FALSE, tcrossprod, x=res$pb, y=fac.x)
    #   x <- mapply(SIMPLIFY = FALSE, function(x, y) x-y, x= x, y = delta.x)
    # }
    # if (dmod == 3) {
    #   # deflat crossprod matrix of x and y, using y loading, the same with 1
    #   delta.cmats <- lapply(cmats, function(b, a) b %*% a, a = tcrossprod(t1))
    #   cmats <- mapply(SIMPLIFY = FALSE, function(x, y) x-y, x= cmats, y = delta.cmats)
    # } 
    # if (dmod == 4) {
    #   # deflat crossprod matrix of x and y, using x loading, the same with 2
    #   delta.cmats <- mapply(SIMPLIFY = FALSE, function(b, a) tcrossprod(a) %*% b, b = cmats, a = res$pb)
    #   cmats <- mapply(SIMPLIFY = FALSE, function(x, y) x-y, x= cmats, y = delta.cmats)
    # }
  }
  
  list(loading.x = loading.x, 
       loading.y = loading.y,
       score.x = score.x,
       score.y = score.y, 
       loading.x.index = f.pb, 
       score.x.index = f.tb,
       eig = eig)
}
