mbpca <-
function(x, ncomp, method, k="all", center=TRUE, 
                  scale=FALSE, option="uniform", maxiter=1000, 
                  moa=TRUE, verbose=TRUE, svd.solver=c("svd", "fast.svd", "propack")) {
  
  
  x <- lapply(x, t)
  
  call <- match.call()
  x <- processOpt(x, center=center, scale=scale, option = option)
  nc <- sapply(x, ncol)
  keepAll <- k[1] == 'all'
  prddata <- lapply(x, t)
  ssl <- match.arg(svd.solver[1], c("svd", "fast.svd", "propack"))
  svdf <- switch(ssl, 
                 "svd" = svd,
                 "fast.svd" = fast.svd,
                 "propack" = function(X) propack.svd(X, neig = 1))
  
  for (i in 1:ncomp) {
    if (verbose)
      cat(paste("calculating component ", i, " ...\n", sep = ""))
    if (keepAll)
      r <- msvd(x, svd.sol=svdf) else {
        if (length(k) < length(x))
          k <- rep(k, length.out = length(x))
        r <- nipalsSoftK(x, maxiter=maxiter, k=k)
      }
    x <- deflat(x, r$t, r$tb, r$pb, method)
    if (i == 1)
      res <- r else {
        res$t <- cbind(res$t, r$t)
        res$w <- cbind(res$w, r$w)
        res$tb <- mapply(cbind, res$tb, r$tb, SIMPLIFY = FALSE)
        res$pb <- mapply(cbind, res$pb, r$pb, SIMPLIFY = FALSE)
      }
  }
  if (moa) 
    res <- toMoa(prddata, res, call=call)
  return(res)
}
