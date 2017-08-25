
biSoftK <- function (x, maxiter, kp, kt, unit.pb = TRUE, unit.tb = FALSE) {
  
  if (length(kp) < length(x)) 
    kp <- rep(kp, length.out = length(x))
  if (length(kt) < length(x))
    kt <- rep(kt, length.out = length(x))
  
  t <- svd(do.call("cbind", x))$u[, 1]
  
  regproj <- function(xb, t, kp, kt, unit.pb, unit.tb) {
    pb <- t(xb) %*% t/c(t(t) %*% t)
    if (!unit.pb)
      pb <- normvec(pb)
    pb <- softK(pb, kp)
    if (unit.pb)
      pb <- normvec(pb)
    tb <- xb %*% pb
    tb <- softK(tb, kt)
    if (unit.tb)
      tb <- normvec(tb)
    list(tb = tb, pb = pb)
  }
  
  for (i in 1:maxiter) { 
    told <- t
    rp <- mapply(SIMPLIFY = FALSE, function(x, y1, y2) 
      regproj(x, t, y1, y2, unit.pb = unit.pb, unit.tb = unit.tb), 
      x = x, y1 = kp, y2 = kt)
    tm <- sapply(rp, "[[", "tb")
    w <- t(tm) %*% t/c(t(t) %*% t)
    w <- w/sqrt(sum(w^2))
    t <- tm %*% w
    if (isTRUE(all.equal(c(t), c(told)))) 
      break
    if (i == maxiter) 
      cat("  Note: maximum number of iterations was reached, algrithm may not converge.\n")
  }
  res <- list(tb = lapply(rp, "[[", "tb"), pb = lapply(rp, "[[", "pb"), t = t, w = w)
  return(res)
  
}