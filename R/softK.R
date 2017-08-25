softK <-
function(x, k) {
  if (k <= 0)
    stop("k should be postive integers or (0, 1)")
  n <- length(x)
  if (k >= n) {
    r <- x
  } else {
    if (k < 1)
      k <- round(n * k)
    r <- rep(0, n)
    sx <- sign(x)
    ax <- abs(x)
    cf <- sort(ax, decreasing = TRUE)[k+1]
    r[ax > cf] <- ax[ax > cf] - cf
    r <- sx*r
  }
  return(r)
}
