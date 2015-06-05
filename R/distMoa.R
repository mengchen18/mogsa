distMoa <-
function(x, nf=Inf, tol=1e-5, method = "euclidean", diag = FALSE, upper = FALSE, p = 2) {
  if (inherits(x, "moa")) {
    nfi <- x@eig > tol
    x <- moaScore(x)
  }
  if (nf > ncol(x) | nf > length(nfi)) {
    nf <- min(ncol(x), length(nfi))
    cat(paste("nf set to ", nf, ".\n", sep = ""))
  }
  nfi <- nfi[1:nf]
  x <- x[, nfi]
  dist(x, method = method, diag = diag, upper = upper, p = p)
}
