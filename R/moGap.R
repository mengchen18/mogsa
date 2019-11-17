#' Gap statistic for clustering latent variables in \code{\link{moa-class}}.
#' 
#' Gap statitistic is a measurement of goodness of clustering result. This is a
#' convenient function to calculate the gap statistic of clustering "moa".
#' 
#' 
#' @param x An object of class \code{\link{moa-class}} returned by
#' \code{\link{mbpca}}.
#' @param K.max The maximum number of clusters to consider, passed to
#' \code{clusGap}
#' @param B The number of bootstrap, passed to \code{clusGap}
#' @param cluster A charater string could be either "kmeans" or "hclust" to
#' specify the clustering algorithm.
#' @param plot Logical; whether return the gap statistic plot.
#' @param dist.method Distance meaurement, passed to function \code{"dist"}.
#' @param dist.diag Passed to function \code{"dist"}.
#' @param dist.upper Passed to function \code{"dist"}.
#' @param dist.p Passed to function \code{"dist"}.
#' @param hcl.method Hierarchical clustering method, passed to \code{"hclust"}
#' @param hcl.members Passed to \code{"hclust"}
#' @param km.iter.max Maximum number of iteration in kmeans, passed to
#' \code{"kmeans"}.
#' @param km.nstart An integer to specify how many random sets should be
#' chosen. passed to \code{"kmeans"}.
#' @param km.algorithm Kmeans algorithm, passed to \code{"kmeans"}.
#' @param km.trace See function \code{"kmeans"}.
#' @return It returns a list consists of five components:
#' 
#' "Tab", "n", "B", "FUNcluster" - see \code{clusGap}
#' 
#' "nClust" - the estimated number of clusters using different method, see
#' \code{maxSE}
#' @author Chen Meng
#' @seealso Function "clusGap" in "cluster" package Function "dist", "hclust",
#' "kmeans"
#' @references Tibshirani, R., Walther, G. and Hastie, T. (2001).  Estimating
#' the number of data clusters via the Gap statistic.  Journal of the Royal
#' Statistical Society B, 63, 411-423.
#' 
#' Maechler, M., Rousseeuw, P., Struyf, A., Hubert, M., Hornik, K.(2015).
#' cluster: Cluster Analysis Basics and Extensions.  R package version 2.0.1.
#' @keywords gap statistic moa
#' @importFrom cluster clusGap
#' @export
#' @examples
#' 
#' # see examples in \code{\link{mbpca}}
#' 
#' 
#' data("NCI60_4arrays")
#' moa <- mbpca(NCI60_4arrays, ncomp = 10, k = "all", method = "globalScore", option = "lambda1", 
#'              center=TRUE, scale=FALSE)
#' gap <- moGap(moa, K.max = 12, cluster = "hcl")
#' 
#' genes <- moaCoef(moa)
#' scr <- moaScore(moa)
#' 
#' moa2 <- moa(NCI60_4arrays, proc.row="center_ssq1", w.data="inertia", w.row=NULL, statis=FALSE)
#' gap2 <- moGap(moa, K.max = 12, cluster = "hcl")
#' 
moGap <-
function(x, K.max, B=100, cluster=c("kmeans", "hclust"), plot=TRUE,
                  dist.method = "euclidean", dist.diag = FALSE, dist.upper = FALSE, dist.p = 2,
                  hcl.method = "complete", hcl.members = NULL,
                  km.iter.max = 10, km.nstart = 10, 
                  km.algorithm = c("Hartigan-Wong", "Lloyd", "Forgy", "MacQueen"), km.trace=FALSE) {
  
  cluster <- cluster[1]
  sr <- moaScore(x)
  fhclust <- function(x, k,
                      dist.method = "euclidean", dist.diag = FALSE, dist.upper = FALSE, dist.p = 2,
                      hcl.method = "complete", hcl.members = NULL) {
    d <- dist(x, method=dist.method, diag=dist.diag, upper=dist.upper, p=dist.p)
    hl <- hclust(d, method=hcl.method, members=hcl.members)
    cls <- cutree(hl, k=k)
    list(cluster=cls)
  }
  
  if (pmatch(cluster, "kmeans", nomatch = 0))
    v <- clusGap(sr, FUNcluster = kmeans, K.max = K.max, B = B, 
                iter.max=km.iter.max, nstart=km.nstart, algorithm=km.algorithm, trace=km.trace) else
      if (pmatch(cluster, "hclust", nomatch = 0))
        v <- clusGap(sr, FUNcluster = fhclust, K.max = K.max, B = B, 
                     dist.method = dist.method, dist.diag = dist.diag, dist.upper = dist.upper, dist.p = dist.p,
                     hcl.method = hcl.method, hcl.members = hcl.members) else
                       stop("Unkown clustering algorithm specified!")
  if (plot)
    plot(v)
  n <- sapply(c("firstSEmax", "Tibs2001SEmax", "globalSEmax","firstmax", "globalmax"), 
              function(x) maxSE(v$Tab[, 3], v$Tab[, 4], method = x))
  v$nClust <- n
  v
}
