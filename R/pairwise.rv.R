#' pairwise RV coefficients.
#' 
#' Calculating pairwise RV coefficients for a list of matrices or data.frame.
#' 
#' The RV coefficient for each pair of matrices is calculated as Rv =
#' trace(XX'YY')/sqrt(trace(XX'XX')*trace(YY'YY'))
#' 
#' @param data.list A list of data.frame or matrix, either rows or columns in
#' each data set should be matched.
#' @param match Whether columns or rows of data.frame/matrix should be matched.
#' @return The function will return a matrix containing the pairwise RV
#' coefficients.
#' @note The variable in matrices are not automatically centered or scaled in
#' this function. So these step may need to be performed before calling this
#' function.
#' @author Chen Meng
#' @references Robert, P.; Escoufier, Y. (1976). A Unifying Tool for Linear
#' Multivariate Statistical Methods: The RV-Coefficient. Applied Statistics 25
#' (3): 257-265.
#' @keywords RV coefficent
#' @examples
#' 
#'     data(NCI60_4arrays)
#'     pairwise.rv(NCI60_4arrays)
#' 
pairwise.rv <- function(data.list, match="col") {
	if (match %in% "row")
		data.list <- lapply(data.list, t)
  ms <- sapply(data.list, function(x) {
    x <- c(crossprod(as.matrix(x)))
    x <- x/sqrt(sum(x^2))})
  m <- crossprod(ms)
  colnames(m) <- rownames(m) <- names(data.list)
  return(m)
}
