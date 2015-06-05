\name{softK}
\alias{softK}

\title{
  Soft-thresholding operator
}
\description{
  Soft-thresholding operator, which is called by \code{\link{mbpca}}.
}
\usage{
softK(x, k)
}

\arguments{
  \item{x}{
  A numerical vector
}
  \item{k}{
  Number of non-zero elements want to keep
}
}

\value{
    A numerical vector
}

\author{
    Chen Meng
}

\examples{
v <- rnorm(10)
softK(v, k = 2)
}

\keyword{soft threshold}