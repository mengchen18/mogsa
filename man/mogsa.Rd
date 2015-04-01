\name{mogsa}
\alias{mogsa}

\title{
    multiple omics data integration and gene set analysis
}

\description{
    The main function called by users, omics data analysis and gene set annotation. 
    A wrapper function of \code{\link{moa}} and \code{\link{sup.moa}}.
}

\usage{
  mogsa(x, sup, nf=NULL, proc.row=NULL, w.data=NULL, w.row=NULL, statis=FALSE)
}

\arguments{
    \item{x}{
      An object of class \code{list} or \code{\link{moa-class}}. A list would be a list of data frame. 
    }
    \item{sup}{
      An object of class \code{list} or \code{\link{moa.sup-class}}. A list would be a list of supplementary data.
    }
    \item{nf}{
      The number of principal components used to reconstruct, only used when x is a an object of \code{list}.
    }
    \item{proc.row}{
      Preprocessing of rows. If x is a object of \code{list}, it is passed \code{moa}
    }
    \item{w.data}{
      Weights of datasets. If x is a object of \code{list}, it is passed \code{moa}
    }
    \item{w.row}{
      Weight of row. If x is a object of \code{list}, it is passed \code{moa}
    }
    \item{statis}{
      A logical indicates if statis algrithm should be used. If x is a object of \code{list}, it is passed \code{moa}
    }
}

\details{
    A wrapper function of \code{\link{moa}} and \code{\link{sup.moa}}. 
}

\value{
    An object of class \code{\link{mgsa-class}}. 
}

\references{
    Not avaiable right now.
}

\author{
    Chen Meng 
}

\note{
  This function will be changed to a generic function for "S4-style" programming.
}

\seealso{
    \code{\link{moa}} and \code{\link{sup.moa}}
}

\examples{
  # library(mogsa)
  # loading gene expression data and supplementary data
  data(NCI60_4array_supdata)
  data(NCI60_4arrays)

  # using a list of data.frame as input
  mgsa1 <- mogsa(x = NCI60_4arrays, sup=NCI60_4array_supdata, nf=9,
                 proc.row = "center_ssq1", w.data = "inertia", statis = TRUE)
  # using moa as input
  ana <- moa(NCI60_4arrays, proc.row = "center_ssq1", w.data = "inertia", statis = TRUE)
  smoa <- sup.moa(ana, sup=NCI60_4array_supdata, nf=3)
  mgsa2 <- mogsa(x = ana, sup=NCI60_4array_supdata, nf=9)
  mgsa3 <- mogsa(x = ana, sup=smoa)
}

\keyword{ MVA }
\keyword{ supplementary data projection }
