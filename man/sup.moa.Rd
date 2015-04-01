\name{sup.moa}
\alias{sup.moa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Projecting supplementary tables on object of class \code{moa-class}.
}
\description{
    Projecting supplementary tables on \code{\link{moa-class}}
}
\usage{
    sup.moa(X, sup, nf = 2)
}

\arguments{
  \item{X}{
    An object of class \code{\link{moa-class}}
  }
  \item{sup}{
    A list of data.frames contains supplementary data.
  }
  \item{nf}{
    The number of principal components used in the projection.
  }
}

\details{
    Projecting supplementary tables on \code{\link{moa-class}}, for details see reference.
}

\value{
    An object of class \code{\link{moa.sup-class}}. 
}

\references{
    Herve Abdi, Lynne J. Williams, Domininique Valentin and Mohammed Bennani-Dosse. STATIS and DISTATIS: optimum multitable principal component analysis and three way metric multidimensional scaling. WIREs Comput Stat 2012. Volume 4, Issue 2, pages 124-167
}
\author{
    Chen Meng
}

\examples{
    # library(mogsa)
    # loading gene expression data and supplementary data
    data(NCI60_4array_supdata)
    data(NCI60_4arrays)
    # check the dimension of each supplementary data to see how many gene set annotated the data
    sapply(NCI60_4array_supdata, dim)
    # run analysis
    ana <- moa(NCI60_4arrays, proc.row = "center_ssq1", w.data = "inertia", statis = TRUE)
    plot(ana, value="eig")
    # projectin supplementary data
    smoa <- sup.moa(ana, sup=NCI60_4array_supdata, nf=3)
    # heatmap visualize the gene set scores
    heatmap(slot(smoa, "score"))
}

\keyword{ data projection }
\keyword{ supplementary data}
