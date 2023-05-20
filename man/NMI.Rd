\name{NMI}
\alias{NMI}
\title{
calculates normalized mutual information
}
\description{
  calculates normalized mutual information, a metric that is commonly
  used to compare clustering results
}
\usage{
NMI(g1, g2)
}
\arguments{
  \item{g1}{
a vector of cluster labels
}
  \item{g2}{
a vector of cluster labels (same length as g1)
}
}

\value{
NMI value
}


\author{
  Tianxi Li, Elizaveta Levina, Ji Zhu\cr
  Maintainer: Tianxi Li  \email{tianxili@virginia.edu}
}

\examples{

dt <- BlockModel.Gen(30,300,K=3,beta=0.2,rho=0.9,simple=FALSE,power=TRUE)


A <- dt$A


ssc <- reg.SSP(A,K=3,lap=TRUE)

NMI(ssc$cluster,dt$g)


}
