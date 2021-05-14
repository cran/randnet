\name{NMI}
\alias{NMI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
calculates normalized mutual information
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
  calculates normalized mutual information, a metric that is commonly
  used to compare clustering results
}
\usage{
NMI(g1, g2)
}
%- maybe also 'usage' for other objects documented here.
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
  %%  ~~who you are~~
  Tianxi Li, Elizaveta Levina, Ji Zhu\cr
  Maintainer: Tianxi Li  \email{tianxili@virginia.edu}
}

\examples{

dt <- BlockModel.Gen(30,300,K=3,beta=0.2,rho=0.9,simple=FALSE,power=TRUE)


A <- dt$A


ssc <- reg.SSP(A,K=3,lap=TRUE)

NMI(ssc$cluster,dt$g)


}
