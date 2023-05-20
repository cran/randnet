\name{SBM.estimate}
\alias{SBM.estimate}
\title{
estimates SBM parameters given community labels
}
\description{
  estimates SBM parameters given community labels

}
\usage{
SBM.estimate(A, g)
}
\arguments{
  \item{A}{
adjacency matrix
}
  \item{g}{
a vector of community labels
}
}
\details{
maximum likelhood is used
}
\value{
a list of
\item{B }{estimated block connection probability matrix}
\item{Phat }{estimated probability matrix}
\item{g }{community labels}
}
\references{
  B. Karrer and M. E. Newman. Stochastic blockmodels and community structure in networks.
Physical Review E, 83(1):016107, 2011.
}
\author{
  Tianxi Li, Elizaveta Levina, Ji Zhu\cr
  Maintainer: Tianxi Li  \email{tianxili@virginia.edu}
}

\seealso{
\code{\link{DCSBM.estimate}}
}
\examples{

dt <- BlockModel.Gen(30,300,K=3,beta=0.2,rho=0)


A <- dt$A


sc <- reg.SP(A,K=3,lap=TRUE)
sbm <- SBM.estimate(A,sc$cluster)
sbm$B

}
\keyword{ SBM}
