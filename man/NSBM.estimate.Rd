\name{NSBM.estimate}
\alias{NSBM.estimate}
\title{
estimates nomination SBM parameters given community labels by the method of moments
}
\description{
  estimates NSBM parameters given community labels

}
\usage{
NSBM.estimate(A,K,g,reg.bound=-Inf)
}
\arguments{
  \item{A}{
adjacency matrix of a directed where Aij = 1 iff i -> j
}
  \item{K}{
number of communities
}
  \item{g}{
a vector of community labels
}
  \item{reg.bound}{
the regularity lower bound of lambda value. By default, -Inf. That means, no constraints. When the network is sparse, using certain constaints may improve stability.
}

}
\details{
The method of moments is used for estimating the edge nomination SBM, so the strategy can be used for both unweighted and weighted networks. The details can be found in Li et. al. (2020).
}
\value{
a list of
\item{B }{estimated block connection probability matrix}
\item{lambda }{estimated lambda values for nomination intensity}
\item{theta }{estimated theta values for nomination preference}
\item{P.tilde }{estimated composiste probability matrix after nomination}
\item{g }{community labels}
}
\references{
  T. Li, E. Levina, and J. Zhu. Community models for networks observed through edge nominations. arXiv preprint arXiv:2008.03652 (2020).
}
\author{
  Tianxi Li, Elizaveta Levina, Ji Zhu\cr
  Maintainer: Tianxi Li  \email{tianxili@virginia.edu}
}

\seealso{
\code{\link{SBM.estimate}}
}
\examples{

dt <- NSBM.Gen(n=200,K=2,beta=0.2,avg.d=10)


A <- dt$A


sc <- RightSC(A,K=3)
est <- NSBM.estimate(A,K=3,g=sc$cluster)

}
\keyword{ NSBM}
