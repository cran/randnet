\name{NSBM.Gen}
\alias{NSBM.Gen}
\title{
Generates networks from nomination stochastic block model
}
\description{
Generates networks from nomination stochastic block model for community structure in edge nomination procedures, proposed in Li et. al. (2020)
}
\usage{
NSBM.Gen( n, K, avg.d,beta,theta.low=0.1,
    theta.p=0.2,lambda.scale=0.2,lambda.exp=FALSE)
}
\arguments{
  \item{n}{
 size of network
}
  \item{K}{
number of communities
}
  \item{avg.d}{
expected average degree of the resuling network (after edge nomination)
}
  \item{beta}{
the out-in ratio of the original SBM
}
  \item{theta.low}{
the lower value of theta's. The theta's are generated as two-point mass at theta.low and 1.
}
  \item{theta.p}{
proportion of lower value of theta's
}
  \item{lambda.scale}{
standard deviation of the lambda (before the exponential, see lambda.exp)
}
  \item{lambda.exp}{
If TRUE, lambda is generated as exponential of uniformation random randomes. Otherwise, they are normally distributed.
}
}

\value{
A list of
 \item{A}{the generated network adjacency matrix}
 \item{g }{community membership}
 \item{P }{probability matrix of the orignal SBM network}
 \item{P.tilde }{probability matrix of the observed network after nomination}
 \item{B }{B parameter}
 \item{lambda }{lambda parameter}
 \item{theta }{theta parameter}
}
\references{
  T. Li, E. Levina, and J. Zhu. Community models for networks observed through edge nominations. arXiv preprint arXiv:2008.03652 (2020).


}
\author{
  Tianxi Li, Elizaveta Levina, Ji Zhu\cr
  Maintainer: Tianxi Li  \email{tianxili@virginia.edu}
}

\examples{

dt <- NSBM.Gen(n=200,K=2,beta=0.2,avg.d=10)



}
\keyword{ BlockModel.Gen }
