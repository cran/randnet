\name{BlockModel.Gen}
\alias{BlockModel.Gen}
\title{
Generates networks from degree corrected stochastic block model
}
\description{
Generates networks from degree corrected stochastic block model, with
various options for node degree distribution
}
\usage{
BlockModel.Gen(lambda, n, beta = 0, K = 3, w = rep(1, K),
 Pi = rep(1, K)/K, rho = 0, simple = TRUE, power = TRUE,
alpha = 5, degree.seed = NULL)
}
\arguments{
  \item{lambda}{
average node degree
}
  \item{n}{
 size of network
}
  \item{beta}{
out-in ratio: the ratio of between-block edges over within-block edges
}
  \item{K}{
number of communities
}
  \item{w}{
not effective
}
  \item{Pi}{
a vector of community proportion
}
  \item{rho}{
proportion of small degrees within each community if the degrees are
from two point mass disbribution. rho >0 gives degree corrected block
model. If rho > 0 and simple=TRUE, then generate the degrees from two
point mass distribution, with rho porition of 0.2 values and 1-rho
proportion of 1 for degree parameters. If rho=0, generate from SBM.
}
  \item{simple}{
Indicator of wether two point mass degrees are used, if rho > 0. If
rho=0, this is not effective
}
  \item{power}{
Whether or not use powerlaw distribution for degrees. If FALSE, generate
from theta from U(0.2,1); if TRUE, generate theta from powerlaw. Only
effective if rho >0, simple=FALSE.
}
  \item{alpha}{
Shape parameter for powerlaw distribution.
}
  \item{degree.seed}{
Can be a vector of a prespecified values for theta. Then the function
will do sampling with replacement from the vector to generate theta. It
can be used to control noise level between different configuration settings.
}
}

\value{
A list of
 \item{A}{the generated network adjacency matrix}
 \item{g }{community membership}
 \item{P }{probability matrix of the network}
 \item{theta }{node degree parameter}
}
\references{
  B. Karrer and M. E. Newman. Stochastic blockmodels and community structure in networks.
  Physical Review E, 83(1):016107, 2011.

A. A. Amini, A. Chen, P. J. Bickel, and E. Levina. Pseudo-likelihood
methods for community detection in large sparse networks. The Annals of
Statistics, 41(4):2097-2122, 2013.

T. Li, E. Levina, and J. Zhu. Network cross-validation by edge sampling. Biometrika, 107(2), pp.257-276, 2020.

}
\author{
  Tianxi Li, Elizaveta Levina, Ji Zhu\cr
  Maintainer: Tianxi Li  \email{tianxili@virginia.edu}
}

\examples{

dt <- BlockModel.Gen(30,300,K=3,beta=0.2,rho=0.9,simple=FALSE,power=TRUE)



}
\keyword{ SBM }
