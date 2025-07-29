\name{k.core}
\alias{k.core}
\title{
identify the K-core component of a network
}
\description{
identify the K-core component of a network
}
\usage{
k.core(A,K)
}
\arguments{
  \item{A}{
adjacency matrix. 
}
  \item{K}{
lower bound of the degree.
}

}
\details{
  The function can be used as a general data processing function to extrat denser component of a network. The algorithm is iteratively filter components with degree below K, following Wang and Rohe (2016).
}
\value{
A list of
\item{A }{The resulting adjacency matrix.}
}
\references{
  Wang S. and Rohe, K, Discussion of "coauthorship and citation networks for statisticians". The Annals of Applied Statistics, 10(4):1820-1826, 2016.

}
\author{
  Tianxi Li, Elizaveta Levina, Ji Zhu, Can M. Le\cr
  Maintainer: Tianxi Li  \email{tianxili@umn.edu}
}

\seealso{
 \code{\link{InformativeCore}}
}
\examples{
set.seed(100)
dt <- BlockModel.Gen(60,1000,K=3,beta=0.2,rho=0.9,simple=FALSE,power=TRUE)
### this is not an interesting case -- only for demonstration of the usage. 
### The network has no periphery nodes, all nodes are in the core.

A <- dt$A


kc <- k.core(A,K=3)
dim(kc)

}
\keyword{ core-periphery}
