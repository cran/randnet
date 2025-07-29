\name{ConsensusClust}
\alias{ConsensusClust}
\title{
clusters nodes by concensus (majority voting) initialized by regularized spectral clustering
}
\description{
  community detection by concensus (majority voting) initialized by regularized spectral clustering
}
\usage{
ConsensusClust(A,K,tau=0.25,lap=TRUE)
}
\arguments{
  \item{A}{
adjacency matrix
}
  \item{K}{
number of communities
}
  \item{tau}{
reguarlization parameter for regularized spectral clustering. Default value is 0.25. Typically set between 0
and 1. If tau=0, no regularization is applied.
}
  \item{lap}{
 indicator. If TRUE, the  Laplacian matrix for initializing clustering. If FALSE, the
 adjacency matrix will be used.
}
}
\details{
Community detection algorithm by majority voting algorithm of Gao
et. al. (2016). When initialized by regularized spectral clustering, it
is shown that the clustering accuracy of this algorithm gives minimax
rate under the SBM. However, it can slow compared with spectral clustering.
}


\value{
cluster labels
}
\references{
Gao, C.; Ma, Z.; Zhang, A. Y. & Zhou, H. H. Achieving optimal misclassification proportion in stochastic block models The Journal of Machine Learning Research, JMLR. org, 2017, 18, 1980-2024

}
\author{
  Tianxi Li, Elizaveta Levina, Ji Zhu\cr

Maintainer: Tianxi Li <tianxili@umn.edu>
}

\seealso{
\code{\link{reg.SP}}
}
\examples{


dt <- BlockModel.Gen(15,50,K=2,beta=0.2,rho=0)


A <- dt$A


cc <- ConsensusClust(A,K=2,lap=TRUE)


}
\keyword{ community detection }
