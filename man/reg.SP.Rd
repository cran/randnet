\name{reg.SP}
\alias{reg.SP}
\title{
clusters nodes by regularized spectral clustering
}
\description{
  community detection by regularized spectral clustering
}
\usage{
reg.SP(A, K, tau = 1, lap = FALSE,nstart=30,iter.max=100)
}
\arguments{
  \item{A}{
adjacency matrix
}
  \item{K}{
number of communities
}
  \item{tau}{
reguarlization parameter. Default value is one. Typically set between 0
and 1. If tau=0, no regularization is applied.
}
  \item{lap}{
 indicator. If TRUE, the  Laplacian matrix for clustering. If FALSE, the
 adjacency matrix will be used.
}
  \item{nstart}{
 number of random initializations for K-means
}
  \item{iter.max}{
 maximum number of iterations for K-means
}
}
\details{
The regularlization is done by adding a small constant to each element
of the adjacency matrix. It is shown by such perturbation helps
concentration in sparse networks. It is shown to give consistent
clustering under SBM.
}


\value{
a list of
 \item{cluster }{cluster labels}
\item{loss }{the loss of Kmeans algorithm}
}
\references{
  K. Rohe, S. Chatterjee, and B. Yu. Spectral clustering and the
high-dimensional stochastic blockmodel. The Annals of Statistics, pages
1878-1915, 2011.

  A. A. Amini, A. Chen, P. J. Bickel, and E. Levina. Pseudo-likelihood
  methods for community detection in large sparse networks. The Annals
  of Statistics, 41(4):2097-2122, 2013.

  J. Lei and A. Rinaldo. Consistency of spectral clustering in stochastic block models. The
Annals of Statistics, 43(1):215-237, 2014.

  C. M. Le, E. Levina, and R. Vershynin. Concentration and regularization of random graphs.
Random Structures & Algorithms, 2017.


}
\author{
  Tianxi Li, Elizaveta Levina, Ji Zhu\cr

Maintainer: Tianxi Li <tianxili@virginia.edu>
}

\seealso{
\code{\link{reg.SP}}
}
\examples{


dt <- BlockModel.Gen(30,300,K=3,beta=0.2,rho=0)


A <- dt$A


sc <- reg.SP(A,K=3,lap=TRUE)


NMI(sc$cluster,dt$g)


}
\keyword{ community detection }
