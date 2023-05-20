\name{RightSC}
\alias{RightSC}
\title{
clusters nodes in a directed network by regularized spectral clustering on right singular vectors
}
\description{
  community detection by regularized spectral clustering on right singular vectors
}
\usage{
RightSC(A, K, normal = FALSE)
}
\arguments{
  \item{A}{
adjacency matrix of a directed adjacecy matrix
}
  \item{K}{
number of communities
}
  \item{normal}{
 indicator. If TRUE, normalization of singular vector rows will be applied, similar to the spectral spherical clustering.
}
}
\details{
This is essentially the spectral clustering applied on right singular vectors. It can be used to handle directed networks where Aij = 1 if and only if i -> j, and the edges tend to have a missing issue specifically depending on the sender i. More details can be found in Li et. al. (2020).
}


\value{
a list of
 \item{cluster }{cluster labels}
\item{loss }{the loss of Kmeans algorithm}
}
\references{

  T. Li, E. Levina, and J. Zhu. Community models for networks observed through edge nominations. arXiv preprint arXiv:2008.03652 (2020).


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


sc <- RightSC(A,K=2)


}
\keyword{ community detection }
