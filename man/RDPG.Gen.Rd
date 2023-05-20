\name{RDPG.Gen}
\alias{RDPG.Gen}
\title{
generates random networks from random dot product graph model
}
\description{
  generates random networks from random dot product graph model

}
\usage{
RDPG.Gen(n, K, directed = TRUE, avg.d = NULL)
}
\arguments{
  \item{n}{
size of the network
}
  \item{K}{
dimension of latent space
}
  \item{directed}{
whether the network is directed or not
}
  \item{avg.d}{
average node degree of the network (in expectation)
}
}
\details{
  The network is generated according to special formulation mentioned in
  ECV paper.
}
\value{
a list of
\item{A }{the adjacency matrix}
\item{P }{the probability matrix}
}
\references{
  S. J. Young and E. R. Scheinerman. Random dot product graph models for
  social networks. In International Workshop on Algorithms and Models for the Web-Graph, pages 138-149. Springer, 2007.
  T. Li, E. Levina, and J. Zhu. Network cross-validation by edge sampling. Biometrika, 107(2), pp.257-276, 2020.
}
\author{
  Tianxi Li, Elizaveta Levina, Ji Zhu\cr
  Maintainer: Tianxi Li  \email{tianxili@virginia.edu}
}

\examples{

dt <- RDPG.Gen(n=600,K=2,directed=TRUE)

A <- dt$A

}
\keyword{ RDPG }
