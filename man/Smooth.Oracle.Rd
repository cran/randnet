\name{smooth.oracle}
\alias{smooth.oracle}
\title{
oracle smooth graphon estimation
}
\description{
oracle smooth graphon estimation according to given latent positions, based on smooth estimation.
}
\usage{
smooth.oracle(Us,A)
}
\arguments{
  \item{Us}{
a vector whose elements are the latent positions of the network nodes under the graphon model. The dimension of the vector equals the number of nodes in the network.
}
  \item{A}{
adjacency matrix. It does not have to be unweighted. 
}

}
\details{
  Note that the latenet positions are unknown in practice, so this estimation is an oracle estimation mainly for evaluation purpose. However, if the latenet positions can be approximated estimated, this function can also be used for estimating the edge probability matrix. This procedure is the M-step of the algorithm used in Sischka & Kauermann (2022). Our implementation is modified from the contribution of an anonymous reviewer, leveraging the tools of the sparseFLMM package.
}
\value{
The estimated probability matrix.
}
\references{
  Sischka, B. and Kauermann, G., 2022. EM-based smooth graphon estimation using MCMC and spline-based approaches. Social Networks, 68, pp.279-295.
}
\author{
  Tianxi Li, Elizaveta Levina, Ji Zhu, Can M. Le\cr
  Maintainer: Tianxi Li  \email{tianxili@umn.edu}
}


\examples{
set.seed(100)
dt <- BlockModel.Gen(10,50,K=2,beta=0.2)

## oracle order
oracle.index <- sort(dt$g,index.return=TRUE)$ix
A <- dt$A[oracle.index,oracle.index]

oracle.est <- smooth.oracle(seq(0,1,length.out=50),A)

}
\keyword{ smooth graphon}
