\name{BHMC.estimate}
\alias{BHMC.estimate}
\title{
Estimates the number of communities under block models by the spectral
methods
}
\description{
Estimates the number of communities under block models by using the
spectral properties of network Beth-Hessian matrix with moment correction.

}
\usage{
BHMC.estimate(A, K.max = 15)
}
\arguments{
  \item{A}{
adjacency matrix of the network
}
  \item{K.max}{
the maximum possible number of communities to check
}
}
\details{
  Note that the method cannot distinguish SBM and DCSBM. But it works
  under either model.
}
\value{
A list of result
  \item{K }{Estimated K}
  \item{values }{eigenvalues of the Beth-Hessian matrix}
}
\references{
C. M. Le and E. Levina. Estimating the number of communities in networks by spectral
methods. arXiv preprint arXiv:1507.00827, 2015.
}
\author{
Tianxi Li, Elizaveta Levina, Ji Zhu\cr
Maintainer: Tianxi Li  \email{tianxili@virginia.edu}
}

\seealso{
\code{LRBIC},\code{ECV.Block}, \code{NCV.select}
}
\examples{

dt <- BlockModel.Gen(30,300,K=3,beta=0.2,rho=0.9,simple=FALSE,power=TRUE)


A <- dt$A


bhmc <- BHMC.estimate(A,15)

bhmc

}
\keyword{ community detection }
