\name{LSM.PGD}
\alias{LSM.PGD}
\title{
estimates inner product latent space model by projected gradient descent
}
\description{
  estimates inner product latent space model by projected gradient descent from the paper of Ma et al. (2020).

}
\usage{
LSM.PGD(A, k,step.size=0.3,niter=500,trace=0)
}
\arguments{
  \item{A}{
adjacency matrix
}
  \item{k}{
the dimension of the latent position
}
  \item{step.size}{
step size of gradient descent
}
  \item{niter}{
maximum number of iterations
}
  \item{trace}{
if trace > 0, the objective will be printed out after each iteration
}


}
\details{
The method is based on the gradient descent of Ma et al (2020), with initialization of the universal singular value thresholding as discussed there. The parameter identifiability constraint is the same as in the paper.
}
\value{
a list of
\item{Z }{latent positions}
\item{alpha }{individual parameter alpha as in the paper}
\item{Phat }{esitmated probability matrix}
\item{obj}{the objective of the gradient method}
}
\references{
  Z. Ma, Z. Ma, and H. Yuan. Universal latent space model fitting for large networks with edge
covariates. Journal of Machine Learning Research, 21(4):1-67, 2020.
}
\author{
  Tianxi Li and Can M. Le\cr
  Maintainer: Tianxi Li  \email{tianxili@umn.edu}
}

\seealso{
\code{\link{DCSBM.estimate}}
}
\examples{

dt <- RDPG.Gen(n=600,K=2,directed=TRUE)


A <- dt$A


fit <- LSM.PGD(A,2,niter=50)

}
\keyword{ Latent space model}
