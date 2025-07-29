\name{NCV.select}
\alias{NCV.select}
\title{
selecting block models by NCV
}
\description{
  selecting block models by NCV of Chen and Lei (2016)
}
\usage{
NCV.select(A, max.K, cv = 3)
}
\arguments{
  \item{A}{
adjacency matrix
}
  \item{max.K}{
largest number of communities to check
}
  \item{cv}{
fold of cross-validation
}
}
\details{
Spectral clustering is used for fitting the block models
}
\value{
a list of
\item{dev }{the binomial deviance loss under SBM for each K}
\item{l2 }{the L_2 loss under SBM for each K}
\item{dc.dev }{the binomial deviance loss under DCSBM for each K}
\item{dc.l2 }{the L_2 loss under DCSBM for each K}
\item{dev.model }{the selected model by deviance loss}
\item{l2.model }{the selected model by L_2 loss}
\item{sbm.l2.mat, sbm.dev.mat,....}{the corresponding matrices of loss
  for each fold (row) and each K value (column)}
}
\references{
Chen, K. & Lei, J. Network cross-validation for determining the number of communities in network data Journal of the American Statistical Association, Taylor & Francis, 2018, 113, 241-251

}
\author{
  Tianxi Li, Elizaveta Levina, Ji Zhu\cr
  Maintainer: Tianxi Li  \email{tianxili@umn.edu}
}

\seealso{
 \code{\link{ECV.block}}
}
\examples{

dt <- BlockModel.Gen(30,300,K=3,beta=0.2,rho=0.9,simple=FALSE,power=TRUE)


A <- dt$A


ncv <- NCV.select(A,6,3)

ncv$l2.model
ncv$dev.model

which.min(ncv$dev)
which.min(ncv$l2)

which.min(ncv$dc.dev)
which.min(ncv$dc.l2)

}
\keyword{ community detection }
