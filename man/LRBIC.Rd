\name{LRBIC}
\alias{LRBIC}
\title{
selecting number of communities by asymptotic likelihood ratio
}
\description{
  selecting number of communities by asymptotic likelihood ratio based
  the methdo of Wang and Bickel 2015
}
\usage{
LRBIC(A, Kmax, lambda = NULL, model = "both")
}
\arguments{
  \item{A}{
adjacency matrix
}
  \item{Kmax}{
the largest possible number of communities to check
}
  \item{lambda}{
a tuning parameter. By default, will use the number recommended in the paper.
}
  \item{model}{
    selecting K under which model. If set to be "SBM", the
    calculation will be done under SBM. If set to be "DCSBM", the
    calculation will be done under DCSBM. The default value is "both" so
will give two selections under SBM and DCSBM respectively.
}
}
\details{
Note that the method cannot distinguish SBM and DCSBM, though different
calculation is done under the two models. So it is not valid to compare
across models. The theoretical analysis of the method is done under
maximum likelhood and variational EM. But as suggested in the paper,
we use spectral clustering for community detection before fitting
maximum likelhood.
}
\value{
a list of
\item{SBM.K }{estimated number of communities under SBM}
\item{DCSBM.K }{estimated number of communities under DCSBM}
\item{SBM.BIC }{the BIC values for the K sequence under SBM}
\item{DCSBM.BIC }{the BIC values for the K sequence under DCSBM}
}
\references{
  Wang, Y. R. & Bickel, P. J. Likelihood-based model selection for stochastic block models The Annals of Statistics, Institute of Mathematical Statistics, 2017, 45, 500-528

}
\author{
  Tianxi Li, Elizaveta Levina, Ji Zhu\cr
  Maintainer: Tianxi Li  \email{tianxili@virginia.edu}
}

\seealso{
 \code{\link{BHMC.estimate}}, \code{\link{ECV.block}}, \code{\link{NCV.select}}
}
\examples{


dt <- BlockModel.Gen(30,300,K=3,beta=0.2,rho=0.9,simple=FALSE,power=TRUE)


A <- dt$A


### test LRBIC

lrbic <- LRBIC(A,6,model="both")

lrbic$SBM.K

lrbic$DCSBM.K

}
\keyword{ community detection }
