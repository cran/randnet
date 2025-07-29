\name{ECV.Rank}
\alias{ECV.Rank}
\title{
estimates optimal low rank model for a network
}
\description{
estimates the optimal low rank model for a network
}
\usage{
ECV.Rank(A, max.K, B = 3, holdout.p = 0.1, weighted = TRUE,mode="directed")
}
\arguments{
  \item{A}{
adjacency matrix
}
  \item{max.K}{
maximum possible rank to check
}
  \item{B}{
number of replications in ECV
}
  \item{holdout.p}{
test set proportion
}
  \item{weighted}{
whether the network is weighted. If TRUE, only sum of squared errors
are computed. If FALSE, then treat the network as
binary and AUC will be computed along with SSE.
}
  \item{mode}{
Selectign the mode of "directed" or "undirected" for cross-validation.
}
}
\details{
  AUC is believed to be more accurate in
many simulations for binary networks. But the computation of AUC is much
slower than SSE, even slower than matrix completion steps.

Note that we do not have to assume the true model is low rank. This
function simply finds a best low-rank approximation to the true model.
}
\value{
A list of
\item{sse.rank }{rank selection by SSE loss}
\item{auc.rank }{rank selection by AUC loss}
\item{auc }{auc sequence for each rank candidate}
\item{sse }{sse sequence for each rank candidate}
}
\references{
  T. Li, E. Levina, and J. Zhu. Network cross-validation by edge sampling. Biometrika, 107(2), pp.257-276, 2020.
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


ecv.rank <- ECV.Rank(A,6,weighted=FALSE,mode="undirected")

ecv.rank

}
\keyword{ rank estimation }
