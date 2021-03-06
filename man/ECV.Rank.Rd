\name{ECV.Rank}
\alias{ECV.Rank}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
estimates optimal low rank model for a network
}
\description{
estimates the optimal low rank model for a network
}
\usage{
ECV.Rank(A, max.K, B = 3, holdout.p = 0.1, weighted = TRUE,mode="directed")
}
%- maybe also 'usage' for other objects documented here.
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
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list of
\item{sse.rank }{rank selection by SSE loss}
\item{auc.rank }{rank selection by AUC loss}
\item{auc }{auc sequence for each rank candidate}
\item{sse }{sse sequence for each rank candidate}
%% ...
}
\references{
  %% ~put references to the literature/web site here ~
  T. Li, E. Levina, and J. Zhu. Network cross-validation by edge sampling. arXiv preprint arXiv:1612.04717, 2016.
}
\author{
  %%  ~~who you are~~
  Tianxi Li, Elizaveta Levina, Ji Zhu\cr
  Maintainer: Tianxi Li  \email{tianxili@umich.edu}
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ rank estimation }
