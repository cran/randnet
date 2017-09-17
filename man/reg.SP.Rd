\name{reg.SP}
\alias{reg.SP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
clusters nodes by regularized spectral clustering
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
  community detection by regularized spectral clustering
}
\usage{
reg.SP(A, K, tau = 1, lap = FALSE)
}
%- maybe also 'usage' for other objects documented here.
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
%% ...
}
\references{
  %% ~put references to the literature/web site here ~
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
  %%  ~~who you are~~
  Tianxi Li, Elizaveta Levina, Ji Zhu\cr

Maintainer: Tianxi Li <tianxili@umich.edu>
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ community detection }
