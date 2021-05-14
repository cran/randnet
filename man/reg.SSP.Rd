\name{reg.SSP}
\alias{reg.SSP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
detects communities by regularized spherical spectral clustering
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
  community detection by regularized spherical spectral clustering
}
\usage{
reg.SSP(A, K, tau = 1, lap = FALSE,nstart=30,iter.max=100)
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
  \item{nstart}{
 number of random initializations for K-means
}
  \item{iter.max}{
 maximum number of iterations for K-means
}
}

\details{
  %%  ~~ If necessary, more details than the description above ~~
  The regularlization is done by adding a small constant to each element
of the adjacency matrix. It is shown by such perturbation helps
concentration in sparse networks. The difference from spectral
clustering (reg.SP) comes from its extra step to normalize the rows of
spectral vectors. It is proved that it gives consistent clustering under DCSBM.
}


\value{
a list of
 \item{cluster }{cluster labels}
\item{loss }{the loss of Kmeans algorithm}
%% ...
}

\references{
  %% ~put references to the literature/web site here ~
    T. Qin and K. Rohe. Regularized spectral clustering under the
  degree-corrected stochastic blockmodel. In Advances in Neural
  Information Processing Systems, pages 3120-3128, 2013.

  J. Lei and A. Rinaldo. Consistency of spectral clustering in stochastic block models. The
Annals of Statistics, 43(1):215-237, 2014.

}
\author{
  %%  ~~who you are~~
    Tianxi Li, Elizaveta Levina, Ji Zhu\cr

Maintainer: Tianxi Li <tianxili@virginia.edu>
}

\seealso{
\code{\link{reg.SP}}
}
\examples{


dt <- BlockModel.Gen(30,300,K=3,beta=0.2,rho=0.9,simple=FALSE,power=TRUE)


A <- dt$A

ssc <- reg.SSP(A,K=3,lap=TRUE)

NMI(ssc$cluster,dt$g)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ community detection }
