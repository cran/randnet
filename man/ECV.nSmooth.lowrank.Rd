\name{ECV.nSmooth.lowrank}
\alias{ECV.nSmooth.lowrank}
\title{
selecting tuning parameter for neighborhood smoothing estimation of
graphon model
}
\description{
  selecting tuning parameter for neighborhood smoothing estimation of
graphon model where the tuning parameter is to control estimation smoothness.
}
\usage{
ECV.nSmooth.lowrank(A, h.seq, K, cv = NULL, B = 3, holdout.p = 0.1)
}
\arguments{
  \item{A}{
adjacency matrix
}
  \item{h.seq}{
a sequence of h values to tune. It is suggested h should be in the order
of sqrt(log(n)/n).
}
  \item{K}{
the optimal rank for approximation. Can be obtained by rank selection of ECV.
}
  \item{cv}{
cross-validation fold. Recomend to use replication number B instead.
}
  \item{B}{
independent replication number of random splitting
}
  \item{holdout.p}{
proportion of test sample
}
}
\details{
The neighborhood smoothing estimation can be slow, so the ECV may take
long even for moderately large networks.
}
\value{
a list object with
\item{err }{average validation error for h.seq}
\item{min.index }{index of the minimum error}
}
\references{
  T. Li, E. Levina, and J. Zhu. Network cross-validation by edge sampling. Biometrika, 107(2), pp.257-276, 2020.

}
\author{
  Tianxi Li, Elizaveta Levina, Ji Zhu\cr
  Maintainer: Tianxi Li  \email{tianxili@virginia.edu}
}

\examples{

set.seed(500)
N <- 300

U = matrix(1:N,nrow=1) / (N+1)
V = matrix(1:N,nrow=1) / (N+1)

W = (t(U))^2%*%matrix(1,nrow=1,ncol=N)+matrix(1,nrow=N,ncol=1)%*%(V^2)
W = W/3*cos(1/(W + 1e-7)) + 0.15



upper.index <- which(upper.tri(W))

A <- matrix(0,N,N)


rand.ind <- runif(length(upper.index))

edge.index <- upper.index[rand.ind < W[upper.index]]

A[edge.index] <- 1

A <- A + t(A)
diag(A) <- 0

\donttest{
h.seq <- sqrt(log(N)/N)*seq(0.5,5,by=0.5)


ecv.nsmooth <- ECV.nSmooth.lowrank(A,h.seq,K=2,B=3) 

h <- h.seq[ecv.nsmooth$min.index]
}

}
\keyword{ graphon model }
