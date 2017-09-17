\name{ECV.nSmooth.lowrank}
\alias{ECV.nSmooth.lowrank}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
selecting tuning parameter for neighborhood smoothing estimation of
graphon model
}
\description{
  selecting tuning parameter for neighborhood smoothing estimation of
graphon model where the tuning parameter is to control estimation smoothness.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ECV.nSmooth.lowrank(A, h.seq, K, cv = NULL, B = 3, holdout.p = 0.1)
}
%- maybe also 'usage' for other objects documented here.
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



#ecv.rank <- ECV.Rank(A,10,B=3,weighted=FALSE,mode="undirected")

#K.hat <- ecv.rank$auc.rank ## first estimate a good rank


#h.seq <- sqrt(log(N)/N)*seq(0.5,5,by=0.5)


#ecv.nsmooth <- ECV.nSmooth.lowrank(A,h.seq,K=2,B=3) ## nSmooth can be slow

#h <- h.seq[ecv.nsmooth$min.index]

#What <- nSmooth(A,h=h)

#par(mfrow=c(1,2))
#image(t(W[N:1,]))
#image(t(What[N:1,]))


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ graphon model }
