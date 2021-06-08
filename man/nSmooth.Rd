\name{nSmooth}
\alias{nSmooth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
estimates probabilty matrix by neighborhood smoothing
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
  estimates probabilty matrix by neighborhood smoothing of Zhang
  et. al. (2017)
}
\usage{
nSmooth(A, h = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
adjacency matrix
}
  \item{h}{
quantile value used for smoothing. Recommended to be in the scale of
sqrt(log(n)/n) where n is the size of the network. The default value
is sqrt(log(n)/n) from the paper.
}
}
\details{
The method assumes a graphon model where the underlying graphon function
is piecewise Lipchitz. However, it may be slow for moderately large
networks, though it is one of the fastest methods for graphon models.
}
\value{
the probability matrix
}
\references{
  %% ~put references to the literature/web site here ~

  Zhang, Y.; Levina, E. & Zhu, J. Estimating network edge probabilities by neighbourhood smoothing Biometrika, Oxford University Press, 2017, 104, 771-783

}
\author{
  %%  ~~who you are~~
    Tianxi Li, Elizaveta Levina, Ji Zhu\cr

Maintainer: Tianxi Li <tianxili@virginia.edu>
}

\examples{


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

#What <- nSmooth(A)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
