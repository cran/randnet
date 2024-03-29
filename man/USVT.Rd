\name{USVT}
\alias{USVT}
\title{
estimates the network probability matrix by the improved universal singular value thresholding
}
\description{
  estimates the network probability matrix by the universal singular value thresholding of Chatterjee (2015), with the improvement mentioned in Zhang et. al. (2017).

}
\usage{
USVT(A)
}
\arguments{
  \item{A}{
adjacency matrix
}



}
\details{
Instead of using the original threshold in Chatterjee (2015), the estimate is generated by taking the n^(1/3) leading spectral components. The method was mentioned in Zhang et. al. (2017) and suggested by an anonymous reviewer.
}
\value{
The estimated probability matrix.
%% ...
}
\references{
  %% ~put references to the literature/web site here ~
  S. Chatterjee. Matrix estimation by universal singular value thresholding. The Annals of Statistics,
43(1):177-214, 2015.
  Y. Zhang, E. Levina, and J. Zhu. Estimating network edge probabilities by neighbourhood smoothing. Biometrika, 104(4):771-783, 2017.
}
\author{
  Tianxi Li and Can M. Le\cr
  Maintainer: Tianxi Li  \email{tianxili@virginia.edu}
}

\seealso{
\code{\link{LSM.PGD}}
}
\examples{

dt <- RDPG.Gen(n=600,K=2,directed=TRUE)


A <- dt$A


fit <- USVT(A)

}
\keyword{ singular value thresholding}
