\name{SBM.estimate}
\alias{SBM.estimate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
estimates SBM parameters given community labels
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
  estimates SBM parameters given community labels

}
\usage{
SBM.estimate(A, g)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
adjacency matrix
}
  \item{g}{
a vector of community labels
}
}
\details{
maximum likelhood is used
}
\value{
a list of
\item{B }{estimated block connection probability matrix}
\item{Phat }{estimated probability matrix}
\item{g }{community labels}
%% ...
}
\references{
  %% ~put references to the literature/web site here ~
  B. Karrer and M. E. Newman. Stochastic blockmodels and community structure in networks.
Physical Review E, 83(1):016107, 2011.
}
\author{
  %%  ~~who you are~~
  Tianxi Li, Elizaveta Levina, Ji Zhu\cr
  Maintainer: Tianxi Li  \email{tianxili@umich.edu}
}

\seealso{
\code{\link{DCSBM.estimate}}
}
\examples{

dt <- BlockModel.Gen(30,300,K=3,beta=0.2,rho=0)


A <- dt$A


sc <- reg.SP(A,K=3,lap=TRUE)
sbm <- SBM.estimate(A,sc$cluster)
sbm$B

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ SBM}
