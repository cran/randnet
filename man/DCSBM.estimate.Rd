\name{DCSBM.estimate}
\alias{DCSBM.estimate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimates DCSBM model
}
\description{
  Estimates DCSBM model by given community labels

}
\usage{
DCSBM.estimate(A, g)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
adjacency matrix
}
  \item{g}{
vector of community labels for the nodes
}
}
\details{
Estimation is based on maximum likelhood.
}
\value{
A list object of
\item{Phat }{estimated probability matrix}
\item{B }{the B matrix with block connection probability, up to a
  scaling constant}
\item{Psi }{vector of of degree parameter theta, up to a scaling constant}
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
  Maintainer: Tianxi Li  \email{tianxili@virginia.edu}
}

\seealso{
\code{SBM.estimate}
}
\examples{

dt <- BlockModel.Gen(30,300,K=3,beta=0.2,rho=0.9,simple=FALSE,power=TRUE)


A <- dt$A


ssc <- reg.SSP(A,K=3,lap=TRUE)

est <- DCSBM.estimate(A,ssc$cluster)

  }

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ DCSBM }
