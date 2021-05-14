\name{NCV.select}
\alias{NCV.select}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
selecting block models by NCV
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
  selecting block models by NCV of Chen and Lei (2016)
}
\usage{
NCV.select(A, max.K, cv = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
adjacency matrix
}
  \item{max.K}{
largest number of communities to check
}
  \item{cv}{
fold of cross-validation
}
}
\details{
Spectral clustering is used for fitting the block models
}
\value{
a list of
\item{dev }{the binomial deviance loss under SBM for each K}
\item{l2 }{the L_2 loss under SBM for each K}
\item{dc.dev }{the binomial deviance loss under DCSBM for each K}
\item{dc.l2 }{the L_2 loss under DCSBM for each K}
\item{dev.model }{the selected model by deviance loss}
\item{l2.model }{the selected model by L_2 loss}
\item{sbm.l2.mat, sbm.dev.mat,....}{the corresponding matrices of loss
  for each fold (row) and each K value (column)}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
  %% ~put references to the literature/web site here ~
Chen, K. & Lei, J. Network cross-validation for determining the number of communities in network data Journal of the American Statistical Association, Taylor & Francis, 2018, 113, 241-251

}
\author{
  %%  ~~who you are~~
  Tianxi Li, Elizaveta Levina, Ji Zhu\cr
  Maintainer: Tianxi Li  \email{tianxili@virginia.edu}
}

\seealso{
 \code{\link{ECV.block}}
}
\examples{

dt <- BlockModel.Gen(30,300,K=3,beta=0.2,rho=0.9,simple=FALSE,power=TRUE)


A <- dt$A


ncv <- NCV.select(A,6,3)

ncv$l2.model
ncv$dev.model

which.min(ncv$dev)
which.min(ncv$l2)

which.min(ncv$dc.dev)
which.min(ncv$dc.l2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ community detection }
