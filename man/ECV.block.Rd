\name{ECV.block}
\alias{ECV.block}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
selecting block models by ECV
}
\description{
Model selection by ECV for SBM and DCSBM. It can be used to select
between the two models or given on model (either SBM or DCSBM) and
select K.
}
\usage{
ECV.block(A, max.K, cv = NULL, B = 3, holdout.p = 0.1, tau = 0, dc.est = 2, kappa = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
adjacency matrix
}
  \item{max.K}{
largest possible K for number of communities
}
  \item{cv}{
cross validation fold. The default value is NULL. We recommend to use
the argument B instead, doing indpendent sampling.
}
  \item{B}{
number of replications
}
  \item{holdout.p}{
testing set proportion
}
  \item{tau}{
constant for numerical stability only. Not useful for current version.
}
  \item{dc.est}{
estimation method for DCSBM. By defaulty (dc.est=2), the maximum
likelihood is used. If dc.est=1, the method used by Chen and Lei (2016)
is used, which is less stable according to our observation.
}
  \item{kappa}{
    %%     ~~Describe \code{kappa} here~~
    constant for numerical stability only. Not useful for current version.
}
}
\details{
  %%  ~~ If necessary, more details than the description above ~~
  The current version is based on a simple matrix completion procedure,
  as described in the paper. The performance can be improved by better
  matrix completion method that will be implemented in next
  version. Moreover, the current implementation is better in
  computational time but less efficient in memory. Another memory
  efficient implementation will be added in next version.
}
\value{
\item{impute.err }{average validaiton imputation error}
\item{l2 }{average validation L_2 loss under SBM}
\item{dev }{average validation binomial deviance loss under SBM}
\item{auc }{average validation AUC}
\item{dc.l2}{average validation L_2 loss under DCSBM}
\item{dc.dev}{average validation binomial deviance loss under DCSBM}
\item{sse}{average validation SSE}
\item{l2.model }{selected model by L_2 loss}
\item{dev.model }{selected model by binomial deviance loss}
\item{l2.mat, dc.l2.mat,... }{cross-validation loss matrix for B replications}
%% ...
}
\references{
  %% ~put references to the literature/web site here ~
  T. Li, E. Levina, and J. Zhu. Network cross-validation by edge sampling. Biometrika, 107(2), pp.257-276, 2020.
}
\author{
  %%  ~~who you are~~
  Tianxi Li, Elizaveta Levina, Ji Zhu\cr
  Maintainer: Tianxi Li  \email{tianxili@virginia.edu}
}

\seealso{
\code{NCV.select}
}
\examples{

dt <- BlockModel.Gen(30,300,K=3,beta=0.2,rho=0.9,simple=FALSE,power=TRUE)


A <- dt$A

ecv <- ECV.block(A,6,B=3)

ecv$l2.model
ecv$dev.model


which.min(ecv$l2)
which.min(ecv$dev)

which.min(ecv$dc.l2)
which.min(ecv$dc.dev)

which.max(ecv$auc)
which.min(ecv$sse)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ECV }
