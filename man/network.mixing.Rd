\name{network.mixing}
\alias{network.mixing}
\title{
estimates network connection probability by network mixing
}
\description{
  estimates network connection probability by network mixing of Li and Le (2021).
}
\usage{
network.mixing(A,index=NULL,rho = 0.1,max.K=15,dcsbm=TRUE, usvt=TRUE,ns=FALSE,
                           lsm=FALSE,lsm.k=4,trace=FALSE)
}
\arguments{
  \item{A}{
adjacency matrix
}
  \item{index}{
a pre-specified hold-out set. If NULL, the set will be randomly generated according to rho.
}


  \item{rho}{
hold-out proportion as validation entries. Only effective when index is NULL.
}
  \item{max.K}{
the maximum number of blocks used for the block model approximation (see details). 
}
  \item{dcsbm}{
whether to include the DCSBM as components, up to max.K. By default, the method will include it.
}
  \item{usvt}{
whether to include the USVT as a component. By default, the method will include it.
}
  \item{ns}{
whether to include the neighborhood smoothing as a component. 
}

  \item{lsm}{
whether to include the gradient estimator of the latent space model as a component. 
}
  \item{lsm.k}{
the dimension of the latent space. Only effective if lsm is TRUE.
}

  \item{trace}{
whether to print the model fitting progress.
}

}
\details{
The basic version of the mixing estimator will include SBM and DCSBM estimates with the number of blocks from 1 to max.K. Users could also specify whether to include additional USVT, neighborhood smoothing and latent space model estimators. If NNL (non-negative linear), exponential, or ECV is used, the mixing is usually robust for a reasonable range of max.K and whether to include the other models. The linear mixing, however, is vulnerable for a large number of base estimates. The NNL is our recommended method. USVT is also recommended. the neighborhood smoothing and latent space model are slower, so are not suitable for large networks. Details can be found in Li and Le (2021).
}
\value{
a list of
\item{linear.Phat }{estimated probability matrix by linear mixing}
\item{linear.weight }{the weights of the indivdiual models in linear mixing}
\item{nnl.Phat }{estimated probability matrix by NNL mixing}
\item{nnl.weight }{the weights of the indivdiual models in NNL mixing}
\item{exp.Phat }{estimated probability matrix by exponential mixing}
\item{exp.weight }{the weights of the indivdiual models in exponential mixing}
\item{ecv.Phat }{estimated probability matrix by ECV mixing (only one nonzero)}
\item{ecv.weight }{the weights of the indivdiual models in ECV mixing (only one nonzero)}
\item{model.names}{the names of all individual models, in the same order as the weights}
}
\references{
  %% ~put references to the literature/web site here ~

  T. Li and C. M. Le, Network Estimation by Mixing: Adaptivity and More. arXiv preprint
arXiv:2106.02803, 2021.

}
\author{
    Tianxi Li and Can M. Le\cr

Maintainer: Tianxi Li <tianxili@virginia.edu>
}

\examples{


dt <- RDPG.Gen(n=500,K=5,directed=TRUE)

A <- dt$A

fit <- network.mixing(A)
fit$model.names

fit$nnl.weight
}
