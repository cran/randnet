\name{network.mixing.Bfold}
\alias{network.mixing.Bfold}
\title{
estimates network connection probability by network mixing with B-fold averaging
}
\description{
  estimates network connection probability by network mixing of Li and Le (2021) with B-fold averaging.
}
\usage{
network.mixing.Bfold(A,B=10,rho = 0.1,max.K=15,dcsbm=TRUE,usvt=TRUE,ns=FALSE,
                       lsm=FALSE,lsm.k=4)
}
\arguments{
  \item{A}{
adjacency matrix
}
  \item{B}{
  number of random replications to average over
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

}
\details{
This is essentially the same procedure as the network.mixing, but repeat it B times and return the average. Use with cautious. Though it can make the estimate more stable, the procedure would increase the computational cost by a factor of B. Based on our limited experience, this is usually not a great trade-off as the improvement might be marginal. 
}
\value{
a list of
\item{linear.Phat }{estimated probability matrix by linear mixing}
\item{nnl.Phat }{estimated probability matrix by NNL mixing}
\item{exp.Phat }{estimated probability matrix by exponential mixing}
\item{ecv.Phat }{estimated probability matrix by ECV mixing (only one nonzero)}
\item{model.names}{the names of all individual models, in the same order as the weights}
}
\references{
  T. Li and C. M. Le, Network Estimation by Mixing: Adaptivity and More. arXiv preprint
arXiv:2106.02803, 2021.

}
\author{
    Tianxi Li and Can M. Le\cr

Maintainer: Tianxi Li <tianxili@virginia.edu>
}
\seealso{
\code{\link{network.mixing}}
}

\examples{


dt <- RDPG.Gen(n=200,K=3,directed=TRUE)
A <- dt$A

fit <- network.mixing.Bfold(A,B=2)
}
