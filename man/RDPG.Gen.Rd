\name{RDPG.Gen}
\alias{RDPG.Gen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
generates random networks from random dot product graph model
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
  generates random networks from random dot product graph model

}
\usage{
RDPG.Gen(n, K, directed = TRUE, avg.d = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
size of the network
}
  \item{K}{
dimension of latent space
}
  \item{directed}{
whether the network is directed or not
}
  \item{avg.d}{
average node degree of the network (in expectation)
}
}
\details{
  %%  ~~ If necessary, more details than the description above ~~
  The network is generated according to special formulation mentioned in
  ECV paper.
}
\value{
a list of
\item{A }{the adjacency matrix}
\item{P }{the probability matrix}
%% ...
}
\references{
  %% ~put references to the literature/web site here ~
  S. J. Young and E. R. Scheinerman. Random dot product graph models for
  social networks. In International Workshop on Algorithms and Models for the Web-Graph, pages 138-149. Springer, 2007.
  T. Li, E. Levina, and J. Zhu. Network cross-validation by edge sampling. arXiv preprint arXiv:1612.04717, 2016.
}
\author{
  %%  ~~who you are~~
  Tianxi Li, Elizaveta Levina, Ji Zhu\cr
  Maintainer: Tianxi Li  \email{tianxili@umich.edu}
}

\examples{

dt <- RDPG.Gen(n=600,K=2,directed=TRUE)

A <- dt$A

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ RDPG }
