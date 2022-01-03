\name{InformativeCore}
\alias{InformativeCore}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
identify the informative core component of a network
}
\description{
identify the informative core component of a network based on the spectral method of Miao and Li (2021). It can be used as a general data processing function for any network modeling purpose. 
}
\usage{
InformativeCore(A,r=3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
adjacency matrix. It does not have to be unweighted. 
}
  \item{r}{
the rank for low-rank denoising. The rank can be selected by ECV or any other methods availale in the package.
}

}
\details{
  The function can be used as a general data processing function for any network modeling purpose. It automatically identify an informative core component with interesting connection structure and a noninformative periphery component with uninterestings structures. Depending on the user's preference, the uninteresting structure can be either the Erdos-Renyi type connections or configuration type of connections, both of which are generally regarded as noninformative structures. Including these additional non-informative structures in network models can potentially lower the modeling efficiency. Therefore, it is preferable to remove them and only focus on the core structure. Details can be found in the reference.
}
\value{
A list of
\item{er.score }{A n dimensional vector of informative scores for ER-type periphery. A larger score indicates the corresponding node is more likely to be in the core.}
\item{config.score }{A n dimensional vector of informative scores for configuration-type periphery. A larger score indicates the corresponding node is more likely to be in the core.}
\item{er.theory.core }{The indices of identified core structure in the ER-type model based on a theoretical threshold of the scores (for large sample size). }
\item{config.theory.core }{The indices of identified core structure in the configuration-type model based on a theoretical threshold of the scores (for large sample size). }
\item{er.kmeans.core }{The indices of identified core structure in the ER-type model based on kmeans clustering of the scores. }
\item{config.kmeans.core }{The indices of identified core structure in the configuration-type model based on kmeans clustering of the scores (for large sample size). }
%% ...
}
\references{
  %% ~put references to the literature/web site here ~
  R. Miao and T. Li. Informative core identification in complex networks. arXiv preprint arXiv:2101.06388, 2021

}
\author{
  %%  ~~who you are~~
  Tianxi Li, Elizaveta Levina, Ji Zhu, Can M. Le\cr
  Maintainer: Tianxi Li  \email{tianxili@virginia.edu}
}

\seealso{
 \code{\link{ECV.Rank}}
}
\examples{
set.seed(100)
dt <- BlockModel.Gen(60,1000,K=3,beta=0.2,rho=0.9,simple=FALSE,power=TRUE)
### this is not an interesting case -- only for demonstration of the usage. 
### The network has no periphery nodes, all nodes are in the core.

A <- dt$A


core.fit <- InformativeCore(A,r=3)
length(core.fit$er.theory.core)
### essentially all nodes are selected as the core.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ core-periphery}
