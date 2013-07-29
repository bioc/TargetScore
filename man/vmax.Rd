\name{vmax}
\alias{vmax}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Variational-Maximimization in VB-EM (Internal function)
}
\description{
The M step in VB-EM iteration.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
vmax(X, model, prior)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
D x N numeric vector or matrix of N observations (columns) and D variables (rows)
}
  \item{model}{
List containing model parameters (see \code{\link{vbgmm}})
}
  \item{prior}{
List containing the hyperparameters defining the prior distributions
}
}

\value{
\item{model}{A list containing the updated model parameters including alpha (Dirichlet), m (Gaussian mean), kappa (Gaussian variance), v (Wishart degree of freedom), M (Wishart precision matrix).
}}

\references{
Mo Chen (2012). Matlab code for Variational Bayesian Inference for Gaussian Mixture Model. http://www.mathworks.com/matlabcentral/fileexchange/35362-variational-bayesian-inference-for-gaussian-mixture-model

Bishop, C. M. (2006). Pattern recognition and machine learning. Springer, Information Science and Statistics. NY, USA. (p474-486)
}

\author{
Yue Li
}

\note{
X is expected to be D x N for N observations (columns) and D variables (rows)
}

\seealso{
\code{\link{vbgmm}}
}
\examples{
X <- c(rnorm(100,mean=2), rnorm(100,mean=3))
tmp <- vbgmm(X, tol=1e-3)
names(tmp$full.model)
}