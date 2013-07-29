\name{vexp}
\alias{vexp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Variational-Expectation in VB-EM (Internal function)
%%  ~~function to do ... ~~
}
\description{
The E step in VB-EM iteration.
}
\usage{
vexp(X, model)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
D x N numeric vector or matrix of N observations (columns) and D variables (rows)
}
  \item{model}{
List containing model parameters (see \code{\link{vbgmm}})
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
dim(tmp$R); head(tmp$R)
}
