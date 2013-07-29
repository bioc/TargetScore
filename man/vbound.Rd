\name{vbound}
\alias{vbound}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Variational Lower Bound Evaluation
}
\description{
Evaluate variational lower bound to determine when to stop VB-EM iteration (convergence).
}
\usage{
vbound(X, model, prior)
}

\arguments{
  \item{X}{
D x N numeric vector or matrix of N observations (columns) and D variables (rows)
}
  \item{model}{
List containing model parameters (see \code{\link{vbgmm}})
}
  \item{prior}{
numeric vector or matrix containing the hyperparameters for the prior distributions
}
}

\value{
A continuous scalar indicating the lower bound (the higher the more converged)
}
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
head(tmp$L) # lower bound should be strictly increasing
}
