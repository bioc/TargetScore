\name{vbgmm}
\alias{vbgmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Variational Bayesian Gaussian mixture model (VB-GMM)
%%  ~~function to do ... ~~
}
\description{
Given a N x D matrix of N observations and D variables, compute VB-GMM via VB-EM.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
vbgmm(data, init = 2, prior, tol = 1e-20, maxiter = 2000, mirprior = TRUE, expectedTargetFreq = 0.01, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
N x D numeric vector or matrix of N observations (rows) and D variables (columns)
}
  \item{init}{
Based on the dimension, init is expected to be one of the followings: scalar: number of components; vector: intial class labels; matrix: initialize with a D x K matrix for D variables and K components.
}
  \item{prior}{
A list containing the hyperparameters including alpha (Dirichlet), m (Gaussian mean), kappa (Gaussian variance), v (Wishart degree of freedom), M (Wishart precision matrix). 
}
  \item{tol}{
Threshold that defines termination/convergence of VB-EM when abs(L[t] - L[t-1])/abs(L[t]) < tol
}
  \item{maxiter}{
Scalar for maximum number of EM iterations
}
  \item{mirprior}{
Boolean to indicate whether to use expectedTargetFreq to initialize alpha0 for the hyperparameters of Dirichlet.
}
  \item{expectedTargetFreq}{
Expected target frequence within the gene population. By default, it is set to 0.01, which is consistent with the widely accepted prior knoweldge that 200/20000 targets per miRNA.
}
  \item{verbose}{
Boolean indicating whether to show progress in terms of lower bound (\code{\link{vbound}}) of VB-EM (default: FALSE)
}
}
\details{
The function implements variation Bayesian multivariate GMM described in Bishop (2006). Please refer to the reference below for more details. This is the workhorse of \code{\link{targetScore}}. Alternatively, user can choose to apply this function to other problems other than miRNA target prediction. 
}
\value{
A list containing:
\item{label}{a vector of maximum-a-posteriori (MAP) assignments of latent discrete values based on the posteriors of latent variables.}
\item{R}{N x D matrix of posteriors of latent variables}
\item{mu}{Gaussian means of the latent components}
\item{full.model}{A list containing posteriors R, logR, and the model parameters including alpha (Dirichlet), m (Gaussian mean), kappa (Gaussian variance), v (Wishart degree of freedom), M (Wishart precision matrix)}
\item{L}{A vector of variational lower bound at each EM iterations (should be strictly increasing)}
}


\references{
Mo Chen (2012). Matlab code for Variational Bayesian Inference for Gaussian Mixture Model. http://www.mathworks.com/matlabcentral/fileexchange/35362-variational-bayesian-inference-for-gaussian-mixture-model

Bishop, C. M. (2006). Pattern recognition and machine learning. Springer, Information Science and Statistics. NY, USA. (p474-486)
}

\author{
Yue Li
}

\seealso{
\code{\link{targetScore}}
}

\examples{
X <- c(rnorm(100,mean=2), rnorm(100,mean=3))
tmp <- vbgmm(X, tol=1e-3)
names(tmp)
}