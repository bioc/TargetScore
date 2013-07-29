\name{initialization}
\alias{initialization}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Initialization of latent variable assignments (responsibility) of the VB-GMM (Internal function)
%%  ~~function to do ... ~~
}
\description{
Initialize latent varaibles based on the number of components. The function is run before the VB-EM iteration in vbgmm.
}
\usage{
initialization(X, init)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
D x N numeric vector or matrix of observations
}
  \item{init}{
%%     ~~Describe \code{init} here~~
Based on the dimension, init is expected to be one of the followings: scalar: number of components; vector: intial class labels; matrix: initialize with a D x K matrix for D variables and K components.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
The function is expected to be used by vbgmm to initialize assignments of latent varaibles before VM-EM iterations.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{R}{N by K matrix for N observations and K latent components (defined by init)}
}
\references{
Mo Chen (2012). Matlab code for Variational Bayesian Inference for Gaussian Mixture Model. http://www.mathworks.com/matlabcentral/fileexchange/35362-variational-bayesian-inference-for-gaussian-mixture-model
}
\author{
Yue Li
}
\seealso{
\code{\link{vbgmm}}
}
\examples{tmp <- initialization(matrix(c(rnorm(100,mean=2), rnorm(100,mean=3)),nrow=1), init=2)}