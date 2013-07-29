\name{logmvgamma}
\alias{logmvgamma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Logarithmic multivariate Gamma function (Internal function)
%%  ~~function to do ... ~~
}
\description{
Compute logarithm multivariate Gamma function.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
logmvgamma(x, d)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
numeric vector or matrix
}
  \item{d}{
dimension
}
}
\details{
Gamma_p(x) = pi^(p(p-1)/4) prod_(j=1)^p Gamma(x+(1-j)/2)

log Gamma_p(x) = p(p-1)/4 log pi + sum_(j=1)^p log Gamma(x+(1-j)/2)
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Matrix of the same dimension as x.
}
\references{
Mo Chen (2012). Matlab code for Variational Bayesian Inference for Gaussian Mixture Model. http://www.mathworks.com/matlabcentral/fileexchange/35362-variational-bayesian-inference-for-gaussian-mixture-model
}
\author{
Yue Li
}
\seealso{
\code{\link[base]{lgamma}}
}
\examples{logmvgamma(matrix(1:6,nrow=3), 2)}