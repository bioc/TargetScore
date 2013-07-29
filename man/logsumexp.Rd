\name{logsumexp}
\alias{logsumexp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute log(sum(exp(x),dim)) while avoiding numerical underflow (Internal function)
%%  ~~function to do ... ~~
}
\description{
Compute log(sum(exp(x),dim)) while avoiding numerical underflow.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
logsumexp(x, margin = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
numeric vector or matrix
}
  \item{margin}{
dimension to apply summation
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
numeric vector or matrix of the same columns or rows (depending on margin) as x
}
\references{
Mo Chen (2012). Matlab code for Variational Bayesian Inference for Gaussian Mixture Model. http://www.mathworks.com/matlabcentral/fileexchange/35362-variational-bayesian-inference-for-gaussian-mixture-model
}
\author{
Yue Li
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{logsumexp(matrix(c(1:5)), 2)}