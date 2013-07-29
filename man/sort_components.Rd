\name{sort_components}
\alias{sort_components}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sort mixture components in increasing order of averaged means (Internal function)
%%  ~~function to do ... ~~
}
\description{
Sort Gaussian mixture components with model paramters in increasing order of averaged means of d variables.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sort_components(model)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
%%     ~~Describe \code{model} here~~
A list containing trained parameters of the Baysian GMM (see Value section in \code{\link{vbgmm}}).
}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
VB-GMM model list in increasing order of averaged means.
}
\author{
Yue Li
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{vbgmm}}
}
\examples{
tmp <- vbgmm(c(rnorm(100,mean=2), rnorm(100,mean=3)), tol=1e-3)
tmp$mu
}