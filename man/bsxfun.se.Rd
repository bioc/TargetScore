\name{bsxfun.se}
\alias{bsxfun.se}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
bsxfun with single expansion (real Matlab style) (Internal function)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Depending on the dimension of x, repeat y in either by row or by column and apply element-wise operation defined by func.
}
\usage{
bsxfun.se(func, x, y, expandByRow = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{func}{
  function with two or more input parameters.
}
  \item{x, y}{
  two vectors, matrices, or arrays
}
  \item{expandByRow}{
  expand by row or by column of x when nrow(x)==ncol(x)==length(y)
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
The function is used by vbgmm.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{func(x,y)}{A matrix of having the same dimension of x.}
}

\author{
%%  ~~who you are~~
Yue Li
}
\note{
%%  ~~further notes~~
Internal function.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[pracma]{bsxfun}}
}
\examples{
bsxfun.se("*", matrix(c(1:10), nrow=2), matrix(c(1:5), nrow=5))
}