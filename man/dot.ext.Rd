\name{dot.ext}
\alias{dot.ext}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Elementwise dot product (modified dot function) (Internal function)
}
\description{
Same as dot but handle single row matrix differently by multiplying each value but not sum them up
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dot.ext(x, y, mydim)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
numeric vector or matrix
}
  \item{y}{
numeric vector or matrix
}
  \item{mydim}{
Elementwise product (if 1); otherwise defined by \code{\link[pracma]{dot}}
}
}
\details{
Returns the 'dot' or 'scalar' product of vectors or columns of matrices. Two vectors must be of same length, two matrices must be of the same size. If x and y are column or row vectors, their dot product will be computed IF mydim is 1 (only difference from \code{\link[pracma]{dot}}).
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A scalar or vector of length the number of columns of x and y.
}

\author{
Yue Li
}

\seealso{
\code{\link[pracma]{dot}}
}
\examples{
dot.ext(1:5, 1:5)
dot.ext(1:5, 1:5, 1)
}