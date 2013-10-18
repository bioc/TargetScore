\name{TargetScore-package}
\alias{TargetScore-package}
\alias{TargetScore}
\docType{package}
\title{
TargetScore: Infer microRNA targets using microRNA-overexpression data and sequence information
}
\description{
Infer the posterior distributions of microRNA targets by probabilistically modeling the likelihood microRNA-overexpression fold-changes and sequence-based scores. Variational Bayesian Gaussian mixture model (VB-GMM) is applied to log fold-changes and sequence scores to obtain the posteriors of latent variable being the miRNA targets. The final targetScore is computed as the sigmoid-transformed fold-change weighted by the averaged posteriors of target components over all of the features.
}
\details{
\tabular{ll}{
Package: \tab TargetScore\cr
Type: \tab Package\cr
Version: \tab 1.1.5\cr
Date: \tab 2013-10-15\cr
License: \tab GPL-2\cr
}
The front-end main function \code{\link{targetScore}} should be used to obtain the probablistic score of miRNA target. The workhourse function is \code{\link{vbgmm}}, which implementates multivariate variational Bayesian Gaussian mixture model.
}
\author{
Yue Li <yueli@cs.toronto.edu>
}
\references{
Lim, L. P., Lau, N. C., Garrett-Engele, P., Grimson, A., Schelter, J. M., Castle, J., Bartel,
D. P., Linsley, P. S., and Johnson, J. M. (2005). Microarray analysis shows that some microRNAs downregulate large numbers of target mRNAs. Nature, 433(7027), 769-773.

Bartel, D. P. (2009). MicroRNAs: Target Recognition and Regulatory Functions. Cell, 136(2), 215-233.

Bishop, C. M. (2006). Pattern recognition and machine learning. Springer, Information Science and Statistics. NY, USA. (p474-486)
}
\keyword{
  package
  microRNA target prediction
}
\seealso{
\code{\link{targetScore}}
}
\examples{
library(TargetScore)
ls("package:TargetScore")
}
