\name{targetScore}
\alias{targetScore}

\title{
Probabilistic score of genes being the targets of an overexpressed microRNA
}
\description{
Given the overexpression fold-change and sequence-scores (optional) of all of the genes, calculate for each gene the TargetScore as a probability of miRNA target.
}
\usage{
targetScore(logFC, seqScores, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{logFC}{
numeric vector of log fold-changes of N genes in treatment (miRNA overexpression) vs control (mock).
}
  \item{seqScores}{
N x D numeric vector or matrix of D sequence-scores for N genes. Each score vector is expected to be equal to or less than 0. The more negative the scores, the more likely the corresponding target.
}
  \item{\dots}{
Paramters passed to \code{\link{vbgmm}}
}
}
\details{
Given expression fold-change (due to miRNA transfection), we use a three-component VB-GMM to infer down-regulated targets accounting for genes with little or positive log fold-change (due to off-target effects (Khan et al., 2009). Otherwise, two-component VB-GMM is applied to unsigned sequence scores (seqScores). The parameters for the VB-GMM are optimized using Variational Bayesian Expectation-Maximization (VB-EM) algorithm. Presumably, the mixture component with the largest absolute means of observed negative fold-change or sequence score is associated with miRNA targets and denoted as "target component". The other components  correspond to the "background component". It follows that inferring  miRNA-mRNA interactions most likely explained by the observed data is equivalent to inferring the posterior distribution of the target component given the observed variables. The targetScore is computed as the sigmoid-transformed fold-change weighted by the averaged posteriors of target components over all of the features. Specifically, we define the targetScore as a composite probabilistic score of a gene being the target t of a miRNA:

sigmoid(-log{FC}) (1/K+1) sum_{x in \{x_f, x_1, ..., x_L\} p(t | x)),

where sigmoid(-log{FC}) = 1/(1 + exp(log{FC})) and p(t | x) is the posterior of the first component computed by \code{\link{vbgmm}}.
}}
\value{
\item{targetScore}{numeric vector of probabilistic targetScores for N genes}
}
\references{
%% ~put references to the literature/web site here ~
Lim, L. P., Lau, N. C., Garrett-Engele, P., Grimson, A., Schelter, J. M., Castle, J., Bartel,
D. P., Linsley, P. S., and Johnson, J. M. (2005). Microarray analysis shows that some microRNAs downregulate large numbers of target mRNAs. Nature, 433(7027), 769-773.

Bartel, D. P. (2009). MicroRNAs: Target Recognition and Regulatory Functions. Cell, 136(2), 215-233.

Bishop, C. M. (2006). Pattern recognition and machine learning. Springer, Information Science and Statistics. NY, USA. (p474-486)
}
\author{
Yue Li
}

\seealso{
\code{\link{vbgmm}}
}
\examples{
# A toy example:
# 10 down-reg, 1000 unchanged, 90 up-reg genes 
# due to overexpressing a miRNA
trmt <- c(rnorm(10,mean=0.01), rnorm(1000,mean=1), rnorm(90,mean=2)) + 1e3

ctrl <- c(rnorm(1100,mean=1)) + 1e3

logFC <- log2(trmt) - log2(ctrl)

# 8 out of the 10 down-reg genes have prominent seq score A
seqScoreA <- c(rnorm(8,mean=-2), rnorm(1092,mean=0))

# 10 down-reg genes plus 10 more genes have prominent seq score B
seqScoreB <- c(rnorm(20,mean=-2), rnorm(1080,mean=0))

seqScores <- cbind(seqScoreA, seqScoreB)              

p.targetScore <- targetScore(logFC, seqScores, tol=1e-3)
}