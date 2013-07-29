\name{getTargetScores}
\alias{getTargetScores}

\title{
Compute targetScore of an overexpressed human microRNA
}
\description{
Obtain for each gene the targetScore using using pre-computed (logFC) TargetScan context score and PCT as sequence score. TargetScanData package is needed.
}
\usage{
getTargetScores(mirID, logFC, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mirID}{
A character string of microRNA ID (e.g., hsa-miR-1)
}
  \item{logFC}{
N x D numeric vector or matrix of logFC with D replicates for N genes.
}
  \item{\dots}{
Paramters passed to \code{\link{vbgmm}}
}
}
\details{
This is a conveinient function for computing targetScore for a human miRNA using user-supplied or pre-computed logFC and (if available) two pre-computed sequence scores namely TargetScan context score and PCT (probibility of conserved targeting). The function also searches for any validated targets from the MirTarBase human validated target list. The function requires TargetScanData to be installed first.
}
\value{
\item{targetScores}{numeric matrix of probabilistic targetScores together with the input variable and a binary vector indicating whether each gene is a valdiated target (if available).}
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
\code{\link{targetScore}}
}
\examples{
if(interactive()) {
  
  library(TargetScoreData)
  library(Biobase)
  library(GEOquery)

  # compute targetScore from pre-computed logFC and sequence socres
  # for hsa-miR-1
  mir1.score <- getTargetScores("hsa-miR-1", tol=1e-3, maxiter=200)

  # download  fold-change data from GEO for hsa-miR-124 overexpression in HeLa
    
  gset <- getGEO("GSE2075", GSEMatrix =TRUE, AnnotGPL=TRUE)

  if (length(gset) > 1) idx <- grep("GPL1749", attr(gset, "names")) else idx <- 1

  gset <- gset[[idx]]

  sampleinfo <- as.character(pData(gset)$title)

  geneInfo <- fData(gset)

  # only 24h data are used (discard 12h data)
  logfc.mir124 <- as.matrix(exprs(gset)[, 
    grep("HeLa transfected with miR-1 versus control transfected HeLa, 24 hours", sampleinfo)])
  
  rownames(logfc.mir124) <- geneInfo$`Gene symbol`
  
  mir124.score <- getTargetScores("hsa-miR-124", logfc.mir124, tol=1e-3, maxiter=200)
  
  head(mir124.score)
}
}




