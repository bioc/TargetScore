# TargetScore
# logFC: N x 1 data matrix for log fold-change of N genes
# seqscores: N x D sequence-based scores matrix for D sequence-based features
# Author: Yue Li (yueli@cs.toronto.edu)

targetScore <- function(logFC, seqScores, ...) {
    
  if(!missing(seqScores)) { 
    
    stopifnot(nrow(logFC) == nrow(seqScores))        
  }
  
  message(sprintf("\n\nCalculate targetScores for %s genes\n", length(logFC)))
  
  stage <- 1
  
  message(sprintf("%s. Infer target component from logFC:\n", as.roman(stage)))
    
  posteriors <- vbgmm(logFC, init=3, ...)$R[,1,drop=F]
  
  if(!missing(seqScores)) {
        
    p.seq <- do.call("cbind", lapply(1:ncol(seqScores),
                                     
            function(j) {
            
            message(sprintf("%s. Infer target component from %s:\n",
                    as.roman(stage+j), colnames(seqScores)[j]))
            
            vbgmm(seqScores[,j,drop=F], init=2, ...)$R[,1,drop=F]
    }))
    
    posteriors <- cbind(posteriors, p.seq)
  }      
  
  sigmoid(-logFC) * apply(posteriors, 1, mean)  
}
