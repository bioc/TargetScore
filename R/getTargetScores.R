getTargetScores <- function(mirID, logFC, ...) {
  
	message(sprintf("Processing miRNA %s", mirID))
		
	targetScanCS <- get_TargetScanHuman_contextScore()

	targetScanPCT <- get_TargetScanHuman_PCT()
	
	mirtarbase <- get_validated_targets()
	
	mirfamily <- get_miRNA_family_info()
    
	if(missing(logFC)) {
		
		load(sprintf("%s/logFC.imputed.RData",system.file("extdata", package = "TargetScoreData")))
		
		logFC <- logFC.imputed[,colnames(logFC.imputed) == mirID, drop=F]
	}
	
	if(ncol(logFC)==1) {
		
		myTargetScores <- data.frame(logFC=as.numeric(as.matrix(logFC)))
		
	} else {
		
		myTargetScores <- as.data.frame(logFC)
	}			 			
		
	external_gene_id <- rownames(logFC)
					
	################## Validated targets ##################
	if(mirID %in% mirtarbase$miRNA) {		
		
		validated_targets <- unique(mirtarbase[mirtarbase$miRNA==mirID &
        mirtarbase$`Species (miRNA)` == "Homo sapiens" &
        mirtarbase$`Species (Target Gene)` == "Homo sapiens", "Target Gene"])    
					
		## get true labels from validated targets
		validated <- toupper(external_gene_id) %in% toupper(validated_targets)
	}
	
	################## Select fold-change ##################
	# for > 1 studies, pick the one that most correlates with the validation
	if(ncol(myTargetScores) > 1) {
		if(exists("validated")) {
			
			mycor <- sapply(1:ncol(myTargetScores), function(j) -cor(myTargetScores[,j], validated))			
			
			myTargetScores <- myTargetScores[,which.max(mycor),drop=F]			
			
		} else {			
			myTargetScores <- data.frame(logFC=apply(myTargetScores,1,mean))
		}
	}
  
	names(myTargetScores) <- "logFC"
			
	myTargetScores$external_gene_id <- external_gene_id
	
	if(exists("validated")) myTargetScores$validated <- validated
		
	#################### TargetScan context score ####################
	if(nrow(targetScanCS[targetScanCS$miRNA == mirID,]) != 0)  {
		targetScanCS.sel <- unique(targetScanCS[targetScanCS$miRNA == mirID,][, -3]) # remove tx column
			
		# take the rows with the most negative context+ score for genes having > 1 tx
		targetScanCS.sel$`context+ score` <- as.numeric(as.matrix(targetScanCS.sel$`context+ score`))
		
		targetScanCS.sel$`context+ score`[is.na(targetScanCS.sel$`context+ score`)] <- 0
			
		targetScanCS.sel <- aggregate(`context+ score`~`Gene Symbol`, targetScanCS.sel, min)
		
		myTargetScores$targetScanCS <- targetScanCS.sel[match(myTargetScores$external_gene_id, targetScanCS.sel$`Gene Symbol`), "context+ score"]
	}
	
			
	#################### TargetScan PCT ####################
	mirfamilyID <- unique(mirfamily[mirfamily$`MiRBase ID`==mirID, "miR family"])
	
	if(nrow(targetScanPCT[targetScanPCT$`miR Family` %in% mirfamilyID,]) != 0) {
	
		targetScanPCT.sel <- unique(targetScanPCT[targetScanPCT$`miR Family` %in% mirfamilyID, c("Gene Symbol", "PCT")])
		
		targetScanPCT.sel$PCT <- as.numeric(as.matrix(targetScanPCT.sel$PCT))
		
		targetScanPCT.sel$PCT[is.na(targetScanPCT.sel$PCT)] <- 0
		
		# take the rows with the most negative context+ score for genes having > 1 tx
		targetScanPCT.sel <- aggregate(PCT~`Gene Symbol`, targetScanPCT.sel, max)
		
		myTargetScores$targetScanPCT <- -targetScanPCT.sel[match(myTargetScores$external_gene_id, targetScanPCT.sel$`Gene Symbol`), "PCT"]
		
	}
			
	myTargetScores[is.na(myTargetScores)] <- 0
		
	#################### TargetScore ####################
    
	if(!is.null(myTargetScores$targetScanCS) & !is.null(myTargetScores$targetScanPCT)) {
								
		seqScores <- data.frame(targetScanCS=myTargetScores$targetScanCS, targetScanPCT=myTargetScores$targetScanPCT)
	}	
  
  if(!is.null(myTargetScores$targetScanCS) & is.null(myTargetScores$targetScanPCT)) {
  	
  		message("TargetScan PCT not found!")
  	
  		seqScores <- data.frame(targetScanCS=myTargetScores$targetScanCS)
  	}
  
  if(is.null(myTargetScores$targetScanCS) & !is.null(myTargetScores$targetScanPCT)) {
  		
  		message("TargetScan context score not found!")
  		
  		seqScores <- data.frame(targetScanPCT=myTargetScores$targetScanPCT)
  	}

  if(is.null(myTargetScores$targetScanCS) & is.null(myTargetScores$targetScanPCT)) {
  		message("TargetScan context score and PCT not found!")
  		seqScores <- NA
  	}	
	
	if(length(seqScores) > 1) {
		
		myTargetScores$targetScore <- targetScore(myTargetScores$logFC, seqScores, ...)
		
	} else {
		
		myTargetScores$targetScore <- targetScore(myTargetScores$logFC, ...)
	}
	
	myTargetScores
}

