###
#   File name : CleanRows.R
#   Author    : Hyunjin Kim
#   Date      : Feb 28, 2018
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Remove genes that have 0 or 1 across all samples
#
#   Instruction
#               1. Source("CleanRows.R")
#               2. Run the function "cleanRows()" - specify the input file (filtered counts) directory and output directory
#               3. Cleaned counts will be generated in the output directory
#
#   Example
#               > source("The_directory_of_CleanRows.R/CleanRows.R")
#               > cleanRows(cntPath="./results/separated_counts/",
#                           outputPath="./results/cleaned_counts/")
###


cleanRows <- function(cntPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/ExpressionMatrices/GTEx_processed/separated_counts/",
                      outputPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/ExpressionMatrices/GTEx_processed/cleaned_counts/") {
  
  ### collect files from the cntPath
  f <- list.files(cntPath)
  
  
  ### a function returns logical value whether a given row has NA, 0 or 1 across all samples
  isRubbish <- function(geneRow) {
    isZeroOrOne <- TRUE
    for(i in 1:length(geneRow)) {
      if((!is.na(geneRow[i])) && (geneRow[i] > 1)) {
        isZeroOrOne <- FALSE
        break
      }
    }
    
    return(isZeroOrOne)
  }
  
  
  ### iteratively clean
  for(i in 1:length(f)) {
    
    ### load filtered counts
    cnt <- read.table(paste0(cntPath, f[i]), sep="\t", row.names = 1, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
    
    ### get rubbish indices
    rubbishCnt <- 0
    rubbishIdx <- 0
    for(j in 1:nrow(cnt)) {
      if(isRubbish(cnt[j,2:ncol(cnt)])) {
        rubbishCnt <- rubbishCnt + 1
        rubbishIdx[rubbishCnt] <- j
      }
    }
    
    ### remove rubbish
    cnt <- cnt[-rubbishIdx,]
    
    ### save cleaned counts
    cnt <- cbind(Entrez_ID=rownames(cnt), cnt)
    cnt <- cnt[order(as.numeric(as.character(cnt$Entrez_ID))),]
    write.table(cnt, paste0(outputPath, substr(f[i], 1, nchar(f[i])-4), "_clean.txt"), sep = "\t", row.names = FALSE)
  }
  
}
