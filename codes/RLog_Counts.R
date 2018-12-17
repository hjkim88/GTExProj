###
#   File name : RLog_Counts.R
#   Author    : Hyunjin Kim
#   Date      : Nov 7, 2017
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Perform rlog() function on counts 
#
#   Instruction
#               1. Source("RLog_Counts.R")
#               2. Run the function "rlog_counts()" - specify the input file (filtered counts) directory and output directory
#               3. The transformed counts will be generated in the output directory
#
#   Example
#               > source("The_directory_of_RLog_Counts.R/RLog_Counts.R")
#               > rlog_counts(cntPath="./results/separated_counts/",
#                             outputPath="./results/transformed_counts/")
###

rlog_counts <- function(cntPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/ExpressionMatrices/GTEx_processed/separated_counts/",
                        outputPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/ExpressionMatrices/GTEx_processed/transformed_counts/") {
  
  ### load library
  if(!require(BiocParallel)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("BiocParallel")
    library(BiocParallel)
  }
  
  
  ### A function to transform RNA-Seq data with rlog() in DESeq2 package
  ### input readCount of this function includes gene symbols in the first column
  RNASEQwithRLog <- function(fileName) {
    
    ### load library
    if(!require(DESeq2)) {
      source("https://bioconductor.org/biocLite.R")
      biocLite("DESeq2")
      library(DESeq2)
    }
    
    ### load filtered counts
    readCount <- read.table(paste0(cntPath, fileName), sep="\t", row.names = 1, header = TRUE, check.names = FALSE)
    
    
    ### make a design matrix for DESeq2 data
    #Coldata <- data.frame(sampleType)
    #rownames(Coldata) <- colnames(readCount)
    condition <- data.frame(factor(rep("OneClass", ncol(readCount)-1)))
    
    
    ### Data preparation for DESeq2 format
    deSeqData <- DESeqDataSetFromMatrix(countData=readCount[,-1], colData=condition, design= ~0)
    ### Remove rubbish rows - this will decrease the number of rows
    #deSeqData <- deSeqData[rowSums(counts(deSeqData))>1,]
    
    ### rlog
    rld <- rlog(deSeqData)
    transCnt <- data.frame(assay(rld))
    
    ### Add gene symbol column and copy column names from original count data
    transCnt <- cbind(readCount$Gene_Symbol, transCnt)
    colnames(transCnt) <- colnames(readCount)
    
    ### save the transformed data
    transCnt <- cbind(Entrez_ID=rownames(transCnt), transCnt)
    write.table(transCnt, paste0(outputPath, substr(fileName, 1, nchar(fileName)-4), "_rlog.txt"), sep = "\t", row.names = FALSE)
    
    ### print current time
    print(Sys.time())
  }
  
  
  ### collect files from the cntPath
  f <- list.files(cntPath)
  
  
  ### perform rlog() parallely
  bplapply(f[1:length(f)], RNASEQwithRLog)
  
}

