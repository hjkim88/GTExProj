###
#   File name : VST_Counts.R
#   Author    : Hyunjin Kim
#   Date      : Nov 7, 2017
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Perform VST(Variance Stabilizing Transformation) on counts 
#
#   Instruction
#               1. Source("VST_Counts.R")
#               2. Run the function "vst_counts()" - specify the input file (filtered counts) directory and output directory
#               3. VS-Transformed counts will be generated in the output directory
#
#   Example
#               > source("The_directory_of_VST_Counts.R/VST_Counts.R")
#               > vst_counts(cntPath="./results/separated_counts/",
#                            outputPath="./results/transformed_counts/")
###

vst_counts <- function(cntPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/ExpressionMatrices/GTEx_processed/separated_counts/",
                       outputPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/ExpressionMatrices/GTEx_processed/transformed_counts/") {
  
  ### load library
  if(!require(DESeq2)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("DESeq2")
    library(DESeq2)
  }
  
  
  ### A function to transform RNA-Seq data with VST in DESeq2 package
  ### input readCount of this function includes gene symbols in the first column
  RNASEQwithVST <- function(readCount) {
    
    ### make a design matrix for DESeq2 data
    #Coldata <- data.frame(sampleType)
    #rownames(Coldata) <- colnames(readCount)
    condition <- data.frame(factor(rep("OneClass", ncol(readCount)-1)))
    
    
    ### Data preparation for DESeq2 format
    deSeqData <- DESeqDataSetFromMatrix(countData=readCount[,-1], colData=condition, design= ~0)
    ### Remove rubbish rows - this will decrease the number of rows
    #deSeqData <- deSeqData[rowSums(counts(deSeqData))>1,]
    
    ### VST
    vsd <- vst(deSeqData)
    transCnt <- data.frame(assay(vsd))
  
    ### Add gene symbol column and copy column names from original count data
    transCnt <- cbind(readCount$Gene_Symbol, transCnt)
    colnames(transCnt) <- colnames(readCount)
    
    return (transCnt)
  }
  
  
  ### collect files from the cntPath
  f <- list.files(cntPath)
  
  
  ### iteratively perform VST
  for(i in 1:length(f)) {
    
    ### load filtered counts
    cnt <- read.table(paste0(cntPath, f[i]), sep="\t", row.names = 1, header = TRUE, check.names = FALSE)
    
    ### perform VST
    tCnt <- RNASEQwithVST(cnt)
    
    ### save the transformed data
    tCnt <- cbind(Entrez_ID=rownames(tCnt), tCnt)
    tCnt <- tCnt[order(as.numeric(as.character(tCnt$Entrez_ID))),]
    write.table(tCnt, paste0(outputPath, substr(f[i], 1, nchar(f[i])-4), "_vst.txt"), sep = "\t", row.names = FALSE)
  }
  
}

