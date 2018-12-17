###
#   File name : RefineCalifanoTCGA.R
#   Author    : Hyunjin Kim
#   Date      : Jun 27, 2018
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Normalize TCGA data of Califano lab in our way 
#
#   Instruction
#               1. Source("RefineCalifanoTCGA.R")
#               2. Run the function "refineTheirTCGA()" - specify the inputs and output path
#               3. A RDA object of normalized TCGA read counts will be generated in the output path
#
#   Example
#               > source("The_directory_of_RefineCalifanoTCGA.R/RefineCalifanoTCGA.R")
#               > refineTheirTCGA(rCntPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/TCGA/ExpressionMatrices/tcga_read_counts/cancer-rawcounts.rda",
#                                 marianoNormDatPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/TCGA/ExpressionMatrices/normalized_by_mariano/",
#                                 outputPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/TCGA/ExpressionMatrices/tcga_read_counts/vst_normalized_counts.rda")
###

refineTheirTCGA <- function(rCntPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/TCGA/ExpressionMatrices/tcga_read_counts/cancer-rawcounts.rda",
                            marianoNormDatPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/TCGA/ExpressionMatrices/normalized_by_mariano/",
                            outputPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/TCGA/ExpressionMatrices/tcga_read_counts/vst_normalized_counts.rda") {
  
  ### load datasets
  load(rCntPath)
  
  ### get file names from Mariano's normalized data
  f <- list.files(marianoNormDatPath)
  f <- f[which(endsWith(f, ".rda"))]
  f <- f[-which(startsWith(f, "cancer"))]
  
  ### a function returns logical value whether a given row has 0 or 1 across all samples
  isRubbish <- function(geneRow) {
    isZeroOrOne <- TRUE
    for(i in 1:length(geneRow)) {
      if(geneRow[i] > 1) {
        isZeroOrOne <- FALSE
        break
      }
    }
    
    return(isZeroOrOne)
  }
  
  ### A function to transform RNA-Seq data with VST in DESeq2 package
  RNASEQwithVST <- function(readCount) {
    
    ### load library
    if(!require(DESeq2)) {
      source("https://bioconductor.org/biocLite.R")
      biocLite("DESeq2")
      library(DESeq2)
    }
    
    ### make a design matrix for DESeq2 data
    condition <- data.frame(factor(rep("OneClass", ncol(readCount))))
    
    
    ### Data preparation for DESeq2 format
    deSeqData <- DESeqDataSetFromMatrix(countData=readCount, colData=condition, design= ~0)
    
    ### VST
    vsd <- vst(deSeqData)
    transCnt <- data.frame(assay(vsd))
    
    ### copy column names from original count data
    colnames(transCnt) <- colnames(readCount)
    
    return (transCnt)
  }
  
  ### normalized count variable names
  norm_counts_vars <- NULL
  
  ### iteratively perform normalization for each tissue
  for(i in 1:length(f)) {
    ### load Mariano's normalized file
    load(paste0(marianoNormDatPath, f[i]))
    
    ### get raw counts of samples which are also in Mariano's normalized file
    if(length(which(colnames(expmat) == "")) > 0) {
      cnt <- rawcounts[,colnames(expmat)[-which(colnames(expmat) == "")]]
    } else {
      cnt <- rawcounts[,colnames(expmat)]
    }
    
    ### change NA values to 0
    cnt[which(is.na(cnt), arr.ind = TRUE)] <- 0
    
    ### get rubbish indices
    rubbishCnt <- 0
    rubbishIdx <- 0
    for(j in 1:nrow(cnt)) {
      if(isRubbish(cnt[j,])) {
        rubbishCnt <- rubbishCnt + 1
        rubbishIdx[rubbishCnt] <- j
      }
    }
    
    ### remove rubbish
    cnt <- cnt[-rubbishIdx,]
    
    ### order based on Entrez ID
    cnt <- cnt[order(as.numeric(rownames(cnt))),]
    
    ### normalize the counts with VST
    norm_cnt <- RNASEQwithVST(cnt)
    
    ### save it with appropriate name
    assign(paste0("emat_tcga2_", strsplit(f[i], split = "-", fixed = TRUE)[[1]][1]), norm_cnt, envir = globalenv())
    
    ### add the var name
    norm_counts_vars <- c(norm_counts_vars, paste0("emat_tcga2_", strsplit(f[i], split = "-", fixed = TRUE)[[1]][1]))
    
    ### progress print
    writeLines(paste(i, "/", length(f)))
  }
  
  ### README function
  README = function(){
    writeLines("The objects are normalized counts of 37 TCGA tissues")
    writeLines("The raw counts were from Califano Lab")
    writeLines("Samples that are only appeared in Mariano's normalized-count-data were used")
    writeLines("All the NAs in the raw count data were changed to 0")
    writeLines("Genes with 0 or 1 across all samples were removed")
    writeLines("VST-normalization of DESeq2 package was used")
    writeLines("norm_counts_vars has all the name of 37 TCGA normalized counts")
    writeLines("For each tissue, the object name is emat_tcga2_<tumor_caronym>")
  }
  
  ### save all the variables in one RDA file
  save(list = c("norm_counts_vars", norm_counts_vars, "README"), file = outputPath)
  
}



