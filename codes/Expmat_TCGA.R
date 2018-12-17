###
#   File name : Expmat_TCGA.R
#   Author    : Hyunjin Kim
#   Date      : Apr 4, 2018
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : load expressions and make emat files for TCGA datasets
#
#   Instruction
#               1. Source("Expmat_TCGA.R")
#               2. Run the function "makeReady()" - specify the input file (TCGA RDAs) directory and output path
#               3. TCGA expression matrices (and their combined RDA file) will be generated
#
#   Example
#               > source("The_directory_of_Expmat_TCGA.R/Expmat_TCGA.R")
#               > ematReady(inputPath="./data/TCGA/ExpressionMatrices/")
###

ematReady <- function(inputPath="./data/TCGA/ExpressionMatrices/", outputPath="./results/transformed_counts/TCGA/tcga_28_emat.rda") {
  
  ### collect files from the cntPath
  f <- list.files(inputPath)
  f <- f[intersect(which(endsWith(f, ".rda")), which(f != "cancer-expmat.rda"))]
  
  emat_tcga_names <- NULL
  
  ### iteratively perform transforming
  for(i in 1:length(f)) {
    
    load(paste0(inputPath, f[i]))
    
    if(length(which(samples[,2] == "normal")) + length(which(samples[,2] == "primary.tumor")) + length(which(samples[,2] == "recurrent.tumor")) + length(which(samples[,2] == "met.tumor")) == nrow(samples)) {
      rIdx <- union(which(samples[,2] == "normal"), which(!startsWith(samples[,1], "TCGA")))
      
      if(length(rIdx) > 0) {
        expmat <- expmat[,-rIdx]
        samples <- samples[-rIdx,]
      }
      
      if(nrow(samples) >= 100) {
        expmat <- data.frame(expmat, check.names = FALSE)
        expmat <- cbind(Gene=rownames(expmat), expmat)
        expmat <- expmat[order(as.numeric(as.character(expmat$Gene))),]
        expmat <- expmat[,-1]
        
        emat_tcga_names <- c(emat_tcga_names, paste0("emat_tcga_", strsplit(f[i], "-", fixed = TRUE)[[1]][1]))
        assign(paste0("emat_tcga_", strsplit(f[i], "-", fixed = TRUE)[[1]][1]), expmat, envir = globalenv())
      } else {
        writeLines(paste0(f[i], " has samples less than 100"))
      }
      
    } else {
      writeLines("sample numbers do not match")
    }
    
  }
  
  assign("emat_tcga_names", emat_tcga_names, envir = globalenv())
  
  vars <- c(emat_tcga_names, "emat_tcga_names")
  save(list = vars, file = outputPath)
  
  rm(expmat)
  rm(samples)
  
}