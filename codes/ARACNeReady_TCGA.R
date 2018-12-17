###
#   File name : ARACNeReady_TCGA.R
#   Author    : Hyunjin Kim
#   Date      : Feb 1, 2018
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Make Aracne-ready file from the TCGA RDAs
#
#   Instruction
#               1. Source("ARACNeReady_TCGA.R")
#               2. Run the function "makeReady()" - specify the input file (TCGA RDAs) directory and output directory
#               3. ARACNe-ready data will be generated in the output directory
#
#   Example
#               > source("The_directory_of_ARACNeReady_TCGA.R/ARACNeReady_TCGA.R")
#               > makeReady(inputPath="./data/TCGA/ExpressionMatrices/", outputPath="./results/aracne_ready/TCGA/")
###

makeReady <- function(inputPath="./data/TCGA/ExpressionMatrices/", outputPath="./results/aracne_ready/TCGA/") {
  
  ### collect files from the cntPath
  f <- list.files(inputPath)
  
  
  ### iteratively perform transforming
  for(i in 1:length(f)) {
    if(endsWith(f[i], ".rda") && (f[i] != "cancer-expmat.rda")) {
      
      load(paste0(inputPath, f[i]))
      
      if(length(which(samples[,2] == "normal")) + length(which(samples[,2] == "primary.tumor")) + length(which(samples[,2] == "recurrent.tumor")) + length(which(samples[,2] == "met.tumor")) == nrow(samples)) {
        rIdx <- union(which(samples[,2] == "normal"), which(!startsWith(samples[,1], "TCGA")))
        
        if(length(rIdx) > 0) {
          expmat <- expmat[,-rIdx]
          samples <- samples[-rIdx,]
        }
        
        if(nrow(samples) >= 100 && nrow(samples) <= 200) {
          expmat <- data.frame(expmat, check.names = FALSE)
          expmat <- cbind(Gene=rownames(expmat), expmat)
          expmat <- expmat[order(as.numeric(as.character(expmat$Gene))),]
          
          ### save the transformed data
          write.table(expmat, paste0(outputPath, "tcga_", strsplit(f[i], "-", fixed = TRUE)[[1]][1], ".dat"), sep = "\t", row.names = FALSE, quote = FALSE)
        } else if(nrow(samples) > 200) {
          set.seed(1234)
          nIdx <- sample(nrow(samples), 200)
          
          expmat <- expmat[,nIdx]
          
          expmat <- data.frame(expmat, check.names = FALSE)
          expmat <- cbind(Gene=rownames(expmat), expmat)
          expmat <- expmat[order(as.numeric(as.character(expmat$Gene))),]
          
          ### save the transformed data
          write.table(expmat, paste0(outputPath, "tcga_", strsplit(f[i], "-", fixed = TRUE)[[1]][1], ".dat"), sep = "\t", row.names = FALSE, quote = FALSE)
        } else {
          writeLines(paste0(f[i], " has samples less than 100"))
        }
        
      } else {
        writeLines("sample numbers do not match")
      }
      
    }
  }
  
}