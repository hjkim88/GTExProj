###
#   File name : MakeTCGARcnt.R
#   Author    : Hyunjin Kim
#   Date      : Jul 9, 2018
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Make appropriate raw count RDA file with TCGA dataset
#
#   Instruction
#               1. Source("MakeTCGARcnt.R")
#               2. Run the function "makeTCGArCnt()" - specify the inputs and output path
#               3. A RDA object of TCGA read counts will be generated in the output path
#
#   Example
#               > source("The_directory_of_MakeTCGARcnt.R/MakeTCGARcnt.R")
#               > makeTCGArCnt(rCntPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/TCGA/ExpressionMatrices/tcga_read_counts/cancer-rawcounts.rda",
#                              marianoNormDatPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/TCGA/ExpressionMatrices/normalized_by_mariano/",
#                              outputPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/TCGA/ExpressionMatrices/tcga_read_counts/TCGA_28_countmat.rda")
###

makeTCGArCnt <- function(rCntPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/TCGA/ExpressionMatrices/tcga_read_counts/cancer-rawcounts.rda",
                         marianoNormDatPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/TCGA/ExpressionMatrices/normalized_by_mariano/",
                         outputPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/TCGA/ExpressionMatrices/tcga_read_counts/TCGA_28_countmat.rda") {
  
  ### load datasets
  load(rCntPath)
  
  ### get file names from Mariano's normalized data
  f <- list.files(marianoNormDatPath)
  f <- f[which(endsWith(f, ".rda"))]
  f <- f[-which(startsWith(f, "cancer"))]
  f <- f[-which(startsWith(f, "meni"))]
  f <- f[-which(startsWith(f, "net"))]
  
  ### raw count variable names
  tcga_raw_counts_vars <- NULL
  
  ### TRUE if a vector is all NA
  isAllNA <- function(v) {
    if(length(which(is.na(v))) == length(v)) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
  
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
    
    ### only with a tissue that has more than 100 samples
    if(ncol(cnt) >= 100) {
      ### get indicies with all NA
      rIdx <- which(apply(cnt, 1, isAllNA))
      
      ### remove all NA rows
      if(length(rIdx) > 0) {
        cnt <- cnt[-rIdx,]
      }
      
      ### order based on Entrez ID
      cnt <- cnt[order(as.numeric(rownames(cnt))),]
      
      ### data frame
      cnt <- data.frame(cnt, check.names = FALSE)
      
      ### save it with appropriate name
      assign(paste0("cmat_tcga_", strsplit(f[i], split = "-", fixed = TRUE)[[1]][1]), cnt, envir = globalenv())
      
      ### add the var name
      tcga_raw_counts_vars <- c(tcga_raw_counts_vars, paste0("cmat_tcga_", strsplit(f[i], split = "-", fixed = TRUE)[[1]][1]))
    }
    
    ### progress print
    writeLines(paste(i, "/", length(f)))
  }
  
  ### README function
  README <- function() {
    writeLines("The objects are raw counts of 28 TCGA tissues")
    writeLines("The raw counts were from Califano Lab")
    writeLines("/ifs/archive/shares/af_lab/TCGA/ExpressionMatrices/tcga_read_counts/cancer-rawcounts.rda")
    writeLines("Samples that are only appeared in Mariano's normalized-count-data were used")
    writeLines("Rows with all NAs were removed")
    writeLines("tcga_raw_counts_vars has all the name of 28 TCGA raw count object names")
    writeLines("For each tissue, the object name is cmat_tcga_<tumor_caronym>")
  }
  
  ### save all the variables in one RDA file
  save(list = c("tcga_raw_counts_vars", tcga_raw_counts_vars, "README"), file = outputPath)
  
}



