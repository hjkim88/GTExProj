###
#   File name : SeparateFiles.R
#   Author    : Hyunjin Kim
#   Date      : Nov 7, 2017
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Some tissue samples have sub-regions. Separate them based on sub-regions.
#
#   Instruction
#               1. Source("SeparateFiles.R")
#               2. Run the function "separate()" - specify the input file (filtered counts) directory, sample information path, and output directory
#               3. Separated counts will be generated in the output directory
#
#   Example
#               > source("The_directory_of_SeparateFiles.R/SeparateFiles.R")
#               > separate(fCntPath="./results/aggregated_counts/",
#                          sampleInfoPath="./data/GTEx_Data_V6_SampleData.csv",
#                          outputPath="./results/separated_counts/")
###

separate <- function(fCntPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/ExpressionMatrices/GTEx_processed/separated_counts/",
                     sampleInfoPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/Annotations/GTEx_Data_V6_SampleData.csv",
                     outputPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/ExpressionMatrices/GTEx_processed/separated_counts/") {
  
  
  ### load & preprocess sample info
  sampleInfo <- read.csv(sampleInfoPath, row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
  sampleInfo <- sampleInfo[,c("SMTSD", "SMNABTCH", "SMNABTCHT", "SMNABTCHD", "SMGEBTCH", "SMGEBTCHD", "SMGEBTCHT")]
  
  
  ### collect files from the fCntPath
  f <- list.files(fCntPath)
  
  
  ### A function to remove spaces in a string
  removeSpace <- function(str) {
    new_str <- unlist(strsplit(str, " ", fixed = TRUE))
    result_str <- ""
    for(i in 1:length(new_str)) {
      result_str <- paste0(result_str, new_str[i])      
    }
    
    return(result_str)
  }
  
  
  ### iteratively perform separation
  for(i in 1:length(f)) {
    
    ### load filtered counts
    fCnt <- read.table(paste0(fCntPath, f[i]), sep="\t", row.names = 1, header=TRUE, check.names = FALSE)
    
    ### set sub-region info
    group <- as.data.frame(as.character(sampleInfo[colnames(fCnt)[2:ncol(fCnt)], "SMTSD"]))
    rownames(group) <- colnames(fCnt)[2:ncol(fCnt)]
    colnames(group) <- "SMTSD"
    
    ### count unique sub-regions
    uniqueSN <- as.character(unique(group$SMTSD))
    uniqueLen <- length(uniqueSN)
    
    ### if there are more than one sub-regions, separate the file
    if(uniqueLen > 1) {
      for(j in 1:uniqueLen) {
        tempCnt <- fCnt[,union(1, which(group$SMTSD == uniqueSN[j])+1)]
        tempCnt <- cbind(Entrez_ID=rownames(tempCnt), tempCnt)
        tempCnt <- tempCnt[order(as.numeric(as.character(tempCnt$Entrez_ID))),]
        write.table(tempCnt, paste0(outputPath, removeSpace(uniqueSN[j]), ".txt"), sep = "\t", row.names = FALSE)
      }
    }
    else {
      fCnt <- cbind(Entrez_ID=rownames(fCnt), fCnt)
      fCnt <- fCnt[order(as.numeric(as.character(fCnt$Entrez_ID))),]
      write.table(fCnt, paste0(outputPath, f[i]), sep = "\t", row.names = FALSE)
    }
    
  }
  
}
