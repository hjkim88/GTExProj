###
#   File name : ARACNeReady_GTEx.R
#   Author    : Hyunjin Kim
#   Date      : Mar 1, 2018
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Make Aracne-ready file from the transformed counts from GTEx
#               Does not create files for tissues that have less than 100 samples
#               If the number of sample exceeds 200, 200 samples will be selected based on RIN
#
#   Instruction
#               1. Source("ARACNeReady_GTEx.R")
#               2. Run the function "makeReady()" - specify the input file (transformed counts) directory and output directory
#               3. ARACNe-ready data will be generated in the output directory
#
#   Example
#               > source("The_directory_of_ARACNeReady_GTEx.R/ARACNeReady_GTEx.R")
#               > makeReady(inputPath="./results/transformed_counts/vst_clean/",
#                           sampleInfoPath="./data/GTEx_Data_V6_SampleData.csv",
#                           outputPath="./results/aracne_ready/GTEx/vst_clean/")
###

makeReady <- function(inputPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/ExpressionMatrices/GTEx_processed/transformed_counts/vst_clean/",
                      sampleInfoPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/Annotations/GTEx_Data_V6_SampleData.csv",
                      outputPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/ExpressionMatrices/GTEx_processed/aracne_ready_counts/vst_clean/") {
  
  ### load sample info
  sampleInfo <- read.csv(sampleInfoPath, check.names = FALSE, stringsAsFactors = FALSE)
  
  
  ### make RIN info
  rin <- sampleInfo$SMRIN
  names(rin) <- sampleInfo$SAMPID
  
  
  ### collect files from the cntPath
  f <- list.files(inputPath)
  
  
  ### A function to change "(" or ")" to "-"
  refineFileName <- function(str_line) {
    result_line <- gsub("\\(", "-", str_line)
    result_line <- gsub("\\)", "-", result_line)
    
    return(result_line)
  }
  
  
  ### iteratively perform transforming
  for(i in 1:length(f)) {
    ### load filtered counts
    cnt <- read.table(paste0(inputPath, f[i]), sep="\t", row.names = 1, header = TRUE, check.names = FALSE)
    
    if((ncol(cnt)-1) > 200) {
      ### get RIN for the samples
      temp <- rin[colnames(cnt)[-1]]
      temp <- temp[order(-temp)]
      
      ### only keep top 200 samples based on RIN
      cnt <- cbind(Gene=rownames(cnt), cnt[,names(temp)[1:200]])
      
      ### save the transformed data
      write.table(cnt, paste0(outputPath, refineFileName(substr(f[i], 1, nchar(f[i])-4)), ".dat"), sep = "\t", row.names = FALSE, quote = FALSE)
    
    } else if((ncol(cnt)-1) >= 100) {
      cnt <- cbind(Gene=rownames(cnt), cnt[,-1])
      
      ### save the transformed data
      write.table(cnt, paste0(outputPath, refineFileName(substr(f[i], 1, nchar(f[i])-4)), ".dat"), sep = "\t", row.names = FALSE, quote = FALSE)
    }
  }
}






