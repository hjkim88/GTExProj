###
#   File name : ARACNeReady.R
#   Author    : Hyunjin Kim
#   Date      : Dec 1, 2017
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Make Aracne-ready file from the transformed counts
#
#   Instruction
#               1. Source("ARACNeReady.R")
#               2. Run the function "makeReady()" - specify the input file (transformed counts) directory and output directory
#               3. ARACNe-ready data will be generated in the output directory
#
#   Example
#               > source("The_directory_of_ARACNeReady.R/ARACNeReady.R")
#               > makeReady(inputPath="./results/transformed_counts/", outputPath="./results/aracne_ready/")
###


makeReady <- function(inputPath="./results/transformed_counts/", outputPath="./results/aracne_ready/") {
  
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
    cnt <- cbind(Gene=rownames(cnt), cnt[,-1])
    
    ### save the transformed data
    write.table(cnt, paste0(outputPath, refineFileName(substr(f[i], 1, nchar(f[i])-4)), ".dat"), sep = "\t", row.names = FALSE, quote = FALSE)
  }
  
}


