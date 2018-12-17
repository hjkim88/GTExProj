###
#   File name : CheckDividable.R
#   Author    : Hyunjin Kim
#   Date      : Mar 22, 2018
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : To check every non-bootstrapped network has complete set of its TF-target pairs
#
#   Instruction
#               1. Source("CheckDividable.R")
#               2. Run the function "check()" - specify the inputs (network dir, exp dir, and hubs list)
#               3. The result info will be printed in console
#
#   Example
#               > source("The_directory_of_CheckDividable.R/CheckDividable.R")
#               > check(netPath="./results/Aracne/GTEx2/MI/", expPath="./results/aracne_ready/GTEx/vst_clean/", hubFilePath="../../ARACNe/hubGenes/hubs.txt")
###

check <- function(netPath="E:/CUMC/GTEx/MI/", expPath="./results/aracne_ready/GTEx/vst_clean/", hubFilePath="../../ARACNe/hubGenes/hubs.txt") {
  
  ### load library
  if(!require(data.table)) {
    install.packages("data.table")
    library(data.table)
  }
  
  
  ### collect non-bootstrapped network files from the netPath
  f1 <- list.files(netPath)
  
  ### collect expression files from expPath
  f2 <- list.files(expPath)
  
  
  ### load hub list
  hubs <- read.table(file = hubFilePath, sep = "\t", header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
  hubs <- as.character(hubs$V1)
  
  
  ### print result in console for each tissue
  for(i in 1:length(f1)) {
    
    ### load network
    net <- fread(input = paste0(netPath, f1[i], "/", substr(f1[i], 1, nchar(f1[i])-3), ".txt"), sep = "\t", header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
    setDF(net)
    
    ### get interaction nums
    interaction_num <- nrow(net)
    
    ### remove net
    rm(net)
    
    ### load exp
    exp <- read.table(file = paste0(expPath, substr(f1[i], 1, nchar(f1[i])-3), ".dat"), sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
    
    ### extract gene symbols from exp
    exp_genes <- as.character(exp$Gene)
    
    ### remove exp
    rm(exp)
    
    ### get shared hubs between expression and hub list
    shared_hubs <- intersect(exp_genes, hubs)
    
    ### print results in console
    cat(rep("-", 20), "\n")
    writeLines(f1[i])
    writeLines(paste0("(# of Hubs in Exp) x (# of Genes in Exp - 1) = ", length(shared_hubs)*(length(exp_genes)-1)))
    writeLines(paste0("Total Number of Interactions in Network = ", interaction_num))
    cat(rep("-", 20), "\n")
  }
  
  
}

