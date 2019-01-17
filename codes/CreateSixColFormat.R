###
#   File name : CreateSixColFormat.R
#   Author    : Hyunjin Kim
#   Date      : Dec 27, 2017
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Create 6-col format from 4-col format
#
#   Instruction
#               1. Source("CreateSixColFormat.R")
#               2. Run the function "createNew()" - specify the input path () and output directory
#               3. The ARACNe scripts will be generated in the output directory
#
#   Example
#               > source("The_directory_of_CreateSixColFormat.R/CreateSixColFormat.R")
#               > createNew(expDir="./results/aracne_ready/GTEx/vst_clean/",
#                           networkDir="./results/Aracne/GTEx2/")
###

createNew <- function(expDir="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/ExpressionMatrices/GTEx_processed/aracne_ready_counts/vst_clean/",
                      networkDir="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/Aracne2/") {
  
  ### load "consolidateNetworks.R"
  source("./codes/consolidateNetworks.R")
  
  
  ### collect files from the cntPath
  f <- list.files(networkDir)
  f <- f[which(f != "MI")]
  
  
  ### iteratively perform conslidateNets() on the datasets
  for(i in 1:length(f)) {
    
    ### set necessary input values
    gexPath <- paste0(expDir, f[i], ".dat")
    netPath <- c(paste0(networkDir, f[i], "/", f[i], "_tf2.tsv"), paste0(networkDir, f[i], "/", f[i], "_cotf2.tsv"), paste0(networkDir, f[i], "/", f[i], "_signal2.tsv"))
    resultPath <- paste0(networkDir, f[i], "/")
    resultName <- paste0(f[i], "_6cols2")
    
    ### consolidate the networks
    consolidateNets(gexPath, netPath, resultPath, resultName)
  }
  
}

