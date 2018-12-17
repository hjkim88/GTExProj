###
#   File name : OverwriteMIs.R
#   Author    : Hyunjin Kim
#   Date      : Mar 19, 2018
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Overwrite MIs in bootstrapped network with right ones
#
#   Instruction
#               1. Source("OverwriteMIs.R")
#               2. Run the function "overWrite()" - specify the input directories (bootstrap and no-bootstrap)
#               3. The overwritten networks will be generated in the bootstrapped input directories
#
#   Example
#               > source("The_directory_of_OverwriteMIs.R/OverwriteMIs.R")
#               > overWrite(bsnPath="./results/Aracne/GTEx2/"
#                           nbsnPath="./results/Aracne/GTEx2/MI/")
###

overWrite <- function(bsnPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/Aracne2/",
                      nbsnPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/Aracne2/MI/") {
  
  ### load library
  if(!require(data.table)) {
    install.packages("data.table")
    library(data.table)
  }
  
  ### collect bootstrapped network files from the bsnPath
  f1 <- list.files(bsnPath)
  f1 <- f1[which(endsWith(f1, "clean_vst"))]
  
  ### iteratively overwrite
  for(i in 1:length(f1)) {
    
    ### load bootstrapped networks
    net1 <- read.table(file = paste0(bsnPath, f1[i], "/", f1[i], "_tf.tsv"), header = FALSE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
    net2 <- read.table(file = paste0(bsnPath, f1[i], "/", f1[i], "_cotf.tsv"), header = FALSE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
    net3 <- read.table(file = paste0(bsnPath, f1[i], "/", f1[i], "_signal.tsv"), header = FALSE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
    
    ### load non-bootstrapped network
    net4 <- fread(input = paste0(nbsnPath, f1[i], "_mi/", f1[i], ".txt"), sep = "\t", header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
    setDF(net4)
    
    ### preprocess the networks
    rownames(net1) <- paste(net1$V1, net1$V2)
    rownames(net2) <- paste(net2$V1, net2$V2)
    rownames(net3) <- paste(net3$V1, net3$V2)
    rownames(net4) <- paste(net4$V1, net4$V2)
    
    ### overwrite & save
    net1$V3 <- as.numeric(net4[rownames(net1), 3])
    write.table(net1, file = paste0(bsnPath, f1[i], "/", f1[i], "_tf2.tsv"), sep = "\t", row.names = FALSE, col.names = FALSE)
    net2$V3 <- as.numeric(net4[rownames(net2), 3])
    write.table(net2, file = paste0(bsnPath, f1[i], "/", f1[i], "_cotf2.tsv"), sep = "\t", row.names = FALSE, col.names = FALSE)
    net3$V3 <- as.numeric(net4[rownames(net3), 3])
    write.table(net3, file = paste0(bsnPath, f1[i], "/", f1[i], "_signal2.tsv"), sep = "\t", row.names = FALSE, col.names = FALSE)
    
    writeLines(paste(i, "/", length(f1)))
    
  }
  
}
