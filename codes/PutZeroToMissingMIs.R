###
#   File name : PutZeroToMissingMIs.R
#   Author    : Hyunjin Kim
#   Date      : Mar 26, 2018
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Put zeros for the omitted TF-target pairs (They are from missing interactions and MIs of them are 0s)
#
#   Instruction
#               1. Source("PutZeroToMissingMIs.R")
#               2. Run the function "putZeroToMissingMIs()" - specify the MI-missing network and output file path
#               3. The corrected network file will be generated in output file path
#
#   Example
#               > source("The_directory_of_PutZeroToMissingMIs.R/PutZeroToMissingMIs.R")
#               > putZeroToMissingMIs(netPath="./results/Aracne/GTEx2/Colon-Transverse_clean_vst/Colon-Transverse_clean_vst_tf2.tsv",
#                                     outputPath="./results/Aracne/GTEx2/Colon-Transverse_clean_vst/Colon-Transverse_clean_vst_tf3.tsv")
###


putZeroToMissingMIs <- function(netPath="./results/Aracne/GTEx2/Colon-Transverse_clean_vst/Colon-Transverse_clean_vst_tf2.tsv",
                                outputPath="./results/Aracne/GTEx2/Colon-Transverse_clean_vst/Colon-Transverse_clean_vst_tf3.tsv") {
  
  ### load network
  net <- read.table(file = netPath, sep = "\t", header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
  colnames(net) <- c("TF", "Target", "MI", "PV")
  
  ### get NA indicies
  naIdx <- which(is.na(net$MI))
  
  ### set NA MIs as 0
  net$MI[naIdx] <- 0
  
  ### write the corrected result
  write.table(net, file = outputPath, sep = "\t", row.names = FALSE, col.names = FALSE)
  
}

