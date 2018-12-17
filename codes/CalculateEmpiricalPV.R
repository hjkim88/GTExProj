###
#   File name : CalculateEmpiricalPV.R
#   Author    : Hyunjin Kim
#   Date      : Mar 28, 2018
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Calculate empirical p-values of no-dpi networks
#
#   Instruction
#               1. Source("CalculateEmpiricalPV.R")
#               2. Run the function "calEmpPV()" - specify the Aracne network file
#               3. Empirical p-values will be added and new Aracne network file will be generated in the same path as the input
#
#   Example
#               > source("The_directory_of_CalculateEmpiricalPV.R/CalculateEmpiricalPV.R")
#               > calEmpPV(inputPath="./results/Aracne/GTEx2/MI/Colon-Sigmoid_clean_vst_mi/Colon-Sigmoid_clean_vst.txt")
###


calEmpPV <- function(inputPath="./results/Aracne/GTEx2/MI/Colon-Sigmoid_clean_vst_mi/Colon-Sigmoid_clean_vst.txt") {
  
  ### load library
  if(!require(data.table)) {
    install.packages("data.table")
    library(data.table)
  }
  
  ### load the network file
  net <- fread(input = inputPath, sep = "\t", header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
  setDF(net)
  
  ### the fastest approach to calculate empirical p-values
  # order the network based on MIs in ascending order
  net <- net[order(net$V3),]
  # give a rank based on MI (smaller rank = larger MI)
  rowNum <- nrow(net)
  net$V4 <- rev(seq(1:rowNum))
  # make an empty MI-pValue map 
  mi_pv_map <- rep(NA, length(unique(net$V3)))
  # remove duplicates (the most important part) - same MI score = larger rank
  # e.g., If the largest MIs are tied among 3, then the rank of them are all 3, not 1
  # e.g., If all the MIs are all the same, then their p-vlaues are all 1, not all 0
  mi_pv_map <- net$V4[which(!duplicated(net$V3))]
  # calculate empirical p-value for the unique ones
  mi_pv_map <- mi_pv_map / rowNum
  # set names for the map
  names(mi_pv_map) <- as.character(unique(net$V3))
  # get p-values from the map
  net$V4 <- mi_pv_map[as.character(net$V3)]
  # order the network based on p-values in ascending order
  net <- net[order(net$V4),]
  
  ### save the result
  fName <- paste0(substr(inputPath, 1, nchar(inputPath)-4), "_pv.txt")
  fwrite(net, file = fName, sep = "\t", row.names = FALSE, col.names = FALSE)
  
}

