###
#   File name : MakeAracneReadyTCGA.R
#   Author    : Hyunjin Kim
#   Date      : Dec 22, 2018
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : We downloaded raw counts of all the TCGA tissues and pre-processed the counts with
#               "PreprocessTCGA.R". We have raw counts and their sample info including RIN and FFPE.
#               The process to make Aracne-ready files are like below:
#               1. If there are duplicated gene symbols, then just add all counts up
#               2. Remove genes that have NA, 0 or 1 across all samples
#               3. Remove FFPE samples
#               4. Perform VST(Variance Stabilizing Transformation) on counts
#               5. Only use tissues with >= 100 samples and if there are more than 200 samples,
#                  choose 200 with highest RIN
#
#   Instruction
#               1. Source("MakeAracneReadyTCGA.R")
#               2. Run the function "makeAracneReady_TCGA()" - specify the pre-processed TCGA RDA file path and output directory
#               3. The Aracne-ready TCGA text files will be generated in the output file directory
#
#   Example
#               > source("The_directory_of_MakeAracneReadyTCGA.R/MakeAracneReadyTCGA.R")
#               > makeAracneReady_TCGA(preprocessedRDAPath="./data/RDA_Files/TCGA_RAW_COUNTS.rda",
#                                      outputDir="./results/aracne_ready/TCGA_our_own/")
###

makeAracneReady_TCGA <- function(preprocessedRDAPath="./data/RDA_Files/TCGA_RAW_COUNTS.rda",
                                 outputDir="./results/aracne_ready/TCGA_our_own/") {
  
  ### load library
  if(!require(org.Hs.eg.db, quietly = TRUE)) {
    if(!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("org.Hs.eg.db", version = "3.8")
    require(org.Hs.eg.db, quietly = TRUE)
  }
  
  
  ### load the pre-processed TCGA RDA file
  load(preprocessedRDAPath)
  
  
  
}



