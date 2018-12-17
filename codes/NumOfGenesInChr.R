###
#   File name : NumOfGenesInChr.R
#   Author    : Hyunjin Kim
#   Date      : Aug 13, 2018
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Count gene numbers in each chromosome
#
#   Instruction
#               1. Source("NumOfGenesInChr.R")
#               2. Run the function "countGeneNum()" - specify the input (gtf file) and output file Path
#               3. A RDA file that contains the number of genes of all the chromsomes will be generated
#
#   Example
#               > source("The_directory_of_NumOfGenesInChr.R/NumOfGenesInChr.R")
#               > countGeneNum(gtfPath="E:/Reference/hg19.gtf",
#                              outputPath="./Gene_num_chr.rda")
###

countGeneNum <- function(gtfPath="E:/Reference/hg19.gtf",
                         outputPath="./Gene_num_chr.rda") {
  
  ### load library
  if(!require(rtracklayer)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("rtracklayer")
    library(rtracklayer)
  }
  
  
  ### load data
  gtf <- import(gtfPath)
  gtf <- as.data.frame(gtf)
  
  
  ### get chromosome info
  chrs <- as.character(gtf$seqnames[which(gtf$type == "gene")])
  
  
  ### create an empty vector
  geneNum <- rep(0, 23)
  names(geneNum) <- paste0(rep("chr", length(geneNum)), 1:length(geneNum))
  
  
  ### count number of genes
  for(i in 1:length(chrs)) {
    x <- substr(chrs[i], 4, nchar(chrs[i]))
    if(x == "X" || x == "Y") {
      geneNum[23] <- geneNum[23] + 1
    } else if(x == "M") {
    } else {
      geneNum[as.integer(x)] <- geneNum[as.integer(x)] + 1
    }
  }
  
  
  ### README function
  README <- function() {
    writeLines("The number of genes in each chromosome")
    writeLines("Only \"gene\" type rows in the gtf file were counted")
  }
  
  
  ### save the result
  save(list = c("geneNum", "README"), file = outputPath)
  
}


