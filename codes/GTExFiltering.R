###
#   File name : GTExFiltering.R
#   Author    : Hyunjin Kim
#   Date      : Nov 2, 2017
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Cleaning up the raw count data (remove rows other than genes)
#
#   Instruction
#               1. Source("GTExFiltering.R")
#               2. Run the function "filtering()" - specify the input file (raw count) and output directory
#               3. The filtered raw counts will be generated in the output path
#
#   Example
#               > source("The_directory_of_GTExFiltering.R/GTExFiltering.R")
#               > filtering(rawCntPath="./data/raw_counts/", annotPath="./data/gencode.v19.annotation.gtf", outputPath="./results/filtered_counts/")
###


filtering <- function(rawCntPath="./data/raw_counts/", annotPath="./data/gencode.v19.annotation.gtf", outputPath="./results/filtered_counts/") {
  
  ### JAVA memory setup - because of xlsx package
  #options(java.parameters = "-Xmx8000m")
  
  ### JAVA garbage collection calling
  #jgc <- function()
  #{
  #  gc()
  #  .jcall("java/lang/System", method = "gc")
  #}
  
  
  ### load library
  if(!require(rtracklayer)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("rtracklayer")
    library(rtracklayer)
  }
  #if(!require(xlsx)) {
  #  install.packages("xlsx")
  #  library(xlsx)
  #}
  
  
  ### load and preprocess annotation file (.GTF)
  annot <- readGFF(annotPath, tags = c("gene_id", "gene_type", "gene_name"))
  annot <- annot[,-which(colnames(annot) %in% c("seqid", "source", "start", "end", "score", "strand", "phase"))]
  ### choose rows which have type "gene"
  annot <- annot[which(annot$type == "gene"),]
  rownames(annot) <- annot$gene_id
  annot <- annot[,-union(which(colnames(annot) == "type"), which(colnames(annot) == "gene_id"))]
  
  
  ### collect files from the rawCntPath
  f <- list.files(rawCntPath)
  
  
  ### A function to iteratively perform filtering process
  filterData <- function(dataPath) {
    ### read raw counts
    rawData <- read.table(dataPath, sep="\t", header=TRUE, stringsAsFactors = FALSE, check.names = FALSE)
    
    ### indicies for protein-coding gene only
    idx <- which(annot[rawData$Name,"gene_type"] == "protein_coding")
    
    ### choose rows which are protein-coding genes
    filteredData <- rawData[idx,]
    rownames(filteredData) <- filteredData$Name
    colnames(filteredData)[which(colnames(filteredData) == "Description")] <- "Gene_Symbol"
    filteredData <- filteredData[,-which(colnames(filteredData) == "Name")]
    
    return(filteredData)
  }
  
  
  ### iteratively perform filterData() function for all the raw count data in the input path
  for(i in 1:length(f)) {
    ### garbage collection because of JAVA heap memory
    #jgc()
    
    ### filter rows with filterData() function
    result <- filterData(paste0(rawCntPath, f[i]))
    
    ### extract file name
    fileName <- substr(f[i], 1, nchar(f[i])-4)
    
    
    ### save the filtered data as csv file
    write.csv(result, paste0(outputPath, fileName, ".csv"))
    
    ### save the filtered data as excel file
    #write.xlsx2(result, paste0(outputPath, fileName, ".xlsx"), fileName)
  }
  
}



