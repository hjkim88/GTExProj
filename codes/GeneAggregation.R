###
#   File name : GeneAggregation.R
#   Author    : Hyunjin Kim
#   Date      : Nov 13, 2017
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : If there are duplicated gene symbols, then just add all counts up
#
#   * First, the data already has Genecode IDs (Ensembl ID) and corresponding gene symbols
#     Several Genecode IDs can be mapped to one gene symbol, so there may be duplicated gene symbols
#     We keep unique gene symbols only, so add up all the counts for any one gene symbol
#     So now there are no duplicated gene symbols
#     Secondly, try to transform the gene symbols to Entrez IDs
#     Very few gene symbols have multiple Entrez IDs mapped
#     In that case, if we use the corresponding Genecode ID, we can get accurate Entrez ID for that
#     This is because both many Genecode IDs and Entrez IDs are mapped to one gene symbol
#
#   Instruction
#               1. Source("GeneAggregation.R")
#               2. Run the function "reconcile()" - specify the input file (raw count)
#
#   Example
#               > source("The_directory_of_GeneAggregation.R/GeneAggregation.R")
#               > reconcile(rawCntPath="./data/raw_counts/",
#                           outputPath="./results/aggregated_counts/")
###

reconcile <- function(rawCntPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/ExpressionMatrices/GTEx_original/Counts/",
                      outputPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/ExpressionMatrices/GTEx_processed/GeneCounts/") {
  
  ### load library
  if(!require(org.Hs.eg.db)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("org.Hs.eg.db")
    library(org.Hs.eg.db)
  }
  
  
  ### collect files from the rawCntPath
  f <- list.files(rawCntPath)
  
  
  ### A function for gene aggregation (Unique gene symbol)
  aggregation <- function(dataPath) {
    
    ### read raw counts
    rawData <- read.table(dataPath, sep="\t", header=TRUE, stringsAsFactors = FALSE, check.names = FALSE)
    colnames(rawData)[1] <- "Gencode_ID"
    colnames(rawData)[2] <- "Gene_Symbol"
    
    ### duplicated values
    dups <- unique(rawData$Gene_Symbol[duplicated(rawData$Gene_Symbol)])
    
    ### initialize indices which should be retained
    retain <- rep(TRUE, nrow(rawData))
    
    ### add up all the expressions for one gene symbol and remove the other rows
    for(dup in dups){
      ind <- which(rawData$Gene_Symbol == dup)
      rawData[ind[1],3:ncol(rawData)] <- apply(rawData[ind,3:ncol(rawData)], 2, sum)
      retain[setdiff(ind, ind[1])] <- FALSE
    }
    
    ### remove the other rows
    rawData <- rawData[retain,]
    
    return(rawData)
  }
  
  
  ### A function to iteratively perform gene name transforming process
  convert <- function(agg_data) {
    
    ### mapping information between Gene Symbol and Entrez ID (NCBI ID)
    map_symbol_eg <- mappedkeys(org.Hs.egSYMBOL2EG)
    list_symbol2eg <- as.list(org.Hs.egSYMBOL2EG[map_symbol_eg])
    
    ### mapping information between Ensembl ID and Entrez ID (NCBI ID)
    map_ensembl_eg <- mappedkeys(org.Hs.egENSEMBL2EG)
    list_ensembl2eg <- as.list(org.Hs.egENSEMBL2EG[map_ensembl_eg])
    
    ### get corresponding Entrez IDs
    entrez_id <- as.character(list_symbol2eg[agg_data$Gene_Symbol])
    
    ### group them all
    agg_data <- cbind(Entrez_ID=entrez_id, agg_data)
    
    ### remove rows with NULL Entrez_ID
    filteredData <- agg_data[-which(agg_data$Entrez_ID == "NULL"),]
    
    ### get Gencode ID transcript version cleaned
    gencode <- apply(data.frame(filteredData$Gencode_ID), 1, function(x) strsplit(x, ".", fixed = TRUE)[[1]][1])
    
    ### get the accurate Entrez ID using Ensembl ID when there are multiple Entrez IDs mapped to one gene symbol
    filteredData$Entrez_ID <- as.character(filteredData$Entrez_ID)
    entrez_dups <- grep("c", filteredData$Entrez_ID)
    for(dup in entrez_dups) {
      filteredData$Entrez_ID[dup] <- list_ensembl2eg[gencode[dup]][[1]]
    }
    
    ### remove Genecode IDs
    filteredData <- filteredData[,-which(colnames(filteredData) == "Gencode_ID")]
    
    ### order based on Entrez ID
    filteredData <- filteredData[order(as.numeric(as.character(filteredData$Entrez_ID))),]
    
    return(filteredData)
  }
  
  
  ### iteratively filter raw count data
  for(i in 1:length(f)) {
    
    ### aggregate duplicated gene symbols
    aggregated <- aggregation(paste0(rawCntPath, f[i]))
    
    ### convert gene symbols to NCBI IDs
    result <- convert(aggregated)
    
    ### extract file name
    fileName <- substr(f[i], 1, nchar(f[i])-4)
    
    ### save the filtered data as tab-separated file
    write.table(result, paste0(outputPath, fileName, ".txt"), sep = "\t", row.names = FALSE)
  }
  
}


