###
#   File name : GTExFiltering2.R
#   Author    : Hyunjin Kim
#   Date      : Nov 10, 2017
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Cleaning up the raw count data (remove rows other than genes)
#
#   Instruction
#               1. Source("GTExFiltering2.R")
#               2. Run the function "filtering()" - specify the input file (raw count) and output directory
#               3. The filtered raw counts will be generated in the output path
#
#   Example
#               > source("The_directory_of_GTExFiltering2.R/GTExFiltering2.R")
#               > filtering(rawCntPath="./data/raw_counts/",
#                           outputPath="./results/filtered_counts/")
###

filtering <- function(rawCntPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/ExpressionMatrices/GTEx_original/Counts/",
                      outputPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/ExpressionMatrices/GTEx_processed/GeneCounts/") {
  
  ### load library
  if(!require(org.Hs.eg.db)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("org.Hs.eg.db")
    library(org.Hs.eg.db)
  }
  
  
  ### collect files from the rawCntPath
  f <- list.files(rawCntPath)
  
  
  ### A function to iteratively perform filtering process
  filterData <- function(dataPath) {
    
    ### read raw counts
    rawData <- read.table(dataPath, sep="\t", header=TRUE, stringsAsFactors = FALSE, check.names = FALSE)
    colnames(rawData)[1] <- "Gencode_ID"
    colnames(rawData)[2] <- "Gene_Symbol"
    
    ### get gene names and clean the transcript versions
    gene_names <- apply(data.frame(rawData$Gencode_ID), 1, function(x) strsplit(x, ".", fixed = TRUE)[[1]][1])
    
    ### mapping information between Ensembl ID and Entrez ID (NCBI ID)
    map_ensembl_eg <- mappedkeys(org.Hs.egENSEMBL2EG)
    list_ensembl2eg <- as.list(org.Hs.egENSEMBL2EG[map_ensembl_eg])
    
    ### initialize entrez_id
    entrez_id <- rep(NA, length(gene_names))
    names(entrez_id) <- gene_names
    
    ### initialize gene_idx
    gene_idx <- c(1:length(gene_names))
    names(gene_idx) <- gene_names
    
    ### iteratively change Ensembl ID to Entrez ID
    for(i in 1:length(entrez_id)) {
      
      ### get Entrez ID
      temp <- list_ensembl2eg[gene_names[i]][[1]][1]
      
      ### If there is corresponding Entrez ID then change it. Otherwise, NA
      if(!is.null(temp)) {
        ### Change the Ensembl ID to Entrez ID
        entrez_id[i] <- temp
        
        ### If there are more than one Entrez IDs for one Ensembl ID, expand the list
        if(length(list_ensembl2eg[gene_names[i]][[1]]) > 1) {
          for(j in 2:length(list_ensembl2eg[gene_names[i]][[1]])) {
            entrez_id <- c(entrez_id, list_ensembl2eg[gene_names[i]][[1]][j])
            names(entrez_id)[length(entrez_id)] <- gene_names[i]
            
            gene_idx <- c(gene_idx, i)
            names(gene_idx)[length(gene_idx)] <- gene_names[i]
          }
        }
      }
    }
    
    ### remove NAs
    gene_idx <- gene_idx[!is.na(entrez_id)]
    entrez_id <- entrez_id[!is.na(entrez_id)]
    
    ### A function implemented by Aris Floratos
    ### Retrun the co-efficient of variation
    cv <- function(x){
      #First, remove 5% of outliers
      PERC = 0.05
      x = x[order(x, decreasing = TRUE)[(as.integer(PERC*length(x))+1):length(x)]]
      m = mean(x)
      if(m == 0)
        return(0)
      else
        return(sd(x)/m)
    }
    
    
    ### remove duplicates
    dups <- unique(entrez_id[duplicated(entrez_id)])
    retain <- rep(TRUE, length(entrez_id))
    for(dup in dups){
      ind <- which(entrez_id == dup)
      cvs <- sapply(ind, function(x) cv(as.numeric(rawData[x,3:ncol(rawData)])))
      ind <- ind[-(which(cvs == max(cvs))[1])]
      retain[ind] = FALSE
    }
    
    ### make non-duplicated entrez_id & gene_idx
    entrez_id <- entrez_id[retain]
    gene_idx <- gene_idx[retain]
    
    ### filter the raw data
    filteredData <- rawData[gene_idx,]
    rownames(filteredData) <- entrez_id
    
    ### sort the rows based on Entrez ID
    filteredData <- filteredData[order(rownames(filteredData)),]
    
    return(filteredData)
  }
  
  
  ### iteratively perform filterData() function for all the raw count data in the input path
  for(i in 1:length(f)) {
    
    ### filter rows with filterData() function
    result <- filterData(paste0(rawCntPath, f[i]))
    
    ### extract file name
    fileName <- substr(f[i], 1, nchar(f[i])-4)
    
    ### save the filtered data as tab-separated file
    write.table(result, paste0(outputPath, fileName, ".txt"), sep = "\t")
  }
  
}



