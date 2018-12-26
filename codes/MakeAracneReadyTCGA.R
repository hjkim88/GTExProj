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
#               > makeAracneReady_TCGA(preprocessedRDAPath="./data/RDA_Files/TCGA_33_RAW_COUNTS.rda",
#                                      outputDir="./results/aracne_ready/TCGA_our_own/")
###

makeAracneReady_TCGA <- function(preprocessedRDAPath="./data/RDA_Files/TCGA_33_RAW_COUNTS.rda",
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
  
  
  ### a function to transfrom Ensembl IDs to Gene symbols
  ensemblIDsToGeneSymbols <- function(ensembl_ids){
    
    ### load library
    if(!require(biomaRt, quietly = TRUE)) {
      if(!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install("biomaRt", version = "3.8")
      require(biomaRt, quietly = TRUE)
    }
    
    ### mapping information between Ensembl ID and Gene symbol
    mart <- useDataset("hsapiens_gene_ensembl", useEnsembl(biomart="ensembl", version=79))
    mart <- getBM( 
      mart = mart, 
      values = ensembl_ids, 
      filter = c("ensembl_gene_id"), 
      attributes = c("ensembl_gene_id", "external_gene_name"),
      verbose = FALSE
    ) 
    
    ### Create a dictionary from ensembl id to gene symbol
    ens_to_gene <- as.character(mart$external_gene_name) 
    names(ens_to_gene) <- as.character(mart$ensembl_gene_id) 
    
    ### return corresponding gene symbols
    return(ens_to_gene[ensembl_ids])
    
  }
  
  
  ### a function returns logical value whether a given row has NA, 0 or 1 across all samples
  isRubbish <- function(geneRow) {
    isZeroOrOne <- TRUE
    for(i in 1:length(geneRow)) {
      if((!is.na(geneRow[i])) && (geneRow[i] > 1)) {
        isZeroOrOne <- FALSE
        break
      }
    }
    
    return(isZeroOrOne)
  }
  
  
  ### A function to transform RNA-Seq data with VST in DESeq2 package
  normalizeRNASEQwithVST <- function(readCount, filtering=TRUE) {
    
    ### load library
    if(!require(DESeq2, quietly = TRUE)) {
      if(!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install("DESeq2", version = "3.8")
      require(DESeq2, quietly = TRUE)
    }
    
    ### make a design matrix for DESeq2 data
    condition <- data.frame(factor(rep("OneClass", ncol(readCount))))
    
    ### Data preparation for DESeq2 format
    deSeqData <- DESeqDataSetFromMatrix(countData=readCount, colData=condition, design= ~0)
    
    if(filtering == TRUE) {
      ### Remove rubbish rows - this will decrease the number of rows
      deSeqData <- deSeqData[rowSums(counts(deSeqData))>1,]
    }
    
    ### VST
    vsd <- vst(deSeqData)
    transCnt <- data.frame(assay(vsd), check.names = FALSE)
    
    return (transCnt)
    
  }
  
  
  ### mapping information between Gene Symbol and Entrez ID (NCBI ID)
  map_symbol_eg <- mappedkeys(org.Hs.egSYMBOL2EG)
  list_symbol2eg <- as.list(org.Hs.egSYMBOL2EG[map_symbol_eg])
  
  
  ### mapping information between Ensembl ID and Entrez ID (NCBI ID)
  map_ensembl_eg <- mappedkeys(org.Hs.egENSEMBL2EG)
  list_ensembl2eg <- as.list(org.Hs.egENSEMBL2EG[map_ensembl_eg])
  
  
  ### make Aracne-ready files for each TCGA tissue
  for(i in 1:length(rcnt_matNames)) {
    ### get raw counts for the given tissue
    df <- get(rcnt_matNames[i])
    
    ### remove FFPE samples
    isFFPE <- which(tcga_sample_info[colnames(df), "is_derived_from_ffpe"] == "YES")
    if(length(isFFPE) > 0) {
      df <- df[,-isFFPE]
    }
    
    ### If the number of samples exceeds 200, 200 samples will be selected based on RIN
    if(ncol(df) > 200) {
      ### collect RIN info
      rin <- tcga_sample_info[colnames(df), "RIN"]
      names(rin) <- colnames(df)
      
      ### order RIN in descending order
      rin <- rin[order(-as.numeric(rin))]
      
      ### only keep the top 200 samples with the highest RIN
      df <- df[,names(rin)[1:200]]
    }
    
    ### if the number of samples is less than 100, we don't run Aracne for that tissue
    if(ncol(df) >= 100) {
      ### get Gencode ID transcript version cleaned
      ensemblIDs <- sapply(rownames(df), function(x) strsplit(x, ".", fixed = TRUE)[[1]][1])
      
      ### annotate gene symbols to the raw counts
      df <- data.frame(Ensembl_ID=ensemblIDs,
                       Gene_Symbol=ensemblIDsToGeneSymbols(ensemblIDs),
                       df,
                       stringsAsFactors = FALSE, check.names = FALSE)
      
      ### remove rows that do not have gene symbols
      naIdx <- which(is.na(df$Gene_Symbol))
      if(length(naIdx) > 0) {
        df <- df[-naIdx,]
      }
      
      ### duplicated gene symbol values
      dups <- unique(df$Gene_Symbol[duplicated(df$Gene_Symbol)])
      
      ### initialize indices which should be retained
      retain <- rep(TRUE, nrow(df))
      
      ### add up all the expressions for one gene symbol and remove the other rows
      for(dup in dups){
        ind <- which(df$Gene_Symbol == dup)
        df[ind[1],3:ncol(df)] <- apply(df[ind,3:ncol(df)], 2, sum)
        retain[setdiff(ind, ind[1])] <- FALSE
      }
      
      ### remove the other rows
      df <- df[retain,]
      
      ### get Entrez IDs correspond to the gene symbols
      df <- data.frame(Entrez_ID=as.character(list_symbol2eg[df$Gene_Symbol]),
                       df,
                       stringsAsFactors = FALSE, check.names = FALSE)
      
      ### remove rows with NULL Entrez_ID
      nullIdx <- union(which(df$Entrez_ID == "NULL"), which(is.null(df$Entrez_ID)))
      if(length(nullIdx) > 0) {
        df <- df[-nullIdx,]
      }
      
      ### get the accurate Entrez ID using Ensembl ID when there are multiple Entrez IDs mapped to one gene symbol
      entrez_dups <- grep("c", df$Entrez_ID)
      for(dup in entrez_dups) {
        df$Entrez_ID[dup] <- list_ensembl2eg[df$Ensembl_ID[dup]][[1]]
      }
      
      ### order based on Entrez ID
      df <- df[order(as.numeric(as.character(df$Entrez_ID))),]
      
      ### cleaning - remove rubbish rows (that have NA, 0 or 1 across all samples)
      df <- df[!apply(df[,4:ncol(df)], 1, isRubbish),]
      
      ### vst-normalization
      df[,4:ncol(df)] <- normalizeRNASEQwithVST(df[,4:ncol(df)], filtering = FALSE)
      
      ### write out the result
      write.table(data.frame(Gene=df$Entrez_ID, df[,4:ncol(df)],
                             stringsAsFactors = FALSE, check.names = FALSE),
                  paste0(outputDir,
                         paste(strsplit(rcnt_matNames[i],
                                        split = "_", fixed = TRUE)[[1]][2:3],
                               collapse = "_"),
                         ".dat"),
                  sep = "\t", row.names = FALSE, quote = FALSE)
    }
  }
  
}
