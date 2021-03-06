###
#   File name : Preprocess_TCGA_Normal.R
#   Author    : Hyunjin Kim
#   Date      : Jun 26, 2019
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : From the TCGA raw count data, select normal samples only and preprocess them for Aracne run
#
#   Instruction
#               1. Source("Preprocess_TCGABRCA_Normal.R")
#               2. Run the function "tcga_normal()" - specify the neccessary input file paths and output file path
#               3. The preprocessed TCGA Normal files will be generated as the output file path
#
#   Example
#               > source("The_directory_of_Preprocess_TCGA_Normal.R/Preprocess_TCGA_Normal.R")
#               > tcga_normal(preprocessedRDAPath="./data/RDA_Files/TCGA_33_RAW_COUNTS.rda",
#                             outputDir="./results/aracne_ready/TCGA_normal/")
###

tcga_normal <- function(preprocessedRDAPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/TCGA_33_RAW_COUNTS.rda",
                        outputDir="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/TCGA_BRCA_Normal/") {
  
  ### load library
  if(!require(org.Hs.eg.db, quietly = TRUE)) {
    if(!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("org.Hs.eg.db", version = "3.8")
    require(org.Hs.eg.db, quietly = TRUE)
  }
  
  ### load the pre-processed TCGA RDA file
  load(preprocessedRDAPath)
  
  ### retain normal tissues only
  tcga_sample_info <- tcga_sample_info[which(tcga_sample_info[,"Sample Type"] == "Solid Tissue Normal"),]
  
  ### remove FFPE samples
  tcga_sample_info <- tcga_sample_info[which(tcga_sample_info[,"is_derived_from_ffpe"] == "NO"),]
  
  ### order the sample info based on Project ID and Case ID
  tcga_sample_info <- tcga_sample_info[order(tcga_sample_info[,"Project ID"],
                                             tcga_sample_info[,"Case ID"]),]
  
  ### if there are multiple samples per one patient in each tissue, select one with the highest RIN
  unique_tissues <- unique(tcga_sample_info[,"Project ID"])
  rIdx <- NULL
  for(i in 1:length(unique_tissues)) {
    ### get indicies for the given tissue
    tempIdx <- which(tcga_sample_info[,"Project ID"] == unique_tissues[i])
    
    ### get duplicated indicies in the given tissue
    dupIdx <- tempIdx[which(duplicated(tcga_sample_info[tempIdx, "Case ID"]))]
    
    ### if there are duplicates, select one with the highest RIN
    ### tie breaker (multiple highest RIN) - select one with the highest lexical order
    if(length(dupIdx) > 0) {
      dups <- unique(tcga_sample_info[dupIdx, "Case ID"])
      
      ### collect indicies except one that will remain
      ### those indicies will be removed away later
      for(j in 1:length(dups)) {
        dIdx <- which(tcga_sample_info[,"Case ID"] == dups[j])
        rIdx <- c(rIdx, dIdx[order(tcga_sample_info[dIdx, "RIN"])])
        rIdx <- rIdx[-length(rIdx)]
      }
    }
  }
  if(length(rIdx) > 0) {
    tcga_sample_info <- tcga_sample_info[-rIdx,]
  }
  
  ### make all the raw count matrices have samples only appeared in the tcga_sample_info
  rIdx <- NULL
  for(i in 1:length(rcnt_matNames)) {
    ### get raw count matrix for the given tissue
    rcnt_mat <- get(rcnt_matNames[i])
    
    ### retain samples only appeared in the tcga_sample_info (which means they are filtered)
    rcnt_mat <- rcnt_mat[,which(colnames(rcnt_mat) %in% rownames(tcga_sample_info)), drop = FALSE]
    
    if(ncol(rcnt_mat) > 0) {
      ### order the samples in lexical order
      rcnt_mat <- rcnt_mat[,order(colnames(rcnt_mat)), drop = FALSE]
      
      ### change the colnames based on the first 15 characters
      colnames(rcnt_mat) <- substr(colnames(rcnt_mat), 1, 15)
      
      ### save the result back to the variable
      assign(rcnt_matNames[i], pos = 1, rcnt_mat, envir = globalenv())
    } else {
      rIdx <- c(rIdx, i)
      rm(list = rcnt_matNames[i])
    }
  }
  if(length(rIdx) > 0) {
    rcnt_matNames <- rcnt_matNames[-rIdx]
  }
  
  ### change the row names based on the first 15 characters
  rownames(tcga_sample_info) <- substr(rownames(tcga_sample_info), 1, 15)
  
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
  
  # ### Test code to see which tissues have samples more than 100
  # for(i in 1:length(rcnt_matNames)) {
  #   writeLines(paste(rcnt_matNames[i], ncol(get(rcnt_matNames[i]))))
  # }
  
  ### make Aracne-ready files for each TCGA normal tissue
  for(i in 1:length(rcnt_matNames)) {
    ### get raw counts for the given tissue
    df <- get(rcnt_matNames[i], pos = 1)
    
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
      
      ### set rownames with Entrez ID
      rownames(df) <- df$Entrez_ID
      
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
                         "_normal.dat"),
                  sep = "\t", row.names = FALSE, quote = FALSE)
    }
  }
  
}
