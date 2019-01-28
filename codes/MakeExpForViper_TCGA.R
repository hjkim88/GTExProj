###
#   File name : MakeExpForViper_TCGA.R
#   Author    : Hyunjin Kim
#   Date      : Jan 25, 2018
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Load all the non-cleaned counts of all the tissues and normalize them all together
#
#   Instruction
#               1. Source("MakeExpForViper_TCGA.R")
#               2. Run the function "norm_counts()" - specify the input file directories (Aracne networks & expressions) and output path
#               3. TCGA expression matrices (and their combined RDA file) will be generated
#
#   Example
#               > source("The_directory_of_MakeExpForViper_TCGA.R/MakeExpForViper_TCGA.R")
#               > norm_counts(expPath="./data/RDA_Files/TCGA_33_RAW_COUNTS.rda",
#                             outputPath="./data/RDA_Files/TCGA_26_EMat_Viper.rda")
###

norm_counts <- function(expPath="./data/RDA_Files/TCGA_33_RAW_COUNTS.rda",
                        outputPath="./data/RDA_Files/TCGA_26_EMat_Viper.rda") {
  
  ### load library
  if(!require(DESeq2)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("DESeq2")
    library(DESeq2)
  }
  if(!require(org.Hs.eg.db, quietly = TRUE)) {
    if(!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("org.Hs.eg.db", version = "3.8")
    require(org.Hs.eg.db, quietly = TRUE)
  }
  
  
  ### load expression RDA
  load(expPath)
  
  
  ### retain only primary tumor & new primary tumor
  tcga_sample_info <- tcga_sample_info[union(union(which(tcga_sample_info[,"Sample Type"] == "Primary Tumor"),
                                                   which(tcga_sample_info[,"Sample Type"] == "Additional - New Primary")),
                                             which(tcga_sample_info[,"Sample Type"] == "Primary Blood Derived Cancer - Peripheral Blood")),]
  
  
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
  tcga_sample_info <- tcga_sample_info[-rIdx,]
  
  
  ### make all the raw count matrices have samples only appeared in the tcga_sample_info
  for(i in 1:length(rcnt_matNames)) {
    ### get raw count matrix for the given tissue
    rcnt_mat <- get(rcnt_matNames[i])
    
    ### retain samples only appeared in the tcga_sample_info (which means they are filtered)
    rcnt_mat <- rcnt_mat[,which(colnames(rcnt_mat) %in% rownames(tcga_sample_info))]
    
    ### order the samples in lexical order
    rcnt_mat <- rcnt_mat[,order(colnames(rcnt_mat))]
    
    ### change the colnames based on the first 15 characters
    colnames(rcnt_mat) <- substr(colnames(rcnt_mat), 1, 15)
    
    ### save the result back to the variable
    assign(rcnt_matNames[i], rcnt_mat, envir = globalenv())
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
  
  
  ### mapping information between Gene Symbol and Entrez ID (NCBI ID)
  map_symbol_eg <- mappedkeys(org.Hs.egSYMBOL2EG)
  list_symbol2eg <- as.list(org.Hs.egSYMBOL2EG[map_symbol_eg])
  
  
  ### mapping information between Ensembl ID and Entrez ID (NCBI ID)
  map_ensembl_eg <- mappedkeys(org.Hs.egENSEMBL2EG)
  list_ensembl2eg <- as.list(org.Hs.egENSEMBL2EG[map_ensembl_eg])
  
  
  ### initialize the count matrix names for a RDA file
  emat_tcga_names <- NULL
  
  
  ### make Aracne-ready files for each TCGA tissue
  for(i in 1:length(rcnt_matNames)) {
    ### get raw counts for the given tissue
    df <- get(rcnt_matNames[i])
    
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
      
      ### save the Aracne-ready expressions to a variable
      emat_tcga_names <- c(emat_tcga_names, paste0("emat_", paste(strsplit(rcnt_matNames[i],
                                                                                 split = "_", fixed = TRUE)[[1]][2:3],
                                                                        collapse = "_")))
      assign(emat_tcga_names[length(emat_tcga_names)], df[,4:ncol(df)], envir = globalenv())
    }
    
    ### remove the input raw count
    rm(list = c(rcnt_matNames[i]))
    
    ### garbage collection
    gc()
  }
  
  
  ### get number of samples of all the tissues
  emat_tcga_sampleNum <- 0
  for(i in 1:length(emat_tcga_names)) {
    emat_tcga_sampleNum[i] <- ncol(get(emat_tcga_names[i]))
  }
  
  
  ### combine all the tissues
  d <- get(emat_tcga_names[1])
  for(i in 2:length(emat_tcga_names)) {
    d <- cbind(d, get(emat_tcga_names[i]))
  }
  
  
  ### A function to transform RNA-Seq data with VST in DESeq2 package
  ### input readCount of this function includes gene symbols in the first column
  RNASEQwithVST <- function(readCount) {
    
    ### make a design matrix for DESeq2 data
    condition <- data.frame(factor(rep("OneClass", ncol(readCount))))
    
    ### Data preparation for DESeq2 format
    deSeqData <- DESeqDataSetFromMatrix(countData=readCount, colData=condition, design= ~0)
    
    ### VST
    vsd <- vst(deSeqData)
    transCnt <- data.frame(assay(vsd))
    
    ### copy column names from original count data
    colnames(transCnt) <- colnames(readCount)
    
    return (transCnt)
  }
  
  
  ### perform VST
  assign("emat_tcga_all", as.matrix(RNASEQwithVST(d)), envir = globalenv())
  
  
  ### split the combined one into tissue-specific data
  for(i in 1:length(emat_tcga_names)) {
    assign(emat_tcga_names[i], emat_tcga_all[,(sum(emat_tcga_sampleNum[0:(i-1)])+1):sum(emat_tcga_sampleNum[0:i])], envir = globalenv())
  }
  
  
  ### save the results
  vars <- c(emat_tcga_names, "emat_tcga_names")
  save(list = vars, file = outputPath)
  
}

