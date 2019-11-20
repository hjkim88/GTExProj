###
#   File name : TCGAMutInfoPerSample.R
#   Author    : Hyunjin Kim
#   Date      : Nov 20, 2019
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : There are preprocessed TCGA MAF files. Use maftools package to extract
#               mutation information from the files for every sample from each TCGA tissue
#
#   Instruction
#               1. Source("TopMutatedRegulons.R")
#               2. Run the function "get_tcga_mutation_info" - specify the necessary input and output paths
#               3. The result RDA file will be generated in the output path
#
#   Example
#               > source("The_directory_of_TCGAMutInfoPerSample.R/TCGAMutInfoPerSample.R")
#               > get_tcga_mutation_info(mafFileDir="C:/Research/CUMC/GTExProj/data/TCGA/MAF/preprocessed/",
#                                        outputPath="C:/Research/CUMC/GTExProj/data/RDA_Files/TCGA_33_maf_per_sample.rda")
###

get_tcga_mutation_info <- function(mafFileDir="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/TCGA/MAF/",
                                   outputPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/TCGA/RDA_Files/TCGA_33_maf_per_sample.rda") {
  
  ### load libraries
  if(!require(maftools, quietly = TRUE)) {
    if(!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("maftools", version = "3.8")
    require(maftools, quietly = TRUE)
  }
  if(!require("annotate", quietly = TRUE)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("annotate")
    require("annotate", quietly = TRUE)
  }
  if(!require("org.Hs.eg.db", quietly = TRUE)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("org.Hs.eg.db")
    require("org.Hs.eg.db", quietly = TRUE)
  }
  
  ### get maf file list
  f <- list.files(mafFileDir, pattern = "\\.maf$")
  
  ### get tcga tissue names
  tissue_names <- tolower(sapply(f, function(x) strsplit(x, ".", TRUE)[[1]][1], USE.NAMES = FALSE))
  
  # *****************************************************************************
  #
  # Map Entrez Gene IDs to Gene Symbols
  #
  # geneID:	A single integer, or a vector of integers, or a single strings, or
  #			a vector of strings. In each case, integers or strings are 
  #			Entrez gene ids.
  # 
  # Returns a vector of the same size as "geneID" where the i-th entry is the 
  # gene symbol corresponsing to the i-th Entrez id in "geneID". The return 
  # vector entries also are named using the Entrez Ids. 
  #
  # -----------------------------------------------------------------------------
  # -------FIXME: Modify the code to use "select", as in method geneSymbolToEntrezId
  # -----------------------------------------------------------------------------
  # *****************************************************************************
  entrezIDtoSymbol <- function(geneID){
    if (is.numeric(geneID))
      geneID = as.character(geneID)
    return(getSYMBOL(geneID, data="org.Hs.eg"))
  }
  
  ### for every tissue, make mutation info and save it to the list
  tcga_maf_per_sample <- vector("list", length(tissue_names))
  names(tcga_maf_per_sample) <- tissue_names
  for(tissue_name in tissue_names) {
    
    ### get maf file name
    maf_file <- paste0(toupper(tissue_name), ".maf")
    
    ### load total MAF file
    maf <- read.maf(maf = paste0(mafFileDir, maf_file))
    
    ### make sample-gene mutation info
    ### here, it is to construct a data frame that has rows as
    ### unique gene mutation - sample pairs
    ### and columns as gene mutation info
    mutation_info <- data.frame(maf@data[which(!duplicated(maf@data[,c("Entrez_Gene_Id", "Tumor_Sample_Barcode")])),
                                         c("Tumor_Sample_Barcode", "Entrez_Gene_Id")],
                                stringsAsFactors = FALSE, check.names = FALSE)
    mutation_info$Tumor_Sample_Barcode <- as.character(mutation_info$Tumor_Sample_Barcode)
    mutation_dups <- data.frame(maf@data[which(duplicated(maf@data[,c("Entrez_Gene_Id", "Tumor_Sample_Barcode")])),
                                         c("Tumor_Sample_Barcode", "Entrez_Gene_Id")],
                                stringsAsFactors = FALSE, check.names = FALSE)
    mutation_dups$Tumor_Sample_Barcode <- as.character(mutation_dups$Tumor_Sample_Barcode)
    mutation_info$Gene_Symbol <- entrezIDtoSymbol(mutation_info$Entrez_Gene_Id)
    mutation_info$Mutation_Count <- 1
    for(dup in 1:nrow(mutation_dups)) {
      idx <- intersect(which(mutation_info$Tumor_Sample_Barcode == mutation_dups$Tumor_Sample_Barcode[dup]),
                       which(mutation_info$Entrez_Gene_Id == mutation_dups$Entrez_Gene_Id[dup]))
      mutation_info$Mutation_Count[idx] <- mutation_info$Mutation_Count[idx] + 1
    }
    
    ### save the info to the list
    tcga_maf_per_sample[[tissue_name]] <- mutation_info
    
  }
  
  ### set README function
  README <- function() {
    writeLines(paste(rep("#", 100), collapse = ""))
    writeLines("The \"tcga_maf_per_sample\" is a list that contains mutation information.")
    writeLines("It has 33 entries which indiciates that there are info for 33 TCGA tissues.")
    writeLines("names(tcga_maf_per_sample) will reveal which mutation info is from which tissue.")
    writeLines("In each of the 33 entries, there is a data frame of the mutation information.")
    writeLines("The rows are unique sample-gene pairs that have at least one mutation.")
    writeLines("E.g., a row indicates that in the given sample and in the given gene,")
    writeLines("it has X (specified in the Mutation_Count column) mutations.")
    writeLines("The columns are:")
    writeLines("\tTumor_Sample_Barcode: they are TCGA sample barcodes indicating each of unique TCGA samples")
    writeLines("\tEntrez_Gene_Id: Entrez gene IDs that the mutations were found")
    writeLines("\tGene_Symbol: gene symbols transformed from the Entrez_Gene_Id column")
    writeLines("\tMutation_Count: the number of mutations found in the given sample and in the given gene")
    writeLines("")
    writeLines("It was generated by using TCGAMutInfoPerSample.R.")
    writeLines("The original MAF files were downloaded from GDC database,")
    writeLines("then mutations were counted for every sample-gene pair for each tissue.")
    writeLines(paste(rep("#", 100), collapse = ""))
  }
  
  ### save as RDA
  save(list = c("tcga_maf_per_sample", "README"), file = outputPath)
  
}

