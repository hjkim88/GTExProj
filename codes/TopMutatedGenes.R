###
#   File name : TopMutatedGenes.R
#   Author    : Hyunjin Kim
#   Date      : Oct 28, 2019
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : There are preprocessed TCGA MAF files. Use maftools package to
#               compute top mutated genes for each TCGA tissue.
#
#   Instruction
#               1. Source("TopMutatedGenes.R")
#               2. Run the function "mutated_genes" - specify the necessary input directory and output directory
#               3. The result RDA file will be generated in the output path
#
#   Example
#               > source("The_directory_of_TopMutatedGenes.R/TopMutatedGenes.R")
#               > mutated_genes(mafFileDir="C:/Research/CUMC/GTExProj/data/TCGA/MAF/preprocessed/",
#                               geneLengthPath="C:/Research/CUMC/GTExProj/data/RDA_Files/Gene_Lengths.rda",
#                               outputPath="C:/Research/CUMC/GTExProj/data/RDA_Files/TCGA_33_Top_Mutated_Genes.rda")
###

mutated_genes <- function(mafFileDir="C:/Research/CUMC/GTExProj/data/TCGA/MAF/preprocessed/",
                          geneLengthPath="C:/Research/CUMC/GTExProj/data/RDA_Files/Gene_Lengths.rda",
                          outputPath="C:/Research/CUMC/GTExProj/data/RDA_Files/TCGA_33_Top_Mutated_Genes.rda") {
  
  ### load libraries
  options(java.parameters = "-Xmx8000m")
  if(!require(xlsx, quietly = TRUE)) {
    install.packages("xlsx")
    require(xlsx, quietly = TRUE)
  }
  if(!require(maftools, quietly = TRUE)) {
    if(!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("maftools", version = "3.8")
    require(maftools, quietly = TRUE)
  }
  
  ### load gene length info
  load(geneLengthPath)
  
  ### get maf file list
  f <- list.files(mafFileDir)
  
  ### create am empty list for RDA
  tcga_highly_mutated_genes <- vector("list", length = length(f))
  
  ### for each TCGA tissue file, compute driver mutation genes
  for(i in 1:length(f)) {
    
    ### get maf file name
    maf_file <- f[i]
    
    ### load total MAF file
    maf <- read.maf(maf = paste0(mafFileDir, maf_file))
    
    ### mutation count info
    mut_cnt <- data.frame(maf@gene.summary, stringsAsFactors = FALSE, check.names = FALSE)
    
    ### compute gene length for all the genes
    symbol_entrez_map <- maf@data[which(!duplicated(maf@data$Hugo_Symbol)),c("Hugo_Symbol", "Entrez_Gene_Id")]
    symbol_entrez_map <- data.frame(symbol_entrez_map, stringsAsFactors = FALSE, check.names = FALSE)
    rownames(symbol_entrez_map) <- symbol_entrez_map$Hugo_Symbol
    mut_cnt$Gene_Length <- gene_lengths[as.character(symbol_entrez_map[mut_cnt$Hugo_Symbol,"Entrez_Gene_Id"])]
    
    ### remove genes that do not have gene length
    mut_cnt <- mut_cnt[!is.na(mut_cnt$Gene_Length),]
    
    ### normalize mutation counts by the gene lengths
    mut_cnt$Normalized_Gene_Count <- mut_cnt$total / mut_cnt$Gene_Length
    
    ### order rows based on normalized mutation counts
    mut_cnt <- mut_cnt[order(-mut_cnt$Normalized_Gene_Count),]
    
    ### get highly mutated genes
    tcga_highly_mutated_genes[[i]] <- mut_cnt$Normalized_Gene_Count
    names(tcga_highly_mutated_genes[[i]]) <- mut_cnt$Hugo_Symbol
    names(tcga_highly_mutated_genes)[i] <- substr(maf_file, 1, nchar(maf_file)-4)
    
    ### garbage collection
    gc()
    
  }
  
  ### set README function
  README <- function(){
    writeLines(paste(rep("#", 100), collapse = ""))
    writeLines("The \"tcga_highly_mutated_genes\" is a list that contains highly mutated gene info")
    writeLines("It has 33 entries which indiciates that there are gene list for 33 TCGA tissues.")
    writeLines("In each of the 33 entries, there is a vector of normalized mutation counts")
    writeLines("that were computed by mutation counts divided by the gene length.")
    writeLines("It is ordered in decending order and names(tcga_highly_mutated_genes[[i]]) will")
    writeLines("retreive top mutated genes of the i-th tissue.")
    writeLines(paste(rep("#", 100), collapse = ""))
  }
  
  ### save as RDA
  save(list = c("tcga_highly_mutated_genes", "README"), file = outputPath)
  
}
