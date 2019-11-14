###
#   File name : TopMutatedRegulons.R
#   Author    : Hyunjin Kim
#   Date      : Nov 4, 2019
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : There are preprocessed TCGA MAF files. Use maftools package to
#               compute top mutated regulons for each TCGA tissue.
#
#   Instruction
#               1. Source("TopMutatedRegulons.R")
#               2. Run the function "mutated_genes" - specify the necessary input directory and output directory
#               3. The result RDA file will be generated in the output path
#
#   Example
#               > source("The_directory_of_TopMutatedRegulons.R/TopMutatedRegulons.R")
#               > mutated_genes(mafFileDir="C:/Research/CUMC/GTExProj/data/TCGA/MAF/preprocessed/",
#                               aracneRDAPath="C:/Research/cumc/GTExProj/data/RDA_Files/All_62_ARACNE.rda",
#                               outputPath="C:/Research/CUMC/GTExProj/data/RDA_Files/TCGA_26_Top_Mutated_Regulons.rda")
###

mutated_genes <- function(mafFileDir="C:/Research/CUMC/GTExProj/data/TCGA/MAF/preprocessed/",
                          aracneRDAPath="C:/Research/cumc/GTExProj/data/RDA_Files/All_62_ARACNE.rda",
                          outputPath="C:/Research/CUMC/GTExProj/data/RDA_Files/TCGA_26_Top_Mutated_Regulons.rda") {
  
  ### load libraries
  if(!require(maftools, quietly = TRUE)) {
    if(!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("maftools", version = "3.8")
    require(maftools, quietly = TRUE)
  }
  if(!require(metaseqR, quietly = TRUE)) {
    if(!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("metaseqR", version = "3.8")
    require(metaseqR, quietly = TRUE)
  }
  
  ### load data
  load(aracneRDAPath)
  
  ### get maf file list
  f <- list.files(mafFileDir)
  
  ### get TCGA Aracne names
  tcga_aracne_names <- varNames[grep("tcga", varNames)]
  
  ### create am empty list for RDA
  tcga_highly_mutated_regulons <- vector("list", length = length(tcga_aracne_names))
  names(tcga_highly_mutated_regulons) <- tcga_aracne_names
  
  ### create temporary directory
  temp_path <- paste0("./temp", sample(10000, 1), "/")
  dir.create(path = temp_path)
  
  ### for each TCGA tissue file, identify mutated regulons
  for(aracne_name in tcga_aracne_names) {
    
    ### get Aracne network
    aracne_net <- get(aracne_name)
    
    ### get maf file name
    maf_file <- paste0(toupper(aracne_name), ".maf")
    
    ### load total MAF file
    maf <- read.maf(maf = paste0(mafFileDir, maf_file))
    
    ### make sample-gene mutation info
    mutation_info <- data.frame(maf@data[which(!duplicated(maf@data[,c("Entrez_Gene_Id", "Tumor_Sample_Barcode")])),
                                         c("Tumor_Sample_Barcode", "Entrez_Gene_Id")],
                                stringsAsFactors = FALSE, check.names = FALSE)
    mutation_info$Tumor_Sample_Barcode <- as.character(mutation_info$Tumor_Sample_Barcode)
    mutation_dups <- data.frame(maf@data[which(duplicated(maf@data[,c("Entrez_Gene_Id", "Tumor_Sample_Barcode")])),
                                         c("Tumor_Sample_Barcode", "Entrez_Gene_Id")],
                                stringsAsFactors = FALSE, check.names = FALSE)
    mutation_dups$Tumor_Sample_Barcode <- as.character(mutation_dups$Tumor_Sample_Barcode)
    mutation_info$Count <- 1
    for(dup in 1:nrow(mutation_dups)) {
      idx <- intersect(which(mutation_info$Tumor_Sample_Barcode == mutation_dups$Tumor_Sample_Barcode[dup]),
                       which(mutation_info$Entrez_Gene_Id == mutation_dups$Entrez_Gene_Id[dup]))
      mutation_info$Count[idx] <- mutation_info$Count[idx] + 1
    }
    
    ### get unique samples
    unique_samples <- unique(mutation_info$Tumor_Sample_Barcode)
    
    ### make an empty list
    tcga_highly_mutated_regulons[[aracne_name]] <- rep(NA, nrow(aracne_net[[1]]))
    names(tcga_highly_mutated_regulons[[aracne_name]]) <- rownames(aracne_net[[1]])
    
    ### make an empty data frame for integration
    pvs.df <- matrix(NA, nrow(aracne_net[[1]]), length(unique_samples))
    rownames(pvs.df) <- rownames(aracne_net[[1]])
    colnames(pvs.df) <- unique_samples
    
    ### for every hub, calculate Cosmic enrichment p-value
    for(hub in rownames(aracne_net[[1]])) {
      ### get total genes in the network
      interactome_total_genes <- as.character(getInteractomeGenes(aracne_name, count = FALSE))
      
      ### get target genes
      target_genes <- rownames(aracne_net[[2]][[hub]])
      
      ### for every sample, calculate individual enrichment p-value
      for(sample in unique_samples) {
        ### get mutation info for a given sample
        mutInfo <- mutation_info[which(mutation_info$Tumor_Sample_Barcode == sample),]
        
        ### compute enriched genes
        enriched_genes <- intersect(target_genes, mutInfo$Entrez_Gene_Id)
        
        ### compute FET p-value
        ### Fisher's Exact Test
        ###
        ###                   regulon   no-regulon
        ###                 -----------------------
        ### mutated gene    |   X           Y
        ### non-mutated gene|   Z           W
        x <- length(enriched_genes)
        y <- length(mutInfo$Entrez_Gene_Id) - x
        z <- length(target_genes) - x
        w <- length(interactome_total_genes) - x - y - z
        pvs.df[hub,sample] <- fisher.test(matrix(c(x, z, y, w), 2, 2), alternative = "greater")$p.value
      }
    }
    
    ### integrate p-values across all the samples using Fisher's method
    tcga_highly_mutated_regulons[[aracne_name]] <- fisher.method(pvs.df, p.corr = "BH")
    tcga_highly_mutated_regulons[[aracne_name]] <- tcga_highly_mutated_regulons[[aracne_name]][order(-tcga_highly_mutated_regulons[[aracne_name]][,"S"]),]
    
    ### save the temp file
    temp <- tcga_highly_mutated_regulons[[aracne_name]]
    save(list = c("temp"), file = paste0(temp_path, aracne_name, ".rda"))
    
    ### garbage collection
    gc()
    
  }
  
  ### set README function
  README <- function(){
    writeLines(paste(rep("#", 100), collapse = ""))
    writeLines("The \"tcga_highly_mutated_regulons\" is a list that contains highly mutated regulon info")
    writeLines("It has 26 entries which indiciates that there are mutated regulon info for 26 TCGA tissues.")
    writeLines("In each of the 26 entries, there is a data frame of mutated regulons.")
    writeLines("The rows are hubs and the columns are results from Fisher's method of each hub.")
    writeLines("")
    writeLines("S: chi-square statistics of the integrated p-values")
    writeLines("num.p: the number of p-values that were integrated")
    writeLines("p.value: integrated p-values using Fisher's method")
    writeLines("p.adj: corrected p-values using Benjamini-Hochberg")
    writeLines("")
    writeLines("It was generated as the following:")
    writeLines("For every TCGA tissue:")
    writeLines("\tFor every hub H:")
    writeLines("\t\tFor every sample S:")
    writeLines("\t\t\tUse the MAF file to identify the genes in S that are mutated.")
    writeLines("\t\t\tUse FET to assess the enrichment of the regulon of H in mutated genes. Let p(H, S) be the p-value of the FET.")
    writeLines("\t\tIntegrate the p-values p(H, S) from all samples S into a master p-value p(H), using Fisher’s method (https://en.wikipedia.org/wiki/Fisher%27s_method)")
    writeLines("\tCorrect the p-values with BH method.")
    writeLines(paste(rep("#", 100), collapse = ""))
  }
  
  ### save as RDA
  save(list = c("tcga_highly_mutated_regulons", "README"), file = outputPath)
  
}
