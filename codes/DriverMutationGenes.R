###
#   File name : DriverMutationGenes.R
#   Author    : Hyunjin Kim
#   Date      : Sep 3, 2019
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : There are preprocessed TCGA MAF files. Use maftools package to
#               compute driver mutation genes for each TCGA tissue.
#
#   Instruction
#               1. Source("DriverMutationGenes.R")
#               2. Run the function "preprocessMAFs" - specify the necessary input directory and output directory
#               3. The result files will be generated in the output directory
#
#   Example
#               > source("The_directory_of_DriverMutationGenes.R/DriverMutationGenes.R")
#               > compute_dm_genes(mafFileDir="C:/Research/CUMC/GTExProj/data/TCGA/MAF/preprocessed/",
#                                  outputDir="C:/Research/CUMC/GTExProj/data/TCGA/MAF/driver_mutation_genes/")
###

compute_dm_genes <- function(mafFileDir="C:/Research/CUMC/GTExProj/data/TCGA/MAF/preprocessed/",
                             outputDir="C:/Research/CUMC/GTExProj/data/TCGA/MAF/driver_mutation_genes/") {
  
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
  
  ### get maf file list
  f <- list.files(mafFileDir)
  
  ### create am empty list for RDA
  tcga_driver_mutation_genes <- vector("list", length = length(f))
  
  ### for each TCGA tissue file, compute driver mutation genes
  for(i in 1:length(f)) {
    
    ### get maf file name
    maf_file <- f[i]
    
    ### load total MAF file
    maf <- read.maf(maf = paste0(mafFileDir, maf_file))
    
    ### create a result directory
    resultDir <- paste0(outputDir, substr(maf_file, 1, nchar(maf_file)-4), "/")
    dir.create(resultDir)
    
    ### create a summary plot
    png(paste0(resultDir, "Summary_Plot_", substr(maf_file, 1, nchar(maf_file)-4), ".png"), width = 2000, height = 1000, res = 200)
    plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = "median", dashboard = TRUE, titvRaw = FALSE)
    dev.off()
    
    ### get driver genes
    driver_genes <- oncodrive(maf = maf, AACol = "HGVSp_Short", minMut = 5)
    
    ### save the result
    tcga_driver_mutation_genes[[i]] <- driver_genes
    names(tcga_driver_mutation_genes)[i] <- substr(maf_file, 1, nchar(maf_file)-4)
    
    ### write out the driver gene result
    write.xlsx2(driver_genes, file = paste0(resultDir, "Driver_Genes_Table_", substr(maf_file, 1, nchar(maf_file)-4), ".xlsx"),
                sheetName = "Driver_Genes", row.names = FALSE)
    
    ### garbage collection
    gc()
    
  }
  
  ### set README function
  README <- function(){
    writeLines(paste(rep("#", 100), collapse = ""))
    writeLines("The \"tcga_driver_mutation_genes\" is a list that contains driver mutation gene info")
    writeLines("There are 33 items in the list and in each list, it has a data frame for the info for the given TCGA tissue")
    writeLines("names(tcga_driver_mutation_genes) shows TCGA tissue names for each list")
    writeLines("Each data frame is a result of OncodriveCLUST method.")
    writeLines("The OncodriveCLUST method is used for the analysis. The concept is based on the fact that")
    writeLines("most of the variants in cancer causing genes are enriched at few specific loci (aka hotspots).")
    writeLines("This method takes advantage of such positions to identify cancer driver genes.")
    writeLines("-")
    writeLines("Hugo_Symbol: Name of driver gene")
    writeLines("Frame_Shift_Del: TMB based on frame shift deletion")
    writeLines("Frame_Shift_Ins: TMB based on frame shift insertion")
    writeLines("In_Frame_Del: TMB based on in frame deletion")
    writeLines("In_Frame_Ins: TMB based on in frame insertion")
    writeLines("Missense_Mutation: TMB based on missense mutations")
    writeLines("Nonsense_Mutation: TMB based on nonsense mutations")
    writeLines("Nonstop_Mutation: TMB based on nonstop mutations")
    writeLines("Splice_Site: TMB based on splice sites")
    writeLines("Translation_Start_Site: TMB based on translation start sites")
    writeLines("total: TMB based on all types of muations")
    writeLines("MutatedSamples: the number of samples that have mutations in the gene")
    writeLines("AlteredSamples: the number of samples that have proteins altered in the gene")
    writeLines("clusters: the number of clusters found in the gene")
    writeLines("muts_in_clusters: the number of mutations in the clusters")
    writeLines("clusterScores: This is a score that represents how well the mutations are clustered")
    writeLines("protLen: total length of the altered proteins in the gene")
    writeLines("zscore: z-score of the clustering scores")
    writeLines("pval: p-value of the zscore")
    writeLines("fdr: adjusted p-value of the p-value – false discovery rate")
    writeLines("fract_muts_in_clusters: the number of fraction of mutations observed in the clusters")
    writeLines(paste(rep("#", 100), collapse = ""))
  }
  
  ### save as RDA
  save(list = c("tcga_driver_mutation_genes", "README"), file = paste0(outputDir, "TCGA_33_Driver_Mutation_Genes.rda"))
  
}
