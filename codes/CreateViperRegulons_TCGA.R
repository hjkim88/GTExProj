###
#   File name : CreateViperRegulons_TCGA.R
#   Author    : Hyunjin Kim
#   Date      : Apr 4, 2017
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Create Viper regulon objects from 4-col format Aracne files
#
#   Instruction
#               1. Source("CreateViperRegulons_TCGA.R")
#               2. Run the function "createRegulons()" - specify the input dirs (expression and 4-col files) and output directory
#               3. The RDA file of Viper regulons will be generated in the output directory
#
#   Example
#               > source("The_directory_of_CreateViperRegulons_TCGA.R/CreateViperRegulons_TCGA.R")
#               > createRegulons(expDir="./results/aracne_ready/GTEx/vst_clean/", networkDir="./results/Aracne/GTEx2/", outputPath="./")
###


createRegulons <- function(expDir="./results/aracne_ready/TCGA/", networkDir="./results/Aracne/TCGA/", outputPath="./TCGA_28_Regulons.rda") {
  
  ### load library
  if(!require(viper)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("viper")
    library(viper)
  }
  
  
  ### collect files from the cntPath
  f <- list.files(networkDir)
  
  
  ### set regNames
  tcgaRegNames <- 0
  for(i in 1:length(f)) {
    tcgaRegNames[i] <- paste0("regulon_", f[i])
  }
  
  
  ### a function to generate regulon object
  consolidateNets <- function(expmat.file, netFiles, resultDir) {
    
    # Turn the 4-col files network files into one merged 3-col file, as needed by the method
    # aracne2regulon, and store it in a temporary directory
    temp_dir = paste(resultDir, "temp_", as.integer(Sys.time()), "/", sep="")
    dir.create(temp_dir)
    temp_combined = paste(temp_dir, "consolidatedNet.tsv", sep="") # Put all interacions here
    master_mat = NULL
    for (i in 1:length(netFiles)){
      mat = as.matrix(read.table(netFiles[i]))
      master_mat = rbind(master_mat, mat)
    }
    
    write.table(master_mat[, 1:3], file = temp_combined, row.names = FALSE, col.names = FALSE, sep="\t")
    
    
    # Run the method that will create the 'regulon' object
    regulon <- aracne2regulon(temp_combined, expmat.file, format="3col")
    
    # order the regulon targets in increasing order of their Entrez id
    for (i in 1:length(regulon)){
      x = as.integer(names(regulon[[i]]$tfmode))
      regulon[[i]]$tfmode = regulon[[i]]$tfmode[order(x)]
      regulon[[i]]$likelihood = regulon[[i]]$likelihood[order(x)]
    }
    regulon = regulon[order(as.integer(names(regulon)))]
    
    # Cleanup
    unlink(paste0(dirname(temp_dir), "/", basename(temp_dir)), recursive=TRUE)
    
    return(regulon)
  }
  
  
  ### iteratively perform conslidateNets() on the datasets
  for(i in 1:length(f)) {
    
    ### set necessary input values
    gexPath <- paste0(expDir, f[i], ".dat")
    netPath <- c(paste0(networkDir, f[i], "/", f[i], "_tf2.tsv"), paste0(networkDir, f[i], "/", f[i], "_cotf2.tsv"), paste0(networkDir, f[i], "/", f[i], "_signal2.tsv"))
    resultPath <- paste0(networkDir, f[i], "/")
    
    ### consolidate the networks
    assign(tcgaRegNames[i], consolidateNets(gexPath, netPath, resultPath), envir = globalenv())
    
  }
  
  assign("tcgaRegNames", tcgaRegNames, envir = globalenv())
  
  ### save the results
  vars <- c(tcgaRegNames, "tcgaRegNames")
  save(list = vars, file = outputPath)
  
}
