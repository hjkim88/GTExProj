###
#   File name : PreprocessMAFs.R
#   Author    : Hyunjin Kim
#   Date      : Aug 30, 2019
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : There are TCGA MAF files from various tissues and from various pipelines.
#               We are trying to preprocess them to make one MAF file for one tissue.
#               For each tissue, 4 different MAF files are combined and with duplicates removed.
#
#   Instruction
#               1. Source("PreprocessMAFs.R")
#               2. Run the function "preprocessMAFs" - specify the necessary input directory and output directory
#               3. The MAF files will be generated in the output directory
#
#   Example
#               > source("The_directory_of_PreprocessMAFs.R/PreprocessMAFs.R")
#               > preprocessMAFs(rawFileDir="C:/Research/CUMC/GTExProj/data/TCGA/MAF/raw/",
#                                outputDir="C:/Research/CUMC/GTExProj/data/TCGA/MAF/preprocessed/")
###

preprocessMAFs <- function(rawFileDir="C:/Research/CUMC/GTExProj/data/TCGA/MAF/raw/",
                           outputDir="C:/Research/CUMC/GTExProj/data/TCGA/MAF/preprocessed/") {
  
  ### get maf.gz file list
  f <- list.files(rawFileDir, pattern = "*.maf.gz$")
  
  ### get unique tissues
  unique_tissues <- unique(sapply(f, function(x) strsplit(x, ".", TRUE)[[1]][2]))
  
  ### for each tissue combine MAF files
  for(tissue in unique_tissues) {
    
    ### get files for a given tissue
    f2 <- f[which(startsWith(f, paste0("TCGA.", tissue)))]
    
    ### load MAFs
    for(i in 1:length(f2)) {
      if(i == 1) {
        maf <- read.table(file = paste0(rawFileDir, f2[i]), header = TRUE, sep = "\t", quote = "",
                          stringsAsFactors = FALSE, check.names = FALSE)
      } else {
        maf <- rbind(maf, read.table(file = paste0(rawFileDir, f2[i]), header = TRUE, sep = "\t", quote = "",
                                     stringsAsFactors = FALSE, check.names = FALSE))
      }
    }
    
    ### remove duplicates in the combined MAF file
    duplicates <- duplicated(maf[,c("Hugo_Symbol", "Entrez_Gene_Id", "Chromosome",
                                    "Start_Position", "End_Position", "Strand", "Reference_Allele",
                                    "Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode")])
    maf <- maf[!duplicates,]
    
    ### these are optional comments at the top of the MAF file
    comments <- "#version gdc-1.0.0\n"
    comments <- paste0(comments, "#filedate 20190903\n")
    comments <- paste0(comments, "#annotation.spec gdc-1.0.1-public\n")
    comments <- paste0(comments, paste("#n.analyzed.samples", length(unique(maf$Tumor_Sample_Barcode))), "\n")
    comments <- paste0(comments, paste0("#tumor.aliquots.submitter_id ", paste(unique(maf$Tumor_Sample_Barcode), collapse = ",")))
    
    ### write out the MAF file
    write.table(comments, file = paste0(outputDir, "TCGA_", tissue, ".maf"),
                sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(maf, file = paste0(outputDir, "TCGA_", tissue, ".maf"),
                sep = "\t", row.names = FALSE, quote = FALSE, append = TRUE)
    
    ### garbage collection
    gc()
    
  }

}
