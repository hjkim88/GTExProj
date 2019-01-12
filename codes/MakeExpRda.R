###
#   File name : MakeExpRda.R
#   Author    : Hyunjin Kim
#   Date      : Feb 16, 2018
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Make one RDA file with GTEx and TCGA expressions
#
#   Instruction
#               1. Source("MakeExpRda.R")
#               2. Run the function "makeRda()" - specify the input (file name & file location) directories and output file path
#               3. One combined expression file in RDA format will be generated in the ouput path
#
#   Example
#               > source("The_directory_of_MakeExpRda.R/MakeExpRda.R")
#               > makeRda(fileNamePath="./results/Aracne/",
#                         fileLocationPath="./results/aracne_ready/",
#                         outputFilePath="./results/aracne_ready/all_64_expmat.rda")
###

makeRda <- function(fileNamePath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/Aracne2/",
                    fileLocationPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/ExpressionMatrices/GTEx_processed/aracne_ready_counts/vst_clean/",
                    outputFilePath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/all_64_expmat.rda") {
  
  ### collect files from the paths
  ### to get file names which were involved in constructing Aracne networks
  gtexFileNames <- list.files(paste0(fileNamePath, "GTEx/"))
  tcgaFileNames <- list.files(paste0(fileNamePath, "TCGA/"))
  
  
  ### shorten the GTEx file names
  ### E.g. "Brain-CerebellarHemisphere_vst" -> "BrainCerHem" 
  getShortGTEx <- function(gtex_file_name) {
    fixed_name <- substr(gtex_file_name, 1, nchar(gtex_file_name)-4)
    fixed_name <- strsplit(fixed_name, "_", fixed = TRUE)[[1]][1]
    
    temp <- strsplit(fixed_name, "-", fixed = TRUE)[[1]]
    
    if(length(temp) > 1) {
      temp2 <- unlist(gregexpr("[A-Z]", temp[2]))
      
      temp3 <- ""
      if((length(temp2) > 1) && (abs(temp2[1] - temp2[2]) > 2)) {
        for(i in 1:2) {
          temp3 <- paste0(temp3, substr(temp[2], temp2[i], temp2[i]+2))
        }  
      } else {
        temp3 <- substr(temp[2], temp2[1], temp2[1]+2)
      }
      
      fixed_name <- paste0(temp[1], temp3)
    } else {
      fixed_name <- temp[1]
    }
    
    return(fixed_name)
  }
  
  
  ### set matNames
  gtexMatNames <- 0
  for(i in 1:length(gtexFileNames)) {
    gtexMatNames[i] <- paste0("emat_gtex_", getShortGTEx(gtexFileNames[i]))
  }
  tcgaMatNames <- paste0("emat_", tcgaFileNames)
  
  ### 1:36 GTExs, 37:64 TCGAs
  matNames <- c(gtexMatNames, tcgaMatNames)
  
  
  ### load exp datasets
  for(i in 1:length(gtexFileNames)) {
    assign(gtexMatNames[i], as.matrix(read.table(paste0(fileLocationPath, "GTEx/", gtexFileNames[i], ".dat"),
                                                 row.names = 1, header = TRUE, sep = "\t", quote = "",
                                                 check.names = FALSE, stringsAsFactors = FALSE), envir = globalenv()))
  }
  for(i in 1:length(tcgaFileNames)) {
    assign(tcgaMatNames[i], as.matrix(read.table(paste0(fileLocationPath, "TCGA/", tcgaFileNames[i], ".dat"),
                                                 row.names = 1, header = TRUE, sep = "\t", quote = "",
                                                 check.names = FALSE, stringsAsFactors = FALSE), envir = globalenv()))
  }
  
  
  ### README function
  README = function(){
    writeLines("- Object names start with \"emat_\" contains normalized gene expressions of the corresponding tissue")
    writeLines("  If a name contains \"gtex\", then it is a gtex data, and \"tcga\", then it is tcga data")
    writeLines("  row names indicate Entrez IDs and column names are sample names")
    writeLines("")
    writeLines("- The variable \"matNames\" has all the exp variable names")
    writeLines("  The 1:36 are GTEx matrix names, and the 37:64 are TCGA matrix names")
  }
  
  
  ### organize all the variables and save them into RDA file
  vars = c(matNames, "matNames", "README")
  save(list = vars, file = outputFilePath)
  
}

