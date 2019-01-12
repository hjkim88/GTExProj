###
#   File name : Expmat_GTEX.R
#   Author    : Hyunjin Kim
#   Date      : Apr 4, 2018
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : load expressions and make emat files for GTEx datasets
#
#   Instruction
#               1. Source("Expmat_GTEX.R")
#               2. Run the function "makeReady()" - specify the input (Aracne networks & expressions) directory and output path
#               3. GTEx expression matrices (and their combined RDA file) will be generated
#
#   Example
#               > source("The_directory_of_Expmat_GTEX.R/Expmat_GTEX.R")
#               > ematReady(fileNamePath="./results/Aracne/GTEx2/",
#                           expPath="./results/transformed_counts/vst_clean/",
#                           outputPath="./results/transformed_counts/vst_clean/gtex_36_emat.rda")
###

ematReady <- function(fileNamePath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/Aracne2/",
                      expPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/ExpressionMatrices/GTEx_processed/transformed_counts/vst_clean/",
                      outputPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/gtex_36_emat.rda") {
  
  ### collect files from the fileNamePath
  f <- list.files(fileNamePath)
  f <- f[which(endsWith(f, "vst"))]
  
  
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
  
  
  exp_f <- list.files(expPath)
  temp <- strsplit(f, split = "-", fixed = TRUE)
  f2 <- 0
  
  for(i in 1:length(temp)) {
    if(length(temp[[i]]) < 3) {
      f2[i] <- paste0(f[i], ".txt")  
    } else {
      f2[i] <- exp_f[which(startsWith(exp_f, paste0(temp[[i]][1], "-", temp[[i]][2])))]
    }
  }
  
  
  ### set matNames
  gtexMatNames <- 0
  for(i in 1:length(f)) {
    gtexMatNames[i] <- paste0("emat_gtex_", getShortGTEx(f[i]))
  }
  
  
  ### load exp datasets
  for(i in 1:length(f)) {
    d <- read.table(paste0(expPath, f2[i]),
                    row.names = 1, header = TRUE, sep = "\t",
                    check.names = FALSE, stringsAsFactors = FALSE)
    d <- d[,-1]
    assign(as.matrix(gtexMatNames[i]), d, envir = globalenv())
  }
  
  assign("emat_gtex_names", gtexMatNames, envir = globalenv())
  
  ### save the results
  vars <- c(gtexMatNames, "gtexMatNames")
  save(list = vars, file = outputPath)
  
}
