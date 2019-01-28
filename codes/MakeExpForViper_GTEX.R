###
#   File name : MakeExpForViper_GTEX.R
#   Author    : Hyunjin Kim
#   Date      : Apr 6, 2018
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Load all the non-cleaned counts of all the tissues and normalize them all together
#
#   Instruction
#               1. Source("MakeExpForViper_GTEX.R")
#               2. Run the function "norm_counts()" - specify the input file directories (Aracne networks & expressions) and output path
#               3. GTEx expression matrices (and their combined RDA file) will be generated
#
#   Example
#               > source("The_directory_of_MakeExpForViper_GTEX.R/MakeExpForViper_GTEX.R")
#               > norm_counts(fileNamePath="./results/Aracne/GTEx2/",
#                             expPath="./results/separated_counts/",
#                             outputPath="./GTEx_36_EMat_Viper.rda")
###

norm_counts <- function(fileNamePath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/Aracne2/",
                        expPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/ExpressionMatrices/GTEx_processed/separated_counts/",
                        outputPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/GTEx_36_EMat_Viper.rda") {
  
  ### load library
  if(!require(DESeq2)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("DESeq2")
    library(DESeq2)
  }
  
  
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
  
  
  ### get exp paths correspond to aracne files
  exp_f <- list.files(expPath)
  temp <- strsplit(f, split = "-", fixed = TRUE)
  f2 <- 0
  
  for(i in 1:length(temp)) {
    if(length(temp[[i]]) < 3) {
      f2[i] <- paste0(substr(f[i], 1, nchar(f[i])-10), ".txt")
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
    assign(gtexMatNames[i], d, envir = globalenv())
  }
  
  ### set matNames as global variable
  assign("emat_gtex_names", gtexMatNames, envir = globalenv())
  
  
  ### get number of samples of all the tissues
  sampleNum <- 0
  for(i in 1:length(gtexMatNames)) {
    sampleNum[i] <- ncol(get(gtexMatNames[i]))
  }
  
  ### set sampleNum as global variable
  assign("emat_gtex_sampleNum", sampleNum, envir = globalenv())
  
  
  ### combine all the tissues
  d <- get(gtexMatNames[1])
  for(i in 2:length(gtexMatNames)) {
    d <- cbind(d, get(gtexMatNames[i]))
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
  assign("emat_gtex_all", as.matrix(RNASEQwithVST(d)), envir = globalenv())
  
  
  ### split the combined one into tissue-specific data
  for(i in 1:length(gtexMatNames)) {
    assign(gtexMatNames[i], emat_gtex_all[,(sum(emat_gtex_sampleNum[0:(i-1)])+1):sum(emat_gtex_sampleNum[0:i])], envir = globalenv())
  }
  
  
  ### save the results
  vars <- c(emat_gtex_names, "emat_gtex_names")
  save(list = vars, file = outputPath)
  
}

