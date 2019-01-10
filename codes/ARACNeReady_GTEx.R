###
#   File name : ARACNeReady_GTEx.R
#   Author    : Hyunjin Kim
#   Date      : Mar 1, 2018
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Make Aracne-ready file from the transformed counts from GTEx
#               Does not create files for tissues that have less than 100 samples
#               If the number of sample exceeds 200, 200 samples will be selected based on RIN
#
#   Instruction
#               1. Source("ARACNeReady_GTEx.R")
#               2. Run the function "makeReady()" - specify the input file (transformed counts) directory and output directory
#               3. ARACNe-ready data will be generated in the output directory
#
#   Example
#               > source("The_directory_of_ARACNeReady_GTEx.R/ARACNeReady_GTEx.R")
#               > makeReady(inputPath="./results/transformed_counts/vst_clean/",
#                           sampleInfoPath="./data/GTEx_Data_V6_SampleData.csv",
#                           outputPath="./results/aracne_ready/GTEx/vst_clean/",
#                           rdaFilePath="./data/RDA_Files/")
###

makeReady <- function(inputPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/ExpressionMatrices/GTEx_processed/transformed_counts/vst_clean/",
                      sampleInfoPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/Annotations/GTEx_Data_V6_SampleData.csv",
                      outputPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/ExpressionMatrices/GTEx_processed/aracne_ready_counts/vst_clean/",
                      rdaFilePath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/") {
  
  ### load sample info
  sampleInfo <- read.csv(sampleInfoPath, check.names = FALSE, stringsAsFactors = FALSE)
  rownames(sampleInfo) <- sampleInfo$SAMPID
  
  
  ### make RIN info
  rin <- sampleInfo$SMRIN
  names(rin) <- sampleInfo$SAMPID
  
  
  ### collect files from the cntPath
  f <- list.files(inputPath)
  
  
  ### A function to change "(" or ")" to "-"
  refineFileName <- function(str_line) {
    result_line <- gsub("\\(", "-", str_line)
    result_line <- gsub("\\)", "-", result_line)
    
    return(result_line)
  }
  
  
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
  
  
  ### initialize the count matrix names for a RDA file
  gtex_expmat_names <- NULL
  
  
  ### iteratively perform transforming
  for(i in 1:length(f)) {
    ### load filtered counts
    cnt <- read.table(paste0(inputPath, f[i]), sep="\t", row.names = 1, header = TRUE, check.names = FALSE)
    
    if((ncol(cnt)-1) > 200) {
      ### get RIN for the samples
      temp <- rin[colnames(cnt)[-1]]
      temp <- temp[order(-temp)]
      
      ### only keep top 200 samples based on RIN
      cnt <- cbind(cnt$Gene_Symbol, cnt[,names(temp)[1:200]])
    }
    
    if((ncol(cnt)-1) >= 100) {
      cnt <- cbind(Gene=rownames(cnt), cnt[,-1])
      
      ### save the transformed data
      write.table(cnt, paste0(outputPath, refineFileName(substr(f[i], 1, nchar(f[i])-4)), ".dat"), sep = "\t", row.names = FALSE, quote = FALSE)
      
      ### save the Aracne-ready expressions to a variable
      gtex_expmat_names <- c(gtex_expmat_names, paste0("expmat_gtex_", getShortGTEx(refineFileName(substr(f[i], 1, nchar(f[i])-4)))))
      assign(gtex_expmat_names[length(gtex_expmat_names)], cnt[,-1], envir = globalenv())
    }
  }
  
  
  ### only keep the sample info of Aracne-ready exps
  gtex_sample_info <- sampleInfo[unlist(sapply(gtex_expmat_names, function(x) {
    return(colnames(get(x)))
  })),]
  
  
  ### set README function
  README <- function(){
    writeLines(paste(rep("#", 100), collapse = ""))
    writeLines("There are 36 Aracne-ready expression matrices of GTEx data")
    writeLines("The \"gtex_expmat_names\" is a vector of length 36 and has all the variable names of the matrices")
    writeLines("The \"expmat_gtex_[TISSUE NAME]\" has the Aracne-ready expressions for the given tissue")
    writeLines("The process to make Aracne-ready files are like below:")
    writeLines("1. If there are duplicated gene symbols, then just add all counts up")
    writeLines("2. Remove genes that have NA, 0 or 1 across all samples")
    writeLines("3. Perform VST(Variance Stabilizing Transformation) on counts")
    writeLines("4. Only use tissues with >= 100 samples and if there are more than 200 samples, choose 200 with highest RIN")
    writeLines("The \"gtex_sample_info\" has all the information related the GTEx samples")
    writeLines("The \"gtex_sample_info\" has 5937 rows and 64 columns")
    writeLines("The rows are samples corresponding to the columns of the \"expmat_gtex_[TISSUE NAME]\"")
    writeLines("The columns represent various attributes of the samples")
    writeLines(paste(rep("#", 100), collapse = ""))
  }
  
  
  ### save the raw counts and the sample info as a RDA file
  save(list = c(gtex_expmat_names, "gtex_expmat_names", "gtex_sample_info", "README"),
       file = paste0(rdaFilePath, "/GTEX_36_ARACNE_READY_EXPMAT.rda"))
  
}
