
### Jul 6, 2018
### Hyunjin Kim
### MakeGTEx_cmat.R

### set raw count path
rCntPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/ExpressionMatrices/GTEx_processed/separated_counts/cleaned/"

### get file names from the directory
f <- list.files(rCntPath)

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

### () removal
f2 <- sapply(f, refineFileName, USE.NAMES = FALSE)

### get short names
f3 <- sapply(f2, getShortGTEx, USE.NAMES = FALSE)

### make an empty variable name
cmatNames <- NULL

### iteratively load matrix and assign it as a variable
for(i in 1:length(f)) {
  ### load a matrix
  mat <- read.table(file = paste0(rCntPath, f[i]),
                    header = TRUE,
                    sep ="\t",
                    row.names = 1,
                    check.names = FALSE)
  
  ### remove gene symbol column
  mat <- mat[,-1]
  
  ### only with more than 100 samples
  if(ncol(mat) >= 100) {
    ### save the mat as a variable
    assign(paste("cmat", "gtex", f3[i], sep = "_"), mat, envir = globalenv())
    
    ### add variable name to cmatNames
    cmatNames <- c(cmatNames, paste("cmat", "gtex", f3[i], sep = "_"))
  }
}

### README function
README = function(){
  writeLines("The objects are raw count matrices of 36 GTEx tissues")
  writeLines("The original files are located in below:")
  writeLines("/ifs/archive/shares/af_lab/GTEx/ExpressionMatrices/GTEx_processed/separated_counts/cleaned")
  writeLines("A matrix variable per file is named as follows: cmat_gtex_<interactome_acronym>")
}

### save as RDA file
save(list = c("cmatNames", cmatNames, "README"), file = "./GTEx_36_countmat.rda")

