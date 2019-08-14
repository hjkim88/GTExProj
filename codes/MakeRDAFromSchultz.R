### There is a dataset of GTEx and TCGA gene expression data
### that the data was preprocessed by the same pipeline and
### batch effect corrected with TCGA normal samples
### The dataset is available from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5903355/
### This code is to combine all the gene expression data into a RDA file

### downloaded data path
dataPath <- "C:/Research/CUMC/GTExProj/data/Schultz_GTEX_TCGA/"

### GTEX - TCGA mapping info path
mapInfoPath <- "./data/RDA_Files/GTEx_TCGA_Map.rda"

### GTEX sample info path
gtexSampleInfoPath <- "./data/GTEx_Data_V6_SampleData.csv"

### get the file names
fileList <- list.files(dataPath)

### remove TCGA normal files from the list
rIdx <- grep(pattern = "tcga.txt", fileList)
if(length(rIdx) > 0) {
  fileList <- fileList[-rIdx]
}

### load map info path
load(mapInfoPath)

### mapping our sample info and the file name
file_tissue <- sapply(fileList, function(x) {
  temp <- strsplit(x, split = "-", fixed = TRUE)[[1]]
  if(grepl("gtex", x)) {
    temp2 <- paste0(toupper(substr(temp[[1]], 1, 1)), substring(temp[[1]], 2))
    if(temp2 == "Esophagus_gas") {
      temp2 <- "EsophagusGasJun"
    } else if(temp2 == "Esophagus_muc") {
      temp2 <- "EsophagusMuc"
    } else if(temp2 == "Esophagus_mus") {
      temp2 <- "EsophagusMus"
    }
    return(temp2)
  } else {
    return(temp[[1]])
  }
  }, USE.NAMES = TRUE)

### fill out the new GTEx_TCGA_Map with the file name notations
GTEx_TCGA_Map <- cbind(GTEx_TCGA_Map, GTEx_File_Name=NA, TCGA_File_Name=NA)
for(i in 1:length(file_tissue)) {
  if(file_tissue[i] == "Colon") {
    idx <- rbind(which(GTEx_TCGA_Map == "ColonSig", arr.ind = TRUE),
                 which(GTEx_TCGA_Map == "ColonTra", arr.ind = TRUE))
  } else {
    idx <- which(GTEx_TCGA_Map == file_tissue[i], arr.ind = TRUE)
  }
  if(nrow(idx) > 0) {
    if(grepl("gtex", names(file_tissue)[i])) {
      GTEx_TCGA_Map[idx[,1],3] <- names(file_tissue)[i]
    } else {
      GTEx_TCGA_Map[idx[,1],4] <- names(file_tissue)[i]
    }
  }
}

### load library
if(!require("org.Hs.eg.db", quietly = TRUE)) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("org.Hs.eg.db")
  require("org.Hs.eg.db", quietly = TRUE)
}

### load the gene expression files and save them as designated objects
unique_files <- unique(c(GTEx_TCGA_Map[,3], GTEx_TCGA_Map[,4]))
unique_files <- unique_files[which(!is.na(unique_files))]
object_names <- sapply(fileList, function(x) {
  temp <- strsplit(x, split = "-", fixed = TRUE)[[1]]
  if(grepl("gtex", x)) {
    temp2 <- paste0(toupper(substr(temp[[1]], 1, 1)), substring(temp[[1]], 2))
    if(temp2 == "Esophagus_gas") {
      temp2 <- "EsophagusGasJun"
    } else if(temp2 == "Esophagus_muc") {
      temp2 <- "EsophagusMuc"
    } else if(temp2 == "Esophagus_mus") {
      temp2 <- "EsophagusMus"
    }
    return(paste0("GTEx_", temp2))
  } else {
    return(paste0("TCGA_", temp[[1]]))
  }
}, USE.NAMES = TRUE)

for(i in 1:length(unique_files)) {
  ### load the file
  temp <- read.table(paste0(dataPath, unique_files[i]), header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  
  ### remove rows that does not have Entrez ID
  temp <- temp[which(temp[,2] > 0),]
  
  ### if there are duplicate Entrez IDs, investigate deeper
  dupIdx <- union(which(duplicated(temp[,2])), which(duplicated(temp[,2], fromLast = TRUE)))
  temp[dupIdx,2] <- mapIds(org.Hs.eg.db, keys = temp[dupIdx,1], column = "ENTREZID", keytype = "SYMBOL")
  
  ### if there are NA for Entrez IDs, use ALIAS
  naIdx <- which(is.na(temp[,2]))
  temp[naIdx,2] <- mapIds(org.Hs.eg.db, keys = temp[naIdx,1], column = "ENTREZID", keytype = "ALIAS")

  ### remove rows that still have NAs or duplicates(retain the first ones)
  temp <- temp[which(!is.na(temp[,2])),]
  temp <- temp[which(!duplicated(temp[,2])),]
  
  ### set row names
  rownames(temp) <- temp[,2]
  
  ### remove unneccessary columns
  temp <- temp[,-c(1,2)]
  
  ### order the rows in accessending order
  temp <- temp[order(as.numeric(rownames(temp))),]
  
  ### set the data as a global object
  assign(paste0("Schultz_Gene_Expression_", object_names[unique_files[i]]), temp, envir = globalenv())
}

### set object names
Schultz_Gene_Expression_Names <- paste0("Schultz_Gene_Expression_", object_names[unique_files])

### load GTEX sample info
sampleInfo <- read.csv(gtexSampleInfoPath, check.names = FALSE, stringsAsFactors = FALSE)
rownames(sampleInfo) <- sampleInfo$SAMPID

### split Colon object into 2 (ColonSig & ColonTra)
sampleInfo <- sampleInfo[which(sampleInfo$SMTS == "Colon"),]
ColonTra <- Schultz_Gene_Expression_GTEx_Colon[,which(sampleInfo[colnames(Schultz_Gene_Expression_GTEx_Colon),"SMTSD"] == "Colon - Transverse")]
ColonSig <- Schultz_Gene_Expression_GTEx_Colon[,which(sampleInfo[colnames(Schultz_Gene_Expression_GTEx_Colon),"SMTSD"] == "Colon - Sigmoid")]
assign("Schultz_Gene_Expression_GTEx_ColonTra", ColonTra, envir = globalenv())
assign("Schultz_Gene_Expression_GTEx_ColonSig", ColonSig, envir = globalenv())
Schultz_Gene_Expression_Names <- Schultz_Gene_Expression_Names[-which(Schultz_Gene_Expression_Names == "Schultz_Gene_Expression_GTEx_Colon")]
Schultz_Gene_Expression_Names <- c(Schultz_Gene_Expression_Names, "Schultz_Gene_Expression_GTEx_ColonTra", "Schultz_Gene_Expression_GTEx_ColonSig")
Schultz_Gene_Expression_Names <- Schultz_Gene_Expression_Names[order(Schultz_Gene_Expression_Names)]

### make a README
README <- function() {
  writeLines(paste(rep("#", 100), collapse = ""))
  writeLines("The RDA file contains normalized gene expression data of GTEx and TCGA.")
  writeLines("The data was preprocessed by the same pipeline and batch effect corrected.")
  writeLines("The dataset is from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5903355/")
  writeLines("Each data frame has gene expression data.")
  writeLines("The rows are Entrez IDs and the columns are sample IDs.")
  writeLines("The \"GTEx_TCGA_Map\" has a comparison info between GTEx and TCGA.")
  writeLines("The thrid and fourth columns of it have the file names of the")
  writeLines("normalized gene expression data.")
  writeLines(paste(rep("#", 100), collapse = ""))
}

### save the results as a RDA file
save(list = c(Schultz_Gene_Expression_Names, "Schultz_Gene_Expression_Names", "GTEx_TCGA_Map", "README"),
     file = "C:/Research/CUMC/GTExProj/data/RDA_Files/All_20_Schultz_Gene_Expressions.rda")
