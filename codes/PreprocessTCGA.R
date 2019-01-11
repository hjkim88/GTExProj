###
#   File name : PreprocessTCGA.R
#   Author    : Hyunjin Kim
#   Date      : Dec 19, 2018
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : We downloaded raw counts of all the TCGA tissues and now we would like to organize them
#               in one RDA file with sample information including RIN
#
#   * GDC.h38 GENCODE v22 GTF (used in RNA-Seq alignment and by HTSeq) - Ensembl 79
#
#   Instruction
#               1. Source("PreprocessTCGA.R")
#               2. Run the function "preprocess_tcga()" - specify the target directory and output directory
#               3. The preprocessed RDA file will be generated in the output directory
#
#   Example
#               > source("The_directory_of_PreprocessTCGA.R/PreprocessTCGA.R")
#               > preprocess_tcga(targetDir="E:/TCGA/",
#                                 targetSampleInfoPath="E:/TCGA/sample_sheet_11093.xlsx",
#                                 targetMetadataPath="E:/TCGA/metadata_11093.json",
#                                 targetClinicalInfoPath="E:/TCGA/clinical_10237.tsv",
#                                 targetExposureInfoPath="E:/TCGA/exposure_10237.tsv",
#                                 rinInfoDir="//isilon.c2b2.columbia.edu/ifs/archive/TCGA/Open_Access/",
#                                 outputDir="./data/RDA_Files/")
###

preprocess_tcga <- function(targetDir="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/TCGA/ExpressionMatrices/tcga_af_lab_raw_data/",
                            targetSampleInfoPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/TCGA/ExpressionMatrices/tcga_af_lab_raw_data/sample_sheet_11093.xlsx",
                            targetMetadataPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/TCGA/ExpressionMatrices/tcga_af_lab_raw_data/metadata_11093.json",
                            targetClinicalInfoPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/TCGA/ExpressionMatrices/tcga_af_lab_raw_data/clinical_10237.tsv",
                            targetExposureInfoPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/TCGA/ExpressionMatrices/tcga_af_lab_raw_data/exposure_10237.tsv",
                            rinInfoDir="//isilon.c2b2.columbia.edu/ifs/archive/TCGA/Open_Access/",
                            outputDir="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/") {
  
  ### JAVA memory setup - because of xlsx package
  options(java.parameters = "-Xmx2048m")
  
  
  ### load libraries
  if(!require(xlsx, quietly = TRUE)) {
    install.packages("xlsx", dependencies = TRUE, quiet = TRUE)
    require(xlsx, quietly = TRUE)
  }
  if(!require(jsonlite, quietly = TRUE)) {
    install.packages("jsonlite", dependencies = TRUE, quiet = TRUE)
    require(jsonlite, quietly = TRUE)
  }
  
  
  ### load TCGA sample files
  tcga_files <- list.files(targetDir, full.names = TRUE, recursive = TRUE, pattern = "*.gz")
  
  
  ### load TCGA sample info
  tcga_sample_info <- read.xlsx2(file = targetSampleInfoPath, sheetIndex = 1,
                                 stringsAsFactors = FALSE, check.names = FALSE)
  
  
  ### order the sample info by the sample file names
  rownames(tcga_sample_info) <- tcga_sample_info$`File Name`
  tcga_sample_info <- tcga_sample_info[basename(tcga_files),]
  
  
  ### load TCGA metadata info
  metadata <- fromJSON(txt = targetMetadataPath)
  rownames(metadata) <- metadata$file_name
  metadata <- metadata[basename(tcga_files),]
  
  
  ### load TCGA clinical info
  clinInfo <- read.table(file = targetClinicalInfoPath, header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE, check.names = FALSE)
  rownames(clinInfo) <- clinInfo$submitter_id
  
  
  ### load TCGA exposure info
  exposureInfo <- read.table(file = targetExposureInfoPath, header = TRUE, sep = "\t",
                             stringsAsFactors = FALSE, check.names = FALSE)
  rownames(exposureInfo) <- exposureInfo$submitter_id
  
  
  ### merge the metadata and the sample info
  tcga_sample_info <- data.frame(tcga_sample_info,
                                 md5sum=metadata$md5sum,
                                 file_size=metadata$file_size,
                                 analysis_id=metadata$analysis$analysis_id,
                                 updated_datetime=metadata$analysis$updated_datetime,
                                 created_datetime=metadata$analysis$created_datetime,
                                 workflow_link=metadata$analysis$workflow_link,
                                 workflow_type=metadata$analysis$workflow_type,
                                 workflow_version=metadata$analysis$workflow_version,
                                 entity_id=sapply(metadata$associated_entities, function(x) x[1,1]),
                                 case_id=sapply(metadata$associated_entities, function(x) x[1,2]),
                                 barcode=sapply(metadata$associated_entities, function(x) x[1,3]),
                                 stringsAsFactors = FALSE, check.names = FALSE)
  rm(metadata)
  
  
  ### merge the clinical & exposure info and the sample info
  tempCommon <- intersect(rownames(clinInfo), rownames(exposureInfo))
  clin_exp_info <- data.frame(clinInfo[tempCommon,], exposureInfo[tempCommon,],
                              stringsAsFactors = FALSE, check.names = FALSE)
  clin_exp_info <- clin_exp_info[,-c(1, 3, 29, 30, 31)]
  tcga_sample_info <- merge.data.frame(tcga_sample_info, clin_exp_info,
                                       by.x = "Case ID", by.y = "submitter_id",
                                       all.x = TRUE)
  rownames(tcga_sample_info) <- tcga_sample_info$`File Name`
  tcga_sample_info <- tcga_sample_info[basename(tcga_files),]
  rm(clin_exp_info)
  rm(clinInfo)
  rm(exposureInfo)
  
  
  ### change all "--" to NA
  tcga_sample_info[which(tcga_sample_info == "--", arr.ind = TRUE)] <- NA
  
  
  ### collect unique project IDs
  project_ids <- unique(tcga_sample_info$`Project ID`)
  project_ids <- substring(project_ids, 6)
  project_ids <- tolower(project_ids)
  
  
  ### add additional info (RIN & FFPE) to the sample info based on "nationwidechildrens.org"
  tcga_sample_info <- data.frame(tcga_sample_info,
                                 is_derived_from_ffpe=NA,
                                 RIN=NA,
                                 stringsAsFactors = FALSE, check.names = FALSE)
  for(i in 1:length(project_ids)) {
    ### load the additional info for the given project ID
    add_info <- read.table(file = paste0(rinInfoDir, project_ids[i],
                                         "/bcr/biotab/clin/nationwidechildrens.org_biospecimen_analyte_",
                                         project_ids[i], ".txt"),
                           header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
    
    ### get indicies of the specific project ID
    idx <- which(tcga_sample_info$`Project ID` == paste0("TCGA-", toupper(project_ids[i])))
    
    ### for every project ID
    for(j in 1:length(idx)) {
      ### get index of the given sample in the add info
      targetIdx <- which(add_info$bcr_analyte_barcode ==
                         paste(strsplit(tcga_sample_info$barcode[idx[j]], split = "-", fixed = TRUE)[[1]][1:5], collapse = "-"))
      
      ### there are identical rows with the same bar code in the add info
      ### this is an issue of TCGA-XX-XXXX-01A-21R-A083-07
      ### only COAD & READ samples have this issue
      ### I found if there are duplicated rows, the second one is always the appropriate one
      ### after manually checking of the TCGA portal
      if(length(targetIdx) == 1) {
        ### fill out the FFPE & RIN column with the add info
        tcga_sample_info$is_derived_from_ffpe[idx[j]] <- add_info$is_derived_from_ffpe[targetIdx]
        tcga_sample_info$RIN[idx[j]] <- add_info$rinvalue[targetIdx]
      } else if(length(targetIdx) == 2) {
        ### fill out the FFPE & RIN column with the add info
        tcga_sample_info$is_derived_from_ffpe[idx[j]] <- add_info$is_derived_from_ffpe[targetIdx[2]]
        tcga_sample_info$RIN[idx[j]] <- add_info$rinvalue[targetIdx[2]]
      } else {
        stop(paste0("ERROR: There are too many duplicated rows in \"add_info\": ",
                    paste0("TCGA-", toupper(project_ids[i])), " ",
                    paste(targetIdx, collapse = "\t")))
      }
    }
  }
  rm(add_info)
  
  
  ### if RIN == "[Not Available]", then change it to NA
  naIdx <- which(tcga_sample_info$RIN == "[Not Available]")
  if(length(naIdx) > 0) {
    tcga_sample_info$RIN[naIdx] <- NA
  }
  
  
  ### set new row names for the sample info
  rownames(tcga_sample_info) <- tcga_sample_info$barcode
  
  
  ### read TCGA sample files
  htseq_raw_counts <- read.table(file = gzfile(tcga_files[1]), sep = "\t", row.names = 1, header = FALSE,
                                 stringsAsFactors = FALSE, check.names = FALSE)
  for(i in 2:length(tcga_files)) {
    htseq_raw_counts <- data.frame(htseq_raw_counts,
                                   read.table(file = gzfile(tcga_files[i]), sep = "\t", row.names = 1, header = FALSE,
                                              stringsAsFactors = FALSE, check.names = FALSE),
                                   stringsAsFactors = FALSE, check.names = FALSE)
  }
  colnames(htseq_raw_counts) <- rownames(tcga_sample_info)
  
  
  ### reorder everything based on the Project ID
  tcga_sample_info <- tcga_sample_info[order(tcga_sample_info$`Project ID`),]
  htseq_raw_counts <- htseq_raw_counts[,rownames(tcga_sample_info)]
  
  
  ### remove rows other than genes - they are at the end of the rows
  rIdx <- which(!startsWith(rownames(htseq_raw_counts), "ENS"))
  if(length(rIdx) > 0) {
    htseq_raw_counts <- htseq_raw_counts[-rIdx,]
  }
  
  
  ### set README function
  README <- function(){
    writeLines(paste(rep("#", 100), collapse = ""))
    writeLines("All the TCGA HT-SEQ raw counts were downloaded from the TCGA portal")
    writeLines("https://portal.gdc.cancer.gov/ on Dec 19, 2018")
    writeLines("There are 11093 files from 33 TCGA tissues")
    writeLines("The sample information are generated by combining info from four sources:")
    writeLines("sample sheet, metadata, clinical info, and biotab analyte info - they are all available from the portal as well")
    writeLines("The \"htseq_raw_counts\" has 60483 rows and 11093 columns")
    writeLines("The rows are Ensembl IDs and the columns are sample IDs")
    writeLines("The \"tcga_sample_info\" has 11093 rows and 21 columns")
    writeLines("The rows are samples corresponding to the columns of the \"htseq_raw_counts\"")
    writeLines("The columns represent various attributes of the samples")
    writeLines(paste(rep("#", 100), collapse = ""))
  }
  
  
  ### save the raw counts and the sample info as a RDA file
  save(list = c("htseq_raw_counts", "tcga_sample_info", "README"), file = paste0(outputDir, "TCGA_RAW_COUNTS.rda"))
  
  
  ### separate raw counts by each tissue
  project_ids <- unique(tcga_sample_info$`Project ID`)
  rcnt_matNames <- NULL
  for(i in 1: length(project_ids)) {
    ### collect rcnt data frame names
    rcnt_matNames <- c(rcnt_matNames, paste("rcnt", tolower(paste(strsplit(project_ids[i], split = "-", fixed = TRUE)[[1]], collapse = "_")), sep = "_"))
    
    ### make data frame object for each tissue
    assign(rcnt_matNames[i],
           htseq_raw_counts[,which(tcga_sample_info$`Project ID` == project_ids[i])],
           envir = globalenv())
  }
  
  
  ### set README function
  README <- function(){
    writeLines(paste(rep("#", 100), collapse = ""))
    writeLines("All the TCGA HT-SEQ raw counts were downloaded from the TCGA portal")
    writeLines("https://portal.gdc.cancer.gov/ on Dec 19, 2018")
    writeLines("There are 11093 files from 33 TCGA tissues")
    writeLines("The sample information are generated by combining info from four sources:")
    writeLines("sample sheet, metadata, clinical info, and biotab analyte info - they are all available from the portal as well")
    writeLines("The \"rcnt_tcga_[TISSUE_NAME]\" has 60483 rows")
    writeLines("There are 33 tissues and corresponding raw count data frames")
    writeLines("The rows are Ensembl IDs and the columns are sample IDs")
    writeLines("The \"tcga_sample_info\" has 11093 rows and 53 columns")
    writeLines("The rows are samples corresponding to the columns of the \"rcnt_tcga_[TISSUE_NAME]\"")
    writeLines("The columns represent various attributes of the samples")
    writeLines("The \"rcnt_matNames\" is a character vector that has the names of 33 rcnt_tcga_[TISSUE_NAME]s")
    writeLines(paste(rep("#", 100), collapse = ""))
  }
  
  
  ### make all the data frames into matrices before creating the RDA file
  for(i in 1:length(rcnt_matNames)) {
    temp <- get(rcnt_matNames[i])
    temp <- as.matrix(temp)
    assign(rcnt_matNames[i], temp, envir = globalenv())
  }
  tcga_sample_info <- as.matrix(tcga_sample_info)
  
  
  ### save the raw counts and the sample info as a RDA file
  save(list = c(rcnt_matNames, "rcnt_matNames", "tcga_sample_info", "README"), file = paste0(outputDir, "TCGA_33_RAW_COUNTS.rda"))
  
}

