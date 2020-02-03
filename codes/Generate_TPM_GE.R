### Generate TPM normalized gene expression data for 36 GTEx and 26 TCGA tissues

### Source the Utils.R file, if not already loaded.
if (Sys.info()["nodename"] %in% c("C2B2AFPL8", "C2B2AFPL9", "C2B2AFPL10", "C2B2AFPD11", "C2B2AFPL7")){
  UTILS_FILE = "C:/Users/floratos/workspace/R/labProjects/Common/Utils.R"
  ARACNE_GO_FILE = "C:/Users/floratos/workspace/R/labProjects/GTEx/R/aracne_go.R"
} else if (Sys.info()["nodename"] == "MDDSB-L004"){
  UTILS_FILE = "C:/Users/af2202/workspace/R/labProjects/Common/Utils.R"
  ARACNE_GO_FILE = "C:/Users/af2202/workspace/R/labProjects/GTEx/R/aracne_go.R"
} else if (Sys.info()["nodename"] == "C2B2AFPD9" || Sys.info()["nodename"] == "DESKTOP-F24420B"){
  UTILS_FILE = "C:/Research/CUMC/GTExProj/codes/Utils.R"
  ARACNE_GO_FILE = "C:/Research/CUMC/GTExProj/codes/aracne_go.R"
} else if (Sys.info()["nodename"] == "afdev5.c2b2.columbia.edu"){
  UTILS_FILE = "/ifs/home/c2b2/af_lab/floratos/cvs/labProjects/Common/Utils.R"
  ARACNE_GO_FILE = "/ifs/home/c2b2/af_lab/floratos/cvs/labProjects/GTEx/R/aracne_go.R"
} else if  (Sys.info()["user"]=='jb3401' || Sys.info()["user"]=='joshuabroyde'){
  source ('~/jb3401/scripts/labProjects/Common/Broyde/rFunctions.R')
  UTILS_FILE='~/jb3401/scripts/Utils/FloratosUtils.R'
} else if (Sys.info()["nodename"] == "C2B2AFPL6"){
  UTILS_FILE = "C:/Users/floratos/eclipse-workspace/R/labProjects/Common/Utils.R"
  ARACNE_GO_FILE = "C:/Users/floratos/eclipse-workspace/R/labProjects/GTEx/R/aracne_go.R"
} else if (any(grepl(Sys.info()["nodename"], c("C2B2AFPD6", "C2B2AFPD10", "C2B2AFPD12"), ignore.case = TRUE))) {
  if (Sys.info()["sysname"] == "Windows") {
    UTILS_FILE = "C:/repository/floratosLabCVS/labProjects/Common/Utils.R"
    ARACNE_GO_FILE = "C:/repository/floratosLabCVS/labProjects/GTEx/R/aracne_go.R"
  } else if (Sys.info()["sysname"] == "Linux"){
    UTILS_FILE = "/mnt/c/repository/floratosLabCVS/labProjects/Common/Utils.R"
    ARACNE_GO_FILE = "/mnt/c/repository/floratosLabCVS/labProjects/GTEx/R/aracne_go.R"
  }
} else {
  stop("No UTILS_FILE defined")
}

source(UTILS_FILE, chdir = TRUE)

### load the raw count data
load("//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/All_62_raw_counts.rda")

### separate indices for GTEx and TCGA
tcga_idx <- grep("tcga", rcntmat_names)
gtex_idx <- which(!grepl("tcga", rcntmat_names))

### make an empty list for the new TPM normalized ones
tpm_norm_cnt <- vector("list", length = length(rcntmat_names))
names(tpm_norm_cnt) <- rcntmat_names

### GTEx - hg19
for(tissue in rcntmat_names[gtex_idx]) {
  rcnt <- get(tissue)
  tpm_norm_cnt[[tissue]] <- normalizeRNASEQwithTPM(rcnt, gene_length_RDA_path = "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/exon_lengths_hg19.rda",
                                             ref = "hg19", filter_thresh = 1)
}

### TCGA - hg38
for(tissue in rcntmat_names[tcga_idx]) {
  rcnt <- get(tissue)
  tpm_norm_cnt[[tissue]] <- normalizeRNASEQwithTPM(rcnt, gene_length_RDA_path = "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/exon_lengths_hg38.rda",
                                                   ref = "hg38", filter_thresh = 1)
}

### polish the list name
names(tpm_norm_cnt) <- sapply(names(tpm_norm_cnt), function(x) {
  substring(x, 9)
})

### README function
README <- function(){
  writeLines(paste(rep("#", 100), collapse = ""))
  writeLines("The list \"tpm_norm_cnt\" contains TPM-normalized gene expression data")
  writeLines("of 36 GTEx and 26 TCGA tissues. The GTEx data were normalized based on")
  writeLines("the gene lengths of hg19 and the TCGA data were normalizd based on")
  writeLines("the gene lengths of hg38 reference human genome.")
  writeLines("Each entry of the list has a normalized data for each tissue, and")
  writeLines("names(tpm_norm_cnt) will show the tissue name of the given data.")
  writeLines(paste(rep("#", 100), collapse = ""))
}

### save the vipersig results as RDA file
save(list = c("tpm_norm_cnt", "README"), file = "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/All_62_TPM_normcnt.rda")
