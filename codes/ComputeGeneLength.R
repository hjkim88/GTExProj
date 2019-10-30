### this script is for computing gene length for all the existing Entrez IDs

### load library
if(!require(org.Hs.eg.db, quietly = TRUE)) {
  if(!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("org.Hs.eg.db", version = "3.8")
  require(org.Hs.eg.db, quietly = TRUE)
}
if(!require(EDASeq, quietly = TRUE)) {
  if(!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("EDASeq", version = "3.8")
  require(EDASeq, quietly = TRUE)
}

### get Entrez ID list
eg_list <- mappedkeys(org.Hs.egSYMBOL)

### make an empty gene length vector
gene_lengths <- rep(NA, length(eg_list))
names(gene_lengths) <- eg_list

### compute gene length
gene_lengths <- getGeneLengthAndGCContent(as.character(eg_list), org = "hg38", mode = "org.db")[,"length"]

### set README function
README <- function(){
  writeLines(paste(rep("#", 100), collapse = ""))
  writeLines("The \"gene_lengths\" is a numeric vector that contains gene lengths (bp)")
  writeLines("of all the existing Entrez IDs. They were computed by using")
  writeLines("\"getGeneLengthAndGCContent()\" function of EDASeq package.")
  writeLines(paste(rep("#", 100), collapse = ""))
}

### save as RDA
save(list = c("gene_lengths", "README"), file = "./data/RDA_Files/Gene_Lengths.rda")
