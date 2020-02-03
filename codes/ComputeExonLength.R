### Compute Exon lengths that will be used in RNA-Seq analysis

### load necessary libraries
if(!require(org.Hs.eg.db, quietly = TRUE)) {
  if(!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("org.Hs.eg.db", version = "3.8")
  require(org.Hs.eg.db, quietly = TRUE)
}
if(!require(TxDb.Hsapiens.UCSC.hg19.knownGene, quietly = TRUE)) {
  if(!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene", version = "3.8")
  require(TxDb.Hsapiens.UCSC.hg19.knownGene, quietly = TRUE)
}
if(!require(TxDb.Hsapiens.UCSC.hg38.knownGene, quietly = TRUE)) {
  if(!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene", version = "3.8")
  require(TxDb.Hsapiens.UCSC.hg38.knownGene, quietly = TRUE)
}

### hg19
### get exon length of a gene (entrez_id) based on hg19 reference
hg19GeneLengths <- function(entrez_ids)
{
  exons.db = exonsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, by='gene')    

  sapply(entrez_ids, function(eg)
  {
    exons = exons.db[[eg]]
    if(is.null(exons)) {
      return(NA)
    } else {
      exons = reduce(exons)
      return(sum(width(exons)))
    }
  })
}

### hg38
### get exon length of a gene (entrez_id) based on hg38 reference
hg38GeneLengths <- function(entrez_ids)
{
  exons.db = exonsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, by='gene')    
  
  sapply(entrez_ids, function(eg)
  {
    exons = exons.db[[eg]]
    if(is.null(exons)) {
      return(NA)
    } else {
      exons = reduce(exons)
      return(sum(width(exons)))
    }
  })
}

### get Entrez ID list
eg_list <- mappedkeys(org.Hs.egSYMBOL)


### hg19
### compute exon lengths of all the existing genes based on hg19
exon_lengths_hg19 <- hg19GeneLengths(eg_list)

### hg38
### compute exon lengths of all the existing genes based on hg38
exon_lengths_hg38 <- hg38GeneLengths(eg_list)

### hg19
### set README function
README <- function(){
  writeLines(paste(rep("#", 100), collapse = ""))
  writeLines("The \"exon_lengths_hg19\" is a numeric vector that contains gene (exon) lengths (bp)")
  writeLines("of all the existing Entrez IDs. They were computed by using")
  writeLines("\"TxDb.Hsapiens.UCSC.hg19.knownGene\" package.")
  writeLines("They were computed based on hg19 genome reference.")
  writeLines(paste(rep("#", 100), collapse = ""))
}

### save as RDA
save(list = c("exon_lengths_hg19", "README"), file = "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/exon_lengths_hg19.rda")


### hg38
### set README function
README <- function(){
  writeLines(paste(rep("#", 100), collapse = ""))
  writeLines("The \"exon_lengths_hg38\" is a numeric vector that contains gene (exon) lengths (bp)")
  writeLines("of all the existing Entrez IDs. They were computed by using")
  writeLines("\"TxDb.Hsapiens.UCSC.hg38.knownGene\" package.")
  writeLines("They were computed based on hg38 genome reference.")
  writeLines(paste(rep("#", 100), collapse = ""))
}

### save as RDA
save(list = c("exon_lengths_hg38", "README"), file = "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/exon_lengths_hg38.rda")
