###
#   File name : PCAWithClusters.R
#   Author    : Hyunjin Kim
#   Date      : Jul 18, 2018
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Make PCA and t-SNE plots based on cluster info
#
#   Instruction
#               1. Source("PCAWithClusters.R")
#               2. Run the function "makePCA()" - specify the inputs (exp and cluster info) and output directory
#               3. PCA and t-SNE plots will be generated in the output directory
#
#   Example
#               > source("The_directory_of_PCAWithClusters.R/PCAWithClusters.R")
#               > makePCA(expPath="./results/aracne_ready/GTEx/vst_clean/Colon-Transverse_clean_vst.dat",
#                   clusterInfoPath="./results/viper/clustering/per_tissue/annotation/vmat_gtex_ColonTra_anno.txt",
#                   outputDir="./results/viper/clustering/per_tissue/annotation/")
###

makePCA <- function(expPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/ExpressionMatrices/GTEx_processed/aracne_ready_counts/vst_clean/Colon-Transverse_clean_vst.dat",
                    clusterInfoPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/viper/clustering/per_tissue/annotation/vmat_gtex_ColonTra_anno.txt",
                    outputDir="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/viper/clustering/per_tissue/annotation/") {
  
  ### load library
  if(!require(Rtsne)) {
    install.packages("Rtsne")
    library(Rtsne)
  }
  
  
  ### load datasets
  exp <- read.table(file = expPath, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
  clusterInfo <- read.table(file = clusterInfoPath, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
  
  
  ### A function to perform PCA and save a plot
  pca_plot <- function(normalizedMat, grp, title, filePath) {
    
    ### load library
    if(!require(ggfortify)) {
      install.packages("ggfortify")
      library(ggfortify)
    }
    
    ### PCA
    pca_result <- prcomp(t(normalizedMat))
    pca_group <- data.frame(pca_result$x, group=grp)
    
    colors = rainbow(length(unique(grp)))
    names(colors) = unique(grp)
    
    ### save as ggplot
    png(filename=filePath, width = 1000, height = 800)
    print(ggplot(pca_group,aes(x=PC1,y=PC2,col=group)) +
            labs(title=title) +
            geom_text(aes(label=colnames(normalizedMat)),hjust=0, vjust=0) +
            scale_color_manual(values = colors) +
            theme_classic(base_size = 16))
    dev.off()
    
  }
  
  
  ### A function to perform t-SNE and save a plot
  tsne_plot <- function(normalizedMat, grp, title, filePath) {
    
    ### load library
    if(!require(Rtsne)) {
      install.packages("Rtsne")
      library(Rtsne)
    }
    
    ### TSNE
    set.seed(1234)
    tryCatch({
      writeLines("Perplexity = 30")
      t <- Rtsne(t(normalizedMat), perplexity = 30)
    }, error = function(err) {
      tryCatch({
        writeLines("Perplexity = 10")
        t <- Rtsne(t(normalizedMat), perplexity = 10)
      }, error = function(err) {
        tryCatch({
          writeLines("Perplexity = 5")
          t <- Rtsne(t(normalizedMat), perplexity = 5)
        }, error = function(err) {
          tryCatch({
            writeLines("Perplexity = 3")
            t <- Rtsne(t(normalizedMat), perplexity = 3)
          }, error = function(err) {
            writeLines("Perplexity = 2")
            t <- Rtsne(t(normalizedMat), perplexity = 2)
          })
        })
      })
    })
    
    colors = rainbow(length(unique(grp)))
    names(colors) = unique(grp)
    
    ### save a plot
    png(filename=filePath, width = 1000, height = 800)
    plot(t$Y, col="white", xlab="tsne1", ylab="tsne2", main = title)
    text(t$Y, labels=colnames(normalizedMat), col=colors[grp])
    legend("topright", legend = unique(grp), col = colors[unique(grp)], pch = 15)
    dev.off()
    
  }
  
  
  ### title
  title <- paste(strsplit(basename(clusterInfoPath), split = "_", fixed = TRUE)[[1]][2:3], collapse = "_")
  title <- paste(title, "ClusterInfo", sep = "_")
  
  
  ### PCA
  pca_plot(normalizedMat = exp,
           grp = as.factor(clusterInfo[colnames(exp), "Group"]),
           title = paste("PCA", title, sep = "_"),
           filePath = paste0(outputDir, "PCA_", title, ".png"))
  
  
  ### t-SNE
  tsne_plot(normalizedMat = exp,
            grp = as.factor(clusterInfo[colnames(exp), "Group"]),
            title = paste("tSNE", title, sep = "_"),
            filePath = paste0(outputDir, "tSNE_", title, ".png"))
  
}



