### This script is to generate t-SNE plots with 36 GTEx and 26 TCGA tissues
### the normalized Aracne-ready gene expressions should be used

### load the Aracne-ready gene expression data - VST-normalized
load("C:/Research/CUMC/GTExProj/data/RDA_Files/ALL_62_ARACNE_READY_EXPMAT.rda")

### load the TPM-normalized gene expression data
load("C:/Research/CUMC/GTExProj/data/RDA_Files/All_62_TPM_normcnt.rda")

### A function to perform t-SNE and save a plot
### normalizedMat: rows are genes and columns are samples
### grp: group information of the samples
### showNames: if TRUE, sample names will be presented in the plot
###            if FALSE, only dots will be presented (default: FALSE)
### num: the number of top genes to be used based on variance (-1 [default]: use all the genes)
### title: title of the plot
### outDir: output directory for the plot
tsne_plot <- function(normalizedMat, grp, showNames = FALSE, num = -1,
                      title="TSNE_Plot", outDir="./") {
  ### load library
  if(!require(Rtsne)) {
    install.packages("Rtsne")
    library(Rtsne)
  }
  
  ### select the top genes based on variance
  if(num >= 0 && num <= nrow(normalizedMat)) {
    v <- apply(normalizedMat, 1, var)
    v <- v[order(-v)]
    top_genes <- names(v)[1:num]
  } else {
    top_genes <- rownames(normalizedMat)
  }
  
  ### TSNE
  set.seed(1234)
  tryCatch({
    writeLines("Perplexity = 30")
    t <- Rtsne(t(normalizedMat[top_genes,]), perplexity = 30)
  }, error = function(err) {
    tryCatch({
      writeLines("Perplexity = 10")
      t <- Rtsne(t(normalizedMat[top_genes,]), perplexity = 10)
    }, error = function(err) {
      tryCatch({
        writeLines("Perplexity = 5")
        t <- Rtsne(t(normalizedMat[top_genes,]), perplexity = 5)
      }, error = function(err) {
        tryCatch({
          writeLines("Perplexity = 3")
          t <- Rtsne(t(normalizedMat[top_genes,]), perplexity = 3)
        }, error = function(err) {
          writeLines("Perplexity = 2")
          t <- Rtsne(t(normalizedMat[top_genes,]), perplexity = 2)
        })
      })
    })
  })
  
  colors = rainbow(length(unique(grp)))
  names(colors) = unique(grp)
  
  ### save a plot
  png(filename=paste0(outDir, title, ".png"), width = 1000, height = 800)
  if(showNames) {
    plot(t$Y, col="white", xlab="tsne1", ylab="tsne2", main = title)
    text(t$Y, labels=colnames(normalizedMat), col=colors[grp])
  } else {
    plot(t$Y, col=colors[grp], pch=19, xlab="tsne1", ylab="tsne2", main = title)
  }
  legend("topleft", legend = unique(grp), col = colors[unique(grp)], pch = 15)
  dev.off()
}

### GTEx tissues
### create the directory
### VST
result_path <- "C:/Research/CUMC/GTExProj/results/TSNE/GTEX/VST/"
dir.create(path = result_path, showWarnings = FALSE, recursive = TRUE)

### refine the GTEx tissue names
gtex_tissue_names <- sapply(gtex_expmat_names, function(x) {
  strsplit(x, split = "_", fixed = TRUE)[[1]][3]
})

### generate tsne for each GTEx tissue
for(i in 1:length(gtex_expmat_names)) {
  
  ### get the gene expressions for the tissue
  ge <- get(gtex_expmat_names[i])
  
  ### generate a t-SNE plot
  tsne_plot(ge, grp = rep(gtex_tissue_names[i], ncol(ge)),
            num = 1000, title = paste0("TSNE_Plot_", gtex_tissue_names[i], "_VST.png"),
            outDir = result_path)
  
}

### TPM
result_path <- "C:/Research/CUMC/GTExProj/results/TSNE/GTEX/TPM/"
dir.create(path = result_path, showWarnings = FALSE, recursive = TRUE)

### generate tsne for each GTEx tissue
for(i in 1:36) {
  
  ### generate a t-SNE plot
  tsne_plot(tpm_norm_cnt[[i]], grp = names(tpm_norm_cnt)[i],
            num = 1000, title = paste0("TSNE_Plot_", names(tpm_norm_cnt)[i], "_TPM.png"),
            outDir = result_path)
  
}



### TCGA tissues
### create the directory
### VST
result_path <- "C:/Research/CUMC/GTExProj/results/TSNE/TCGA/VST/"
dir.create(path = result_path, showWarnings = FALSE, recursive = TRUE)

### refine the TCGA tissue names
tcga_tissue_names <- sapply(tcga_expmat_names, function(x) {
  strsplit(x, split = "_", fixed = TRUE)[[1]][3]
})

### generate tsne for each TCGA tissue
for(i in 1:length(tcga_expmat_names)) {
  
  ### get the gene expressions for the tissue
  ge <- get(tcga_expmat_names[i])
  
  ### generate a t-SNE plot
  tsne_plot(ge, grp = rep(tcga_tissue_names[i], ncol(ge)),
            num = 1000, title = paste0("TSNE_Plot_", tcga_tissue_names[i], "_VST.png"),
            outDir = result_path)
  
}

### TPM
result_path <- "C:/Research/CUMC/GTExProj/results/TSNE/TCGA/TPM/"
dir.create(path = result_path, showWarnings = FALSE, recursive = TRUE)

### generate tsne for each GTEx tissue
for(i in 37:62) {
  
  ### generate a t-SNE plot
  tsne_plot(tpm_norm_cnt[[i]], grp = names(tpm_norm_cnt)[i],
            num = 1000, title = paste0("TSNE_Plot_", names(tpm_norm_cnt)[i], "_TPM.png"),
            outDir = result_path)
  
}
