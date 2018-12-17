###
#   File name : IdentifyingBatchEffect.R
#   Author    : Hyunjin Kim
#   Date      : Nov 3, 2017
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Performing PCA and TSNE to the datsets to identify batch effect
#
#   Instruction
#               1. Source("IdentifyingBatchEffect.R")
#               2. Run the function "batchEff()" - specify the input file (filtered counts) directory, sample information path, and output directory
#               3. PCA and t-SNE plots will be generated in the output directory
#
#   Example
#               > source("The_directory_of_IdentifyingBatchEffect.R/IdentifyingBatchEffect.R")
#               > batchEff(fCntPath="./results/aggregated_counts/",
#                          sampleInfoPath="./data/GTEx_Data_V6_SampleData.csv",
#                          outputPath="./results/batch_effect/")
###

batchEff <- function(fCntPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/ExpressionMatrices/GTEx_processed/separated_counts/",
                     sampleInfoPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/Annotations/GTEx_Data_V6_SampleData.csv",
                     outputPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/batch_effect/") {
  
  ### load library
  if(!require(stats)) {
    install.packages("stats")
    library(stats)
  }
  if(!require(Rtsne)) {
    install.packages("Rtsne")
    library(Rtsne)
  }
  
  
  ### load & preprocess sample info
  sampleInfo <- read.csv(sampleInfoPath, row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
  sampleInfo <- sampleInfo[,c("SMTSD", "SMNABTCH", "SMNABTCHT", "SMNABTCHD", "SMGEBTCH", "SMGEBTCHD", "SMGEBTCHT")]
  
  
  ### collect files from the fCntPath
  f <- list.files(fCntPath)
  
  
  ### A function for PCA plot
  plot_pca <- function(pc, group, outPath) {
    
    ### make plot list
    my_plots <- vector(ncol(group), mode = 'list')
    
    ### iteratively plot PCA with all the batches
    for(i in 1:ncol(group)) {
      
      ### set colors for PCA plot
      colors = rainbow(length(unique(group[,i])))
      names(colors) = unique(group[order(group[,i]),i])
      
      ### PCA plot
      plot(pc$x[,1], pc$x[,2], xlab="pc1", ylab="pc2", main = paste0("PCA_", colnames(group)[i]), col=colors[group[,i]])
      #text(pc$x[,1], pc$x[,2], labels=rownames(group), col=colors[group[,i]])
      if(length(unique(group[,i])) < 20) {
        legend("topright", legend = unique(group[,i]), col = colors[unique(group[,i])], pch = 15, cex = 0.7)
      }
      my_plots[[i]] <- recordPlot()
      graphics.off()
    }
    
    ### print the plots
    pdf(outPath, onefile = TRUE)
    for(my_plots in my_plots) {
      replayPlot(my_plots)
    }
    graphics.off()
  }
  
  
  ### A function for t-SNE plot
  plot_tsne <- function(tsne, group, outPath) {
    
    ### make plot list
    my_plots <- vector(ncol(group), mode = 'list')
    
    ### iteratively plot t-SNE with all the batches
    for(i in 1:ncol(group)) {
      
      ### set colors for t-SNE plot
      colors = rainbow(length(unique(group[,i])))
      names(colors) = unique(group[order(group[,i]),i])
      
      ### t-SNE plot
      plot(tsne$Y, xlab="tsne1", ylab="tsne2", main = paste0("t-SNE_", colnames(group)[i]), col=colors[group[,i]])
      #text(tsne$Y, labels=rownames(group), col=colors[group[,i]])
      if(length(unique(group[,i])) < 20) {
        legend("topright", legend = unique(group[,i]), col = colors[unique(group[,i])], pch = 15, cex = 0.7)
      }
      my_plots[[i]] <- recordPlot()
      graphics.off()
    }
    
    ### print the plots
    pdf(outPath, onefile = TRUE)
    for(my_plots in my_plots) {
      replayPlot(my_plots)
    }
    graphics.off()
  }
  
  
  ### iteratively perform PCA and t-SNE for all types of raw count data in the input path
  for(i in 1:length(f)) {
    
    ### load filtered counts
    fCnt <- read.table(paste0(fCntPath, f[i]), sep="\t", row.names = 1, header=TRUE, check.names = FALSE)
    
    ### set batch info
    group <- as.data.frame(as.character(sampleInfo[colnames(fCnt)[2:ncol(fCnt)], "SMTSD"]))
    rownames(group) <- colnames(fCnt)[2:ncol(fCnt)]
    colnames(group) <- "SMTSD"
    group$SMNABTCH <- as.character(sampleInfo[colnames(fCnt)[2:ncol(fCnt)], "SMNABTCH"])
    group$SMNABTCHT <- as.character(sampleInfo[colnames(fCnt)[2:ncol(fCnt)], "SMNABTCHT"])
    group$SMNABTCHD <- as.character(sampleInfo[colnames(fCnt)[2:ncol(fCnt)], "SMNABTCHD"])
    group$SMGEBTCH <- as.character(sampleInfo[colnames(fCnt)[2:ncol(fCnt)], "SMGEBTCH"])
    group$SMGEBTCHD <- as.character(sampleInfo[colnames(fCnt)[2:ncol(fCnt)], "SMGEBTCHD"])
    group$SMGEBTCHT <- as.character(sampleInfo[colnames(fCnt)[2:ncol(fCnt)], "SMGEBTCHT"])
    
    ### PCA with shifted log transformed data
    pca_result <- prcomp(t(log2(fCnt[,-1]+1)))
    
    ### plot PCA with the data
    plot_pca(pca_result, group, paste0(outputPath, substr(f[i], 1, nchar(f[i])-4), "_PCA.pdf"))
    cat(paste0(outputPath, substr(f[i], 1, nchar(f[i])-4), "_PCA.pdf\n"))
    
    
    tryCatch({
      
      ### t-SNE with shifted log transformed data
      set.seed(2990)
      tsne_result <- Rtsne(t(log2(fCnt[,-1]+1)), perplexity = 5)
      
      ### plot t-SNE with the data
      plot_tsne(tsne_result, group, paste0(outputPath, substr(f[i], 1, nchar(f[i])-4), "_TSNE.pdf"))
      
    }, warning = function(war) {
      cat(paste0(war, "\n"))
    }, error = function(err) {
      cat(paste0(err, "\n"))
    }, finally = {
      cat(paste0(outputPath, substr(f[i], 1, nchar(f[i])-4), "_TSNE.pdf\n"))
    })
    
  }
  
}

