###
#   File name : TCGA_AracneReadyExp_DensityCheck.R
#   Author    : Hyunjin Kim
#   Date      : Jan 11, 2019
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Before Aracne run, we would like to mare sure that there is nothing wrong with
#               the Aracne-ready expression files. We already created PCA/TSNE plots and now
#               we also want to see density plots of the expressions
#
#   Instruction
#               1. Source("TCGA_AracneReadyExp_DensityCheck.R")
#               2. Run the function "densityCheck()" - specify the input file (Aracne-ready exp RDA)
#               3. Analysis results will be shown in the console
#
#   Example
#               > source("The_directory_of_TCGA_AracneReadyExp_DensityCheck.R/TCGA_AracneReadyExp_DensityCheck.R")
#               > densityCheck(AracneReadyRDAPath="./data/RDA_Files/TCGA_26_ARACNE_READY_EXPMAT.rda",
#                              rawCntRDAPath="./data/RDA_Files/TCGA_33_RAW_COUNTS.rda")
###

densityCheck() <- function(AracneReadyRDAPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/TCGA_26_ARACNE_READY_EXPMAT.rda",
                           rawCntRDAPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/TCGA_33_RAW_COUNTS.rda") {
  
  ### load library
  if(!require(ggplot2)) {
    install.packages("ggplot2")
    library(ggplot2)
  }
  if(!require(ggbeeswarm)) {
    install.packages("ggbeeswarm")
    library(ggbeeswarm)
  }
  if(!require(ggpubr)) {
    install.packages("ggpubr")
    library(ggpubr)
  }
  
  
  ### load data
  load(rawCntRDAPath)
  load(AracneReadyRDAPath)
  
  
  ### a function to perform both PCA and t-SNE
  pca_tsne_plot <- function(normalizedMat, grp, title, filePath=NULL) {
    
    ### load libraries
    if(!require(ggfortify)) {
      install.packages("ggfortify")
      library(ggfortify)
    }
    if(!require(Rtsne)) {
      install.packages("Rtsne")
      library(Rtsne)
    }
    if(!require(gridExtra)) {
      install.packages("gridExtra")
      library(gridExtra)
    }
    
    ### PCA
    pca_result <- prcomp(t(normalizedMat))
    pca_group <- data.frame(pca_result$x, group=grp)
    
    ### TSNE
    env <- new.env()
    set.seed(1234)
    tryCatch({
      writeLines("Perplexity = 30")
      t <- Rtsne(t(normalizedMat), perplexity = 30)
      assign("t", t, envir = env)
    }, error = function(err) {
      tryCatch({
        writeLines("Perplexity = 10")
        t <- Rtsne(t(normalizedMat), perplexity = 10)
        assign("t", t, envir = env)
      }, error = function(err) {
        tryCatch({
          writeLines("Perplexity = 5")
          t <- Rtsne(t(normalizedMat), perplexity = 5)
          assign("t", t, envir = env)
        }, error = function(err) {
          tryCatch({
            writeLines("Perplexity = 3")
            t <- Rtsne(t(normalizedMat), perplexity = 3)
            assign("t", t, envir = env)
          }, error = function(err) {
            writeLines("Perplexity = 2")
            t <- Rtsne(t(normalizedMat), perplexity = 2)
            assign("t", t, envir = env)
          })
        })
      })
    })
    t <- get("t", envir = env)
    tsne_group <- data.frame(Dimension1=t$Y[,1], Dimension2=t$Y[,2], group=grp)
    
    ### set colors for the samples
    colors = rainbow(length(unique(grp)))
    names(colors) = unique(grp)
    
    ### print out the results
    pca_plot <- ggplot(pca_group,aes(x=PC1,y=PC2,col=group)) +
      labs(title=paste("PCA", title, sep = "_")) +
      geom_text(aes(label=colnames(normalizedMat)),hjust=0, vjust=0) +
      scale_color_manual(values = colors) +
      theme_classic(base_size = 16)
    tsne_plot <- ggplot(tsne_group, aes(x=Dimension1,y=Dimension2,col=group)) +
      labs(title=paste("TSNE", title, sep = "_")) +
      geom_text(aes(label=colnames(normalizedMat)),hjust=0, vjust=0) +
      scale_color_manual(values = colors) +
      theme_classic(base_size = 16)
    if(is.null(filePath)) {
      grid.arrange(pca_plot, tsne_plot, ncol = 2)
    } else {
      ggsave(filename = filePath, arrangeGrob(pca_plot, tsne_plot, ncol = 2), width = 20, height = 10)
    }
    
  }
  
  
  ### print density plots
  for(i in 1:length(tcga_expmat_names)) {
    mat <- get(tcga_expmat_names[i])
    plot(density(mat[,1]), main = tcga_expmat_names[i])
    apply(mat[, 2:ncol(mat)], 2, function(x){lines(density(x))})
  }
  
  
  ### we found tcga_esca and tcga_stad have deviating samples in their density plots
  ### we want to know what they are
  
  ### TCGA_ESCA
  
  ### variances check
  variances <- apply(expmat_tcga_esca, 2, var)
  barplot(variances, main = "TCGA_ESCA")
  
  ### deviating samples have relatively low variances than the others
  plot(density(expmat_tcga_esca[,1]), main = "TCGA_ESCA")
  apply(expmat_tcga_esca[, 2:ncol(expmat_tcga_esca)], 2, function(x){lines(density(x))})
  apply(expmat_tcga_esca[, order(variances)[1:10]], 2, function(x){lines(density(x), col = "RED")})
  
  ### compare sample info between deviating and normal samples
  target_sample_info <- tcga_sample_info[colnames(expmat_tcga_esca)[order(variances)[1:10]],]
  esca_sample_info <- tcga_sample_info[colnames(expmat_tcga_esca),]
  
  
  ### now check depth of coverage (total number of reads) of the deviating samples
  target_no.reads <- apply(rcnt_tcga_esca[,target_sample_info[,"barcode"]], 2, sum)
  other_no.reads <- apply(rcnt_tcga_esca[,setdiff(esca_sample_info[,"barcode"],
                                                  target_sample_info[,"barcode"])], 2, sum)
  avg.no.reads <- mean(c(target_no.reads, other_no.reads))
  
  ### beeswarm plot
  df <- data.frame(Type=c(rep("Deviating_Samples", length(target_no.reads)),
                          rep("The_Others", length(other_no.reads))),
                   No_Reads=c(target_no.reads, other_no.reads),
                   stringsAsFactors = FALSE, check.names = FALSE)
  ggplot(df, aes(x=Type, y=No_Reads)) +
    labs(title="TCGA_ESCA", y="Total number of reads") +
    theme_classic(base_size = 16) +
    geom_boxplot() +
    geom_beeswarm(aes(color=Type)) +
    stat_compare_means()
  
  ### PCA/TSNE plot with the deviating samples
  grp <- rep("The_Others", ncol(expmat_tcga_esca))
  grp[order(variances)[1:10]] <- "Deviating_Samples"
  pca_tsne_plot(normalizedMat = expmat_tcga_esca,
                grp = grp,
                title = "TCGA_ESCA",
                filePath = NULL)
   
  
  ### TCGA_STAD
  
  ### variances check
  variances <- apply(expmat_tcga_stad, 2, var)
  barplot(variances, main = "TCGA_STAD")
  
  ### deviating samples have relatively low variances than the others
  plot(density(expmat_tcga_stad[,1]), main = "TCGA_STAD")
  apply(expmat_tcga_stad[, 2:ncol(expmat_tcga_stad)], 2, function(x){lines(density(x))})
  apply(expmat_tcga_stad[, order(variances)[1:20]], 2, function(x){lines(density(x), col = "RED")})
  
  ### compare sample info between deviating and normal samples
  target_sample_info <- tcga_sample_info[colnames(expmat_tcga_stad)[order(variances)[1:20]],]
  stad_sample_info <- tcga_sample_info[colnames(expmat_tcga_stad),]
  
  ### now check depth of coverage (total number of reads) of the deviating samples
  target_no.reads <- apply(rcnt_tcga_stad[,target_sample_info[,"barcode"]], 2, sum)
  other_no.reads <- apply(rcnt_tcga_stad[,setdiff(stad_sample_info[,"barcode"],
                                                  target_sample_info[,"barcode"])], 2, sum)
  avg.no.reads <- mean(c(target_no.reads, other_no.reads))
  
  ### beeswarm plot
  df <- data.frame(Type=c(rep("Deviating_Samples", length(target_no.reads)),
                          rep("The_Others", length(other_no.reads))),
                   No_Reads=c(target_no.reads, other_no.reads),
                   stringsAsFactors = FALSE, check.names = FALSE)
  ggplot(df, aes(x=Type, y=No_Reads)) +
    labs(title="TCGA_STAD", y="Total number of reads") +
    theme_classic(base_size = 16) +
    geom_boxplot() +
    geom_beeswarm(aes(color=Type)) +
    stat_compare_means()
  
  ### PCA/TSNE plot with the deviating samples
  grp <- rep("The_Others", ncol(expmat_tcga_stad))
  grp[order(variances)[1:20]] <- "Deviating_Samples"
  pca_tsne_plot(normalizedMat = expmat_tcga_stad,
                grp = grp,
                title = "TCGA_STAD",
                filePath = NULL)
  

  ### Concluded that both in TCGA_ESCA and TCGA_STAD, there are no difference in sample info
  ### between the deviating and normal samples
  
}



