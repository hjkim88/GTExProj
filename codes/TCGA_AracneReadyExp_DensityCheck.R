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
#               > densityCheck(AracneReadyRDAPath="./data/RDA_Files/TCGA_26_ARACNE_READY_EXPMAT.rda")
###

densityCheck() <- function(AracneReadyRDAPath="./data/RDA_Files/TCGA_26_ARACNE_READY_EXPMAT.rda") {
  
  ### load data
  load(AracneReadyRDAPath)
  
  
  ### print density plots
  for(i in 1:length(tcga_expmat_names)) {
    mat = get(tcga_expmat_names[i])
    plot(density(mat[,1]), main = tcga_expmat_names[i])
    apply(mat[, 2:ncol(mat)], 2, function(x){lines(density(x))})
  }
  
  
  ### we found tcga_esca and tcga_stad have deviating samples in their density plots
  ### we want to know what they are
  
  ### TCGA_ESCA
  tum_no <- 5
  mat = get(tcga_expmat_names[tum_no])
  
  ### variances check
  variances <- apply(mat, 2, var)
  barplot(variances, main = tcga_expmat_names[tum_no])
  
  ### deviating samples have relatively low variances than the others
  plot(density(mat[,1]), main = tcga_expmat_names[tum_no])
  apply(mat[, 2:ncol(mat)], 2, function(x){lines(density(x))})
  apply(mat[, order(variances)[1:10]], 2, function(x){lines(density(x), col = "RED")})
  
  ### compare sample info between deviating and normal samples
  target_sample_info <- tcga_sample_info[colnames(mat)[order(variances)[1:10]],]
  esca_sample_info <- tcga_sample_info[colnames(mat),]
  
  
  ### TCGA_STAD
  tum_no <- 22
  mat = get(tcga_expmat_names[tum_no])
  
  ### variances check
  variances <- apply(mat, 2, var)
  barplot(variances, main = tcga_expmat_names[tum_no])
  
  ### deviating samples have relatively low variances than the others
  plot(density(mat[,1]), main = tcga_expmat_names[tum_no])
  apply(mat[, 2:ncol(mat)], 2, function(x){lines(density(x))})
  apply(mat[, order(variances)[1:20]], 2, function(x){lines(density(x), col = "RED")})
  
  ### compare sample info between deviating and normal samples
  target_sample_info <- tcga_sample_info[colnames(mat)[order(variances)[1:20]],]
  stad_sample_info <- tcga_sample_info[colnames(mat),]
  
  
  ### Concluded that both in TCGA_ESCA and TCGA_STAD, there are no difference in sample info
  ### between the deviating and normal samples
  
}



