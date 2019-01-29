### This code is for checking the old and the new TCGA Viper activity profiles
### Hyunjin Kim
### Jan 29, 2019

### load the old Viper profiles
load("./data/RDA_Files/TCGA_28_ViperMats.rda")

### change the variable name of the old data
varNamesVP_old <- NULL
for(i in 1:length(varNamesVP)) {
  varNamesVP_old <- c(varNamesVP_old, paste0("old_", varNamesVP[i]))
  assign(varNamesVP_old[i], get(varNamesVP[i]), envir = globalenv())
}

### remove the old loaded data
rm(list = varNamesVP)
rm(varNamesVP)

### garbage collection
gc()


### load the new Viper profiles
load("./data/RDA_Files/TCGA_26_ViperMats.rda")

### change the variable name of the new data
varNamesVP_new <- NULL
for(i in 1:length(varNamesVP)) {
  varNamesVP_new <- c(varNamesVP_new, paste0("new_", varNamesVP[i]))
  assign(varNamesVP_new[i], get(varNamesVP[i]), envir = globalenv())
}

### remove the old loaded data
rm(list = varNamesVP)
rm(varNamesVP)

### garbage collection
gc()


### only retain the same tissues between the old and the new
rm(list = varNamesVP_old[c(1, 12)])
varNamesVP_old <- varNamesVP_old[-c(1, 12)]


### generate density plots of correlation between old and new
for(i in 1:length(varNamesVP_new)) {
  X <- get(varNamesVP_new[i])
  Y <- get(varNamesVP_old[i])
  
  common_hubs <- intersect(rownames(X), rownames(Y))
  common_samples <- intersect(colnames(X), colnames(Y))
  
  X <- X[common_hubs, common_samples]
  Y <- Y[common_hubs, common_samples]
  
  correlation <- NULL
  for(j in 1:ncol(X)) {
    correlation <- c(correlation, cor(X[,j], Y[,j]))
  }
  names(correlation) <- common_samples
  
  png(paste0("./results/viper/check/", substring(varNamesVP_new[i], 10), ".png"),
      width = 1500, height = 1000, res = 130)
  plot(density(correlation),
       main = paste0("Correlations between old and new - ", substring(varNamesVP_new[i], 10)))
  dev.off()
}
