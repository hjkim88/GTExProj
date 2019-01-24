### one-time script for combine tfNetEnrich of each tissue

### set working directory
wdPath="C:/Research/CUMC/GTExProj/data/RDA_Files/tfNetEnrich/"

### get RDA file nams in the directory
f <- list.files(wdPath)
f <- f[which(endsWith(f, ".rda"))]

### create the results list (a named list),  with one entry per interactome
tfNetEnrich <- vector("list", length = length(f))
names(tfNetEnrich) <- sapply(f, function(x) {
  substr(x, 13, nchar(x)-4)
})
names(f) <- names(tfNetEnrich)

### load RDA (if there is only one object in RDA)
loadRDA <- function(file) {
  env <- new.env()
  load(file = file, envir = env)
  env[[ls(env)[1]]]
}

### load partial tfNetEnrich RDA files and save only PVs
for(net in names(tfNetEnrich)) {
  
  ### print progress
  writeLines(paste("\tProcess ->", net))
  
  ### load the RDA
  rda_result <- loadRDA(paste0(wdPath, f[net]))
  
  ### sort the RDA file
  rda_result <- rda_result[order(as.integer(names(rda_result)))]
  rda_result <- lapply(rda_result, function(x) {
    return(x[order(as.integer(rownames(x))),])
  })
  
  ### create an empty matrix for current interactome
  tfNetEnrich[[net]] <- matrix(NA, nrow = length(rda_result), ncol = length(rda_result))
  rownames(tfNetEnrich[[net]]) <- names(rda_result)
  colnames(tfNetEnrich[[net]]) <- names(rda_result)
  
  ### save the results in the matrix from the RDA file
  for(i in 1:nrow(tfNetEnrich[[net]])) {
    tfNetEnrich[[net]][i,-i] <- rda_result[[i]][,"-log10_pval"]
  }
  
  ### tfNetEnrich[[net]][i,i] = Inf
  diag(tfNetEnrich[[net]]) <- Inf
  
  ### garbage collection
  gc()
  
}

### save the result as a RDA file
save(list = c("tfNetEnrich"), file = "C:/Research/CUMC/GTExProj/data/RDA_Files/All_62_tfNetEnrich.rda")
