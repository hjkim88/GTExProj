
### Jul 11, 2018
### Hyunjin Kim
### AnnotateViperClusters.R

### Needs GTEx_36_ViperMats.rda
load("./GTEx_36_ViperMats.rda")

### Needs
annoTable <- read.csv(file = "./data/GTEx_Data_V6_SampleData.csv", header = TRUE,
                      row.names = 1, check.names = FALSE)

### set tissues of interest
interest <- c("vmat_gtex_CellsTra", "vmat_gtex_ColonTra", "vmat_gtex_Stomach")

### set output directory
outputDir <- "./results/viper/clustering/per_tissue/annotation/"

### iteratively perform the annotation
for(i in 1:length(interest)) {
  ### get viper matrix of each tissue
  viperMat <- get(interest[i])
  
  ### get Euclidean distance
  d <- dist(t(viperMat))
  
  ### hierarchical clustering
  h <- hclust(d, method = "average")
  
  # ### get clustered sample names
  # clustered_samples <- h$labels[h$order]
  # 
  # ### get annotation
  # annoResult <- annoTable[clustered_samples,]
  # 
  # ### set sample name column
  # annoResult <- cbind(Sample_Name=rownames(annoResult), annoResult)
  # 
  # ### save the annoResult
  # write.table(annoResult, file = paste0(outputDir, interest[i], "_anno.txt"),
  #             sep = "\t", row.names = FALSE)
  
  ### group samples into two
  clusterCut <- cutree(h,2)
  
  ### get annotation
  annoResult <- annoTable[names(clusterCut[c(which(clusterCut == 1), which(clusterCut == 2))]),]
  
  ### add group column
  annoResult <- cbind(Group=c(rep(1, length(which(clusterCut == 1))), rep(2, length(which(clusterCut == 2)))), annoResult)
  
  ### set sample name column
  annoResult <- cbind(Sample_Name=rownames(annoResult), annoResult)
  
  ### save the annoResult
  write.table(annoResult, file = paste0(outputDir, interest[i], "_anno.txt"),
              sep = "\t", row.names = FALSE)
  
  ### find which column has the most differential values between the two clusters
  ### get the numeric column
  numericIdx <- which(sapply(annoResult, is.numeric))
  ### remove NA columns
  temp <- sapply(annoResult, is.na)
  naIdx <- which(apply(matrix(as.numeric(temp), nrow(temp), ncol(temp)), 2, sum) == nrow(temp))
  numericIdx <- setdiff(numericIdx, naIdx)
  
  ### perform t-test for each numeric column
  result <- apply(annoResult[,numericIdx], 2, function(x) {
    try(t.test(as.numeric(x[which(annoResult$Group == 1)]), as.numeric(x[which(annoResult$Group == 2)])), silent = TRUE)
  })
  
  ### organize the result for print
  # col 1: t
  # col 2: pv
  # col 3: mean of x
  # col 4: mean of y
  # rows : colnames
  printResult <- matrix(NA, length(numericIdx), 4)
  colnames(printResult) <- c("t-statistics", "p-value", "group1-mean", "group2-mean")
  rownames(printResult) <- colnames(annoResult)[numericIdx]
  for(j in 1:length(numericIdx)) {
    if(length(result[[j]]) > 1) {
      printResult[j,1] <- result[[j]]$statistic
      printResult[j,2] <- result[[j]]$p.value
      printResult[j,3] <- result[[j]]$estimate[1]
      printResult[j,4] <- result[[j]]$estimate[2]
    }
  }
  
  ### order the printResult based on p-value - asscending order
  printResult <- printResult[order(printResult[,2]),]
  
  ### set row name column
  printResult <- cbind(ColName=rownames(printResult), printResult)
  
  ### save the printResult
  write.table(printResult, file = paste0(outputDir, interest[i], "_anno_diff.txt"),
              sep = "\t", row.names = FALSE)
}



