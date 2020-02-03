### We found that mutual information values are smaller in TCGA Aracne networks than in GTEx networks
### and we think that this is because TCGA gene expressions are more heterogeneous
### so we want to measure heterogeneity of every tissue and compare the order
### with mutatual information

### load the Aracne network data
load("//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/All_62_ARACNE.rda")

### Aris' code
### get a list with one entry per interactome, each entry being a vector of interactome's
### interaction MI values
getNetMis <- function(nets = varNames){
  res = lapply(nets, function(net){
    return(unlist(lapply(get(net)[[2]], function(reg){
      return(reg[reg[, "Target"] > 0, "MI"])
    })))
  })
}

mis = getNetMis(varNames) 		# a little time consuming
# Remove run-through expression
mis = lapply(mis, function(x){return(x[x < 2])})
names(mis) = varNames
# Find the 5% quantiles for each MI series
mi_q = t(sapply(mis, function(m){quantile(m, probs = seq(0.05, 1, 0.05))}))
mi_q1 = apply(mi_q, 1, function(qs){
  res = rep(qs[1], length(qs))
  for (i in 2:length(qs))
    res[i] = qs[i] - qs[i-1]
  return(res)
})
rownames(mi_q1) = colnames(mi_q)
barplot(mi_q1, las=2)
# Sort by 95% percentile
ind = order(mi_q[, "95%"])
par(mar= c(10, 5, 5, 3))
barplot(mi_q1[, ind], las=2, main = "Mutual Information - 5% Quantiles")


### VST normalized data
load("//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/ALL_62_ARACNE_READY_EXPMAT.rda")

### mean of Euclidean distance
ed_dispersion <- rep(0, 62)
names(ed_dispersion)[1:36] <- gtex_expmat_names
names(ed_dispersion)[37:62] <- tcga_expmat_names
for(gtex_name in gtex_expmat_names) {
  ed_dispersion[gtex_name] <- mean(dist(t(get(gtex_name))))
  gc()
}
for(tcga_name in tcga_expmat_names) {
  ed_dispersion[tcga_name] <- mean(dist(t(get(tcga_name))))
  gc()
}
# save(list = c("ed_dispersion"), file = "./ed_dispersion.rda")

names(ed_dispersion) <- sapply(names(ed_dispersion), function(x) {
  if(grepl("tcga", x)) {
    return(substring(x, 8))
  } else {
    return(substring(x, 13))
  }
})

ed_dispersion <- ed_dispersion[colnames(mi_q1[, ind])]

par(mar= c(10, 5, 5, 3))
barplot(ed_dispersion, las=2, main = "Euclidean Distance of Gene expressions")


### determinant of covariance matrix
det_cov <- rep(0, 62)
names(det_cov)[1:36] <- gtex_expmat_names
names(det_cov)[37:62] <- tcga_expmat_names
for(gtex_name in gtex_expmat_names) {
  det_cov[gtex_name] <- det(cov(get(gtex_name)))
  gc()
}
for(tcga_name in tcga_expmat_names) {
  det_cov[tcga_name] <- det(cov(get(tcga_name)))
  gc()
}

names(det_cov) <- sapply(names(det_cov), function(x) {
  if(grepl("tcga", x)) {
    return(substring(x, 8))
  } else {
    return(substring(x, 13))
  }
})

det_cov <- det_cov[colnames(mi_q1[, ind])]

par(mar= c(10, 5, 5, 3))
barplot(-log2(det_cov), las=2, main = "-log2(Determinant of Covariance Matrix)")


### Euclidean on the t-SNE

### load library
if(!require(Rtsne)) {
  install.packages("Rtsne")
  library(Rtsne)
}

### a function to get mean distance on TSNE
getMeanD_tsne <- function(normalizedMat, num = 1000) {
  
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
  
  return(mean(dist(t$Y)))
  
}

ed_tsne <- rep(0, 62)
names(ed_tsne)[1:36] <- gtex_expmat_names
names(ed_tsne)[37:62] <- tcga_expmat_names
for(gtex_name in gtex_expmat_names) {
  ed_tsne[gtex_name] <- getMeanD_tsne(get(gtex_name))
  gc()
}
for(tcga_name in tcga_expmat_names) {
  ed_tsne[tcga_name] <- getMeanD_tsne(get(tcga_name))
  gc()
}

names(ed_tsne) <- sapply(names(ed_tsne), function(x) {
  if(grepl("tcga", x)) {
    return(substring(x, 8))
  } else {
    return(substring(x, 13))
  }
})

ed_tsne <- ed_tsne[colnames(mi_q1[, ind])]

par(mar= c(10, 5, 5, 3))
barplot(ed_tsne, las=2, main = "Euclidean Distance on the t-SNE")



# ### The Weighted Quantile-Adjusted Conditional Maximum Likelihood (wqCML)
# library(edgeR)
# 
# edger_dispersion <- rep(0, 62)
# names(edger_dispersion)[1:36] <- gtex_expmat_names
# names(edger_dispersion)[37:62] <- tcga_expmat_names
# for(gtex_name in gtex_expmat_names) {
#   edger_dispersion[gtex_name] <- estimateCommonDisp(get(gtex_name))
#   gc()
# }
# for(tcga_name in tcga_expmat_names) {
#   edger_dispersion[tcga_name] <- estimateCommonDisp(get(tcga_name))
#   gc()
# }
# 
# barplot(edger_dispersion, las=2)


### TPM normalized data
load("//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/All_62_TPM_normcnt.rda")

### mean of Euclidean distance
ed_dispersion <- rep(0, length(tpm_norm_cnt))
names(ed_dispersion) <- names(tpm_norm_cnt)
for(tissue_name in names(ed_dispersion)) {
  ed_dispersion[tissue_name] <- mean(dist(t(tpm_norm_cnt[[tissue_name]])))
  gc()
}
# save(list = c("ed_dispersion"), file = "./ed_dispersion_tpm.rda")


ed_dispersion <- ed_dispersion[colnames(mi_q1[, ind])]

par(mar= c(10, 5, 5, 3))
barplot(ed_dispersion, las=2, main = "Euclidean Distance of Gene expressions (TPM)")


# ### determinant of covariance matrix
# det_cov <- rep(0, length(tpm_norm_cnt))
# names(det_cov) <- names(tpm_norm_cnt)
# for(tissue_name in names(det_cov)) {
#   det_cov[tissue_name] <- det(cov(tpm_norm_cnt[[tissue_name]]))
#   gc()
# }
# 
# det_cov <- det_cov[colnames(mi_q1[, ind])]
# 
# par(mar= c(10, 5, 5, 3))
# barplot(-log2(det_cov), las=2, main = "-log2(Determinant of Covariance Matrix) (TPM)")


### Euclidean on the t-SNE
ed_tsne <- rep(0, length(tpm_norm_cnt))
names(ed_tsne) <- names(tpm_norm_cnt)
for(tissue_name in names(ed_tsne)) {
  ed_tsne[tissue_name] <- getMeanD_tsne(tpm_norm_cnt[[tissue_name]])
  gc()
}

ed_tsne <- ed_tsne[colnames(mi_q1[, ind])]

par(mar= c(10, 5, 5, 3))
barplot(ed_tsne, las=2, main = "Euclidean Distance on the t-SNE (TPM)")



# ### The Weighted Quantile-Adjusted Conditional Maximum Likelihood (wqCML)
# library(edgeR)
# 
# edger_dispersion <- rep(0, 62)
# names(edger_dispersion)[1:36] <- gtex_expmat_names
# names(edger_dispersion)[37:62] <- tcga_expmat_names
# for(gtex_name in gtex_expmat_names) {
#   edger_dispersion[gtex_name] <- estimateCommonDisp(get(gtex_name))
#   gc()
# }
# for(tcga_name in tcga_expmat_names) {
#   edger_dispersion[tcga_name] <- estimateCommonDisp(get(tcga_name))
#   gc()
# }
# 
# barplot(edger_dispersion, las=2)
