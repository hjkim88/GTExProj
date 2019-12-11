### This is a test script for sample run of network similarity computation
### AdiposeVis

### load the igraphs
load("C:/Research/CUMC/GTExProj/data/RDA_Files/All_62_Aracne_igraphs.rda")

### load library
if(!require(graphkernels, quietly = TRUE)) {
  install.packages("graphkernels")
  require(graphkernels, quietly = TRUE)
}
if(!require(igraph, quietly = TRUE)) {
  install.packages("igraph")
  require(igraph, quietly = TRUE)
}

### only retain the igraphs of GTEx tissues
igs <- igs[-c(37:62)]
gc()

### set target network - AdiposeVis
target_net <- igs[["AdiposeVis"]]

### 1
start <- 1
end <- 9
result <- vector("list", length = end - start + 1)
names(result) <- names(igs)[seq(start, end)]
for(tissue in names(result)) {
  result[[tissue]] <- CalculateGeometricRandomWalkKernel(G = list(target_net, igs[[tissue]]), par = 0.1)
}
result1 <- result
save(list = c("result1"), file = "C:/Research/CUMC/GTExProj/result1.rda")

### 2
start <- 10
end <- 18
result <- vector("list", length = end - start + 1)
names(result) <- names(igs)[seq(start, end)]
for(tissue in names(result)) {
  result[[tissue]] <- CalculateGeometricRandomWalkKernel(G = list(target_net, igs[[tissue]]), par = 0.1)
}
result2 <- result
save(list = c("result2"), file = "C:/Research/CUMC/GTExProj/result2.rda")

### 3
start <- 19
end <- 27
result <- vector("list", length = end - start + 1)
names(result) <- names(igs)[seq(start, end)]
for(tissue in names(result)) {
  result[[tissue]] <- CalculateGeometricRandomWalkKernel(G = list(target_net, igs[[tissue]]), par = 0.1)
}
result3 <- result
save(list = c("result3"), file = "C:/Research/CUMC/GTExProj/result3.rda")

### 4
start <- 28
end <- 36
result <- vector("list", length = end - start + 1)
names(result) <- names(igs)[seq(start, end)]
for(tissue in names(result)) {
  result[[tissue]] <- CalculateGeometricRandomWalkKernel(G = list(target_net, igs[[tissue]]), par = 0.1)
}
result4 <- result
save(list = c("result4"), file = "C:/Research/CUMC/GTExProj/result4.rda")


### combine the result RDA files
load("C:/Research/CUMC/GTExProj/result1.rda")
load("C:/Research/CUMC/GTExProj/result2.rda")
load("C:/Research/CUMC/GTExProj/result3.rda")
load("C:/Research/CUMC/GTExProj/result4.rda")
result <- c(sapply(result1, function(x) x[1,2]),
            sapply(result2, function(x) x[1,2]),
            sapply(result3, function(x) x[1,2]),
            sapply(result4, function(x) x[1,2]))
minus_log2_result <- -log2(result)

### write the result
write.table(data.frame(Tissue=names(result), Similarity=result),
            file = "C:/Research/CUMC/GTExProj/AdiposeVis_RandomWalk_Similarity_Among_GTEx.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

### load library
if(!require(ggrepel, quietly = TRUE)) {
  install.packages("ggrepel")
  library(ggrepel, quietly = TRUE)
}

### make an 1D similarity plot
df <- data.frame(x=minus_log2_result, y=0)
color <- rep("blue", nrow(df))
names(color) <- rownames(df)
color["AdiposeVis"] <- "red"
ggplot(data = df, aes(x=x, y=y)) +
  geom_point(color = "black", size = 1) +
  geom_label_repel(aes(x, y, label = rownames(df)), color = color, box.padding = unit(0.45, "lines")) +
  labs(title="Similarity to AdiposeVis") +
  xlab("-log2(Similarity)") +
  theme_classic(base_size = 16) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        axis.line.y=element_blank(),
        axis.text.y=element_blank())

### save the plot
ggsave(filename = "C:/Research/CUMC/GTExProj/AdiposeVis_RandomWalk_Similarity_Among_GTEx.png", width = 20, height = 10)
