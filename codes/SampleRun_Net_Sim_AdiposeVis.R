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


### This is another sample run for getting a sense of impact of regulon size on the similarity

### load Aracne network
load("C:/Research/CUMC/GTExProj/data/RDA_Files/All_62_ARACNE.rda")

### get the lung cacer network
aracne_net <- Lung[[2]]

### remove all the networks
rm(list = c(varNames, "varNames", "tfPairEnrich", "tfPairProb", "netSizes", "pairWise"))
gc()

set.seed(1234)
### network with 3 hubs
aracne_net3 <- aracne_net[sample(length(aracne_net), 3)]

### network with 5 hubs
aracne_net5 <- aracne_net[sample(length(aracne_net), 5)]

### network with 10 hubs
aracne_net10 <- aracne_net[sample(length(aracne_net), 10)]

### network with 50 hubs
aracne_net50 <- aracne_net[sample(length(aracne_net), 50)]

### network with 100 hubs
aracne_net100 <- aracne_net[sample(length(aracne_net), 100)]

### a function to
transform_to_igraph <- function(aracne_network) {
  ### gather all interactions for the given network
  edge_list <- data.frame(Hub=names(aracne_network)[1],
                          Target=as.character(abs(aracne_network[[1]][,1])),
                          weight=aracne_network[[1]][,2],
                          stringsAsFactors = FALSE)
  for(i in 2:length(aracne_network)) {
    edge_list <- rbind(edge_list,
                       data.frame(Hub=names(aracne_network)[i],
                                  Target=as.character(abs(aracne_network[[i]][,1])),
                                  weight=aracne_network[[i]][,2],
                                  stringsAsFactors = FALSE))
  }
  
  ### transform into an igraph
  result_igs <- graph.data.frame(edge_list, directed = FALSE)
  
  ### simplify the igraph (remove duplicate edges)
  result_igs <- simplify(result_igs,
                         remove.multiple = TRUE,
                         remove.loops = TRUE,
                         edge.attr.comb = "first")
  
  return(result_igs)
}

### complete igs
igs1 <- igs[["Lung"]]

### remove igs object
rm(list = c("igs"))
gc()

### igs with 3 hubs
igs3 <- transform_to_igraph(aracne_net3)

### igs with 5 hubs
igs5 <- transform_to_igraph(aracne_net5)

### igs with 10 hubs
igs10 <- transform_to_igraph(aracne_net10)

### igs with 50 hubs
igs50 <- transform_to_igraph(aracne_net50)

### igs with 100 hubs
igs100 <- transform_to_igraph(aracne_net100)

### exp1
exp1 <- CalculateGeometricRandomWalkKernel(G = list(igs1, igs3), par = 0.1)

### exp2
exp2 <- CalculateGeometricRandomWalkKernel(G = list(igs1, igs5), par = 0.1)

### exp3
exp3 <- CalculateGeometricRandomWalkKernel(G = list(igs1, igs10), par = 0.1)

### exp4
exp4 <- CalculateGeometricRandomWalkKernel(G = list(igs1, igs50), par = 0.1)

### exp5
exp5 <- CalculateGeometricRandomWalkKernel(G = list(igs1, igs100), par = 0.1)

### exp6
exp6 <- CalculateGeometricRandomWalkKernel(G = list(igs3, igs100), par = 0.1)

### exp7
exp7 <- CalculateGeometricRandomWalkKernel(G = list(igs1, igs3, igs5, igs10, igs50, igs100), par = 0.1)

### combine all the results
exp_results <- list(exp1, exp2, exp3, exp4, exp5, exp6, exp7)

### save the exp results
save(list = c("exp_results"), file = "C:/Research/CUMC/GTExProj/data/RDA_Files/exp_results.rda")
