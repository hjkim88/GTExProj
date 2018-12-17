###
#   File name : MakeRandomNetworkMap.R
#   Author    : Hyunjin Kim
#   Date      : Oct 11, 2018
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Make a gene map for 1000 random network
#
#   Instruction
#               1. Source("MakeRandomNetworkMap.R")
#               2. Run the function "makeRandomNetMap()" - specify the input (64 Aracne Networks) and output file Path
#               3. A map will be generated as the output file path
#
#   Example
#               > source("The_directory_of_MakeRandomNetworkMap.R/MakeRandomNetworkMap.R")
#               > makeRandomNetMap(aracnePath="./All_64_Aracne_MI_Fixed.rda",
#                                  permutation=1000,
#                                  outputPath="./Random_Network_Gene_Map.rda")
###

makeRandomNetMap <- function(aracnePath="./All_64_Aracne_MI_Fixed.rda",
                             permutation=1000,
                             outputPath="./Random_Network_Gene_Map.rda") {
  
  ### load aracne networks
  load(aracnePath)
  
  
  ### get all genes from the networks
  ### needs "aracne.R" loaded
  all_genes <- getInteractomeGenes(varNames, count = FALSE)
  
  
  ### make an empty map
  randomNetMap <- matrix(NA, length(all_genes), permutation)
  rownames(randomNetMap) <- all_genes
  colnames(randomNetMap) <- paste0("Random", 1:permutation)
  
  
  ### put random values
  set.seed(1234)
  for(i in 1:permutation) {
    randomNetMap[,i] <- rownames(randomNetMap)[sample(length(all_genes), length(all_genes))]
  }
  
  
  # ### testing
  # temp <- apply(randomNetMap, 2, function(x) {
  #           return(length(unique(x)))
  #         })
  # print(length(which(temp != length(all_genes))))
  
  
  ### set README function
  README <- function(){
    writeLines(paste(rep("#", 100), collapse = ""))
    writeLines("The \"randomNetMap\" has mapping info between original genes and new random genes")
    writeLines("The rows are genes which is an unique set of all the existing genes in all the networks")
    writeLines("The columns represent 1000 permutation")
    writeLines("The rownames(randomNetMap) has the original Entrez IDs")
    writeLines(paste(rep("#", 100), collapse = ""))
  }
  
  
  ### save the map as a RDA file
  save(list = c("randomNetMap", "README"), file = outputPath)
  
}

