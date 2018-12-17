###
#   File name : MakeFetRDA.R
#   Author    : Hyunjin Kim
#   Date      : Jul 13, 2018
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Make a hub FET matrix for each tissue and save them in a RDA file
#
#   Instruction
#               1. Source("MakeFetRDA.R")
#               2. Run the function "fetRDA()" - specify the input aracne RDA and output path
#               3. A RDA object of hub FET matrices will be generated in the output path
#
#   Example
#               > source("The_directory_of_MakeFetRDA.R/MakeFetRDA.R")
#               > fetRDA(aracneFile="./All_64_Aracne_MI_Fixed.rda",
#                        outputPath="./all_64_regulon_fet.rda")
###

fetRDA <- function(aracneFile="./All_64_Aracne_MI_Fixed.rda",
                   outputPath="./all_64_regulon_fet.rda") {
  
  ### load libraries
  if(!require(BiocParallel)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("BiocParallel")
    library(BiocParallel)
  }
  
  
  ### load input
  load(aracneFile)
  
  
  ### a function to make hub fet matrix
  makeFetMat <- function(net) {
    # *****************************************************************************
    # Return the genes (both hubs and targets; or just hubs) in all interactomes 
    # passed as argument. A gene can be returned if it appears in at least one
    # interactome; or if it appears in all interactomes.
    #
    # ARGUMENTS:
    # nets:	One of the following:
    #		- string vector (e.g., "Spleen", or c("Spleen", "tcga_gbm")) containing
    #		variable names of interactome objects
    #		to an ARACNe network.
    #		- a single interactome variable.
    #		- a list of interactome variables, e.g., list(Spleen, tcga_gbm)
    # count:	If TRUE, return the number of genes in the interactomes. Otherwise
    #		return a vector comprising the Entrez IDs of the genes.
    # hubs_only:	if TRUE, return only hub genes. Otherwise, return all genes,
    #		both hubs and not.
    # common:	If TRUE, return only genes that appear in all interactomes.
    #		Otherwise return genes that appear in at least one interactome.
    #
    # RETURN VALUE:
    # It generates a vector containing all the unique genes encountered as hubs or 
    # targets in all interactomes specified in the arguments "nets". if count=TRUE
    # it ruturns the slze of that set. If count=FALSE it returns the full vector
    # of the genes (in the form of integer Entrez ids). If "hubs_only" is TRUE it 
    # only considers hub genes. If "common" is TRUE it only considers genes that 
    # are seeon in at least one interactions in every interactome in "nets"
    # *****************************************************************************	
    getInteractomeGenes <- function(nets, count = TRUE, hubs_only = FALSE, common = FALSE){
      netList = list()
      if (is.character(nets)){
        for (i in 1:length(nets))
          netList[[i]] = get(nets[i])
      } else if (is.list(nets) && !is.matrix(nets[[1]]))
        netList = nets
      else
        netList[[1]] = nets
      allGenes = list()
      for (i in 1:length(netList)){
        net = netList[[i]][[2]]
        genes = NULL
        if (!hubs_only)
          genes = unique(unlist(sapply(net, function(r){
            return(abs(r[,1]))
          })))
        genes = unique(c(genes, as.integer(names(net))))
        allGenes[[i]] = genes
      }
      
      if (length(allGenes) == 1)
        genes = allGenes[[1]]
      else{
        if (!common)
          genes = unique(Reduce(union, allGenes))
        else
          genes = unique(Reduce(intersect, allGenes))
      }
      if (count)
        return(length(genes))
      return(genes)
    }
    
    # ******************************************************************************************
    # Fisher's exact test for the intesection of two sets.
    #
    # ARGUMENTS:
    # * s1:		The first of the two sets, in the form of a vector of objects (objects are expected
    #		to be comparable through the == operator). E.g., in the case of a regulons this could 
    #		be a vector of the entrez ids of the regulon target genes.
    # * s2:		The second of the two sets, in the form of a vector of objects.
    # * total:	The total number of objects in the universe the sets are drawn from. E.g., in the
    #		case where s1 and s2 represent regulons from the same interactome this could be the
    #		total number of genes in the interactome.
    # * alternative:	Same as the argment 'alternative' of the method fisher.test(). The default
    #		value is "greater", indicating that we only care where the intersection of s1 and s2
    #		is greater than expected by chance. In this case, an intersection whose size is much
    #		lower than what expected by chance will be treated as non-signficant.
    # ******************************************************************************************
    fet.set <- function(s1, s2, total, alternative = "greater"){
      common = length(intersect(s1, s2))
      s1_minus_s2 = max(length(s1), common) - common
      s2_minus_s1 = max(length(s2), common) - common
      remainder = total - (common + s1_minus_s2 + s2_minus_s1)
      return(fisher.test(rbind(c(common, s1_minus_s2), c(s2_minus_s1, remainder)), alternative = alternative))
    }
    
    ### get number of hubs
    N <- nrow(net[[1]])
    
    ### get number of genes
    all_num <- getInteractomeGenes(net) 
    
    ### get hub gene names
    G <- rownames(net[[1]])
    
    ### make an empty matrix
    M <- matrix(NA, N, N)
    colnames(M) <- G
    rownames(M) <- G
    
    ### calculate hub FETs
    for(j in 1:(N-1)) {
      for(k in (j+1):N) {
        fet <- fet.set(rownames(net[[2]][[net[[1]][G[j],2]]]), rownames(net[[2]][[net[[1]][G[k],2]]]), total=all_num)
        M[j, k] <- M[k, j] <- log10(fet$p)
      }
    }
    
    ### -Inf to diagonal
    diag(M) <- -Inf
    
    return(M)
  }
  
  
  ### make a one big list
  input <- list()
  for(i in 1:length(varNames)) {
    input[[i]] <- get(varNames[i])
  }
  
  
  ### get FET matrices
  set.seed(1234)
  results <- bplapply(input, makeFetMat)
  
  
  ### make fetMatNames and result matrices
  fetMatNames <- NULL
  for(i in 1:length(varNames)) {
    ### add FET matrix variable name
    if(grepl("tcga", varNames[i])) {
      fetMatNames <- c(fetMatNames, paste0("fet_", varNames[i]))
    } else {
      fetMatNames <- c(fetMatNames, paste0("fet_gtex_", varNames[i]))
    }
    
    assign(fetMatNames[i], results[[i]])
  }
  
  
  ### README function
  README <- function() {
    writeLines("There are 64 matrices and each matrix contains FET scores of hubs x hubs of the corresponding tissue.")
    writeLines("Each matrix has same number of rows and columns: the number of hubs in the tissue.")
    writeLines("The values in a matrix indicate log10(p-value of FET).")
    writeLines("The FET basically means how close the two hubs are.")
    writeLines("The matrices' name should be either fet_gtex_<TissueAcronym> or fet_tcga_<TumorAcronym>")
    writeLines("fetMatNames has all the names of the matrices")
  }
  
  
  ### save the result
  save(list = c("fetMatNames", fetMatNames, "README"), file = outputPath)
  
}
