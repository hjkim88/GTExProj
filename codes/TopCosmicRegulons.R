###
#   File name : TopCosmicRegulons.R
#   Author    : Hyunjin Kim
#   Date      : Nov 11, 2019
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Compute p-values of Cosmic enrichment in all the regulons in all the tissues
#
#   Instruction
#               1. Source("TopCosmicRegulons.R")
#               2. Run the function "mutated_genes" - specify the necessary input directory and output directory
#               3. The result RDA file will be generated in the output path
#
#   Example
#               > source("The_directory_of_TopCosmicRegulons.R/TopCosmicRegulons.R")
#               > mutated_genes(aracneRDAPath="C:/Research/CUMC/GTExProj/data/RDA_Files/All_62_ARACNE.rda",
#                               cosmicFilePath="C:/Research/CUMC/GTExProj/data/Cosmic/Cosmic_Census_100419_all.tsv",
#                               outputPath="C:/Research/CUMC/GTExProj/data/RDA_Files/All_62_Top_Cosmic_Regulons.rda")
###

mutated_genes <- function(aracneRDAPath="C:/Research/CUMC/GTExProj/data/RDA_Files/All_62_ARACNE.rda",
                          cosmicFilePath="C:/Research/CUMC/GTExProj/data/Cosmic/Cosmic_Census_100419_all.tsv",
                          outputPath="C:/Research/CUMC/GTExProj/data/RDA_Files/All_62_Top_Cosmic_Regulons.rda") {
  
  ### load data
  load(aracneRDAPath)
  cgc <- read.table(file = cosmicFilePath, header = TRUE, sep = "\t",
                    stringsAsFactors = FALSE, check.names = FALSE)
  
  ### make an empty list for saving top cosmic enriched regulons
  cosmic_enriched_regulons <- vector("list", length = length(varNames))
  names(cosmic_enriched_regulons) <- varNames
  
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
  
  ### A function to correct p-values with Benjamini-Hochberg approach
  correct_bh <- function(pvs) {
    
    temp <- cbind(p=pvs, No=seq(1:length(pvs)))
    
    if(length(which(is.na(temp[,"p"]))) > 0) {
      temp[which(is.na(temp[,"p"])),"p"] <- 1
    }
    
    temp <- temp[order(temp[,"p"]),]
    temp <- cbind(temp, Rank=seq(1:length(pvs)), BH=1)
    
    
    temp[length(pvs), "BH"] <- temp[length(pvs), "p"]
    for(i in (length(pvs)-1):1) {
      temp[i,"BH"] <- min(temp[i+1, "BH"], temp[i,"p"]*length(pvs)/temp[i,"Rank"])
    }
    
    temp <- temp[order(temp[,"No"]),]
    
    return(as.numeric(temp[,"BH"]))
  }
  
  ### perform analysis for all the tissues
  for(aracne_name in varNames) {
    ### get Aracne network
    aracne <- get(aracne_name)
    
    ### get total genes in the network
    interactome_total_genes <- as.character(getInteractomeGenes(aracne_name, count = FALSE))
    
    ### get cancer genes in the network
    interactome_cancer_genes <- intersect(interactome_total_genes, as.character(cgc$`Entrez GeneId`))
    
    ### compute cosmic enrichment (Fisher's exact test p-values) for all the hubs in the TCGA Aracne network
    all_hubs <- rownames(aracne[[1]])
    all_enrichment_pvs <- rep(0, length(all_hubs))
    names(all_enrichment_pvs) <- all_hubs
    for(hub in all_hubs) {
      ### get target genes
      target_genes <- rownames(aracne[[2]][[hub]])
      
      ### compute enriched genes
      enriched_genes <- intersect(target_genes, interactome_cancer_genes)
      
      ### calculate p-value
      ### Fisher's exact test
      ###
      ###                 regulon   no-regulon
      ###               -----------------------
      ### cancer gene   |   X           Y
      ### no-cancer gene|   Z           W
      x <- length(enriched_genes)
      y <- length(interactome_cancer_genes) - x
      z <- length(target_genes) - x
      w <- length(interactome_total_genes) - x - y - z
      
      ### Fisher's exact test p-value
      all_enrichment_pvs[hub] <- fisher.test(matrix(c(x, z, y, w), 2, 2), alternative = "greater")$p.value
    }
    
    ### get top cosmic enriched hubs
    all_enrichment_pvs <- all_enrichment_pvs[order(all_enrichment_pvs)]
    
    ### save the cosmic enriched regulons
    cosmic_enriched_regulons[[aracne_name]] <- data.frame(PVal=all_enrichment_pvs, FDR=correct_bh(all_enrichment_pvs),
                                                          stringsAsFactors = FALSE, check.names = FALSE)
  }
  
  ### set README function
  README <- function(){
    writeLines(paste(rep("#", 100), collapse = ""))
    writeLines("The \"cosmic_enriched_regulons\" is a list that contains cosmic-enriched regulon info")
    writeLines("It has 62 entries which indiciates that there are cosmic-enriched regulon info for 62 GTEx+TCGA tissues.")
    writeLines("In each of the 62 entries, there is a data frame of Fisher's Exact Test p-values.")
    writeLines("rownames(cosmic_enriched_regulons[[\"TISSUE_NAME\"]]) represents the cosmic-enriched hubs.")
    writeLines("")
    writeLines("PVal: p-values from Fisher's Exact Test of Cosmic enrichment test")
    writeLines("FDR: corrected p-values using Benjamini-Hochberg")
    writeLines("")
    writeLines("It was generated as the following:")
    writeLines("For every TCGA tissue:")
    writeLines("\tFor every hub H:")
    writeLines("\t\tCalculate p-value of Cosmic enrichment using Fisher's Exact Test.")
    writeLines("\tCorrect the p-values with BH method.")
    writeLines(paste(rep("#", 100), collapse = ""))
  }
  
  ### save as RDA
  save(list = c("cosmic_enriched_regulons", "README"), file = outputPath)
  
}
