# *****************************************************************************
# This file contains functions related to Go term analysis on regulons
# *****************************************************************************

# ******************** cluster_regulons_with_pathways *****************************
# We ran pathway analysis on every regulon of every tissue of GTEx and TCGA.
# The pathways indicate that they are associated with target genes of the hubs.
# If we cluster the regulons based on their pathways, we may get hubs/regulons that have
# similar biolocial functions. And we could also see the biological functional differences
# among tissues or between GTEx and TCGA.
# 
# This function needs [RegulonPathwayAnnotation.rda] file which contains pathway analysis
# results of all the regulons of all the tissues of GTExa and TCGA.
#
# The distance is calculated based on Jaccard distance
# = 1 - Intersection over Union
#
# params[[1]]: The file path of the "RegulonPathwayAnnotation.rda" file
#              (a character vector of length 1)
# params[[2]]: A Jaccard distance threshold for selecting top hub pairs
#              Should bigger than 0 and smaller or equal than 1
#              If NA, then select all the hub pairs
#              (a numeric value between 0-1 or NA)
# params[[3]]: An integer threshold that determines the most appeared pathways
#              e.g., if 50, the top 50 most appeared pathways will be presented
#              if NA, then select all the pathways
#              (an integer value or NA)
# params[[4]]: The output directory that results will be printed out
#              (a character vector of length 1)
#
# e.g., params <- list("//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/RegulonPathwayAnnotation.rda",
#                      0.1, 50,
#                      "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/regulon_pathway/")
# e.g., params <- list("./data/RDA_Files/RegulonPathwayAnnotation.rda", 0.1, 50, "./results/regulon_pathway/")

cluster_regulons_with_pathways <- function(params) {
  
  ### argument checking
  assertString(params[[1]])
  assertNumeric(params[[2]])
  assertIntegerish(params[[3]])
  assertString(params[[4]])
  
  ### load the pathway analysis results
  load(params[[1]])
  
  ### create an empty matrix for distances among regulons
  distance_mats <- vector("list", length(varGOnames))
  names(distance_mats) <- varGOnames
  
  ### calculate distances among regulons based on their associated pathways for each tissue
  for(tissue in names(distance_mats)) {
    ### write a log to present progress
    logLines(paste("\nProcessing tissue -> ", tissue))
    
    ### get the pathway analysis results for the given tissue
    pathRes <- get(tissue)
    
    ### create an empty distance matrix for the given tissue
    distance_mats[[tissue]] <- matrix(NA, length(pathRes), length(pathRes))
    rownames(distance_mats[[tissue]]) <- names(pathRes)
    colnames(distance_mats[[tissue]]) <- names(pathRes)
    
    ### compuate distances based on the Jaccard distance
    for(i in 1:(nrow(distance_mats[[tissue]])-1)) {
      for(j in (i+1):ncol(distance_mats[[tissue]])) {
        distance_mats[[tissue]][i, j] <- 1 - (length(intersect(pathRes[[i]][,1], pathRes[[j]][,1])) /
                                                length(union(pathRes[[i]][,1], pathRes[[j]][,1])))
      }
    }
    
    ### distance_mats[[tissue]][i,i] = 0
    diag(distance_mats[[tissue]]) <- 0
    
    ### there are some already-calculated stuffs since distance_mats[[tissue]][i,j] == distance_mats[[tissue]][j,i] 
    ### just copy them into the appropriate places
    distance_mats[[tissue]][lower.tri(distance_mats[[tissue]])] <- t(distance_mats[[tissue]])[lower.tri(distance_mats[[tissue]])]
  }
  
  ### set README function
  README <- function(){
    writeLines(paste(rep("#", 100), collapse = ""))
    writeLines("The \"distance_mats\" is a list with 62 length.")
    writeLines("In each element of the list, there is a distance matrix for each tissue from GTEx and TCGA.")
    writeLines("Each matrix represents closeness between two regulons based on shared pathways.")
    writeLines("The distance is calculated by Jaccard distance.")
    writeLines("For example, if two regulons share many associated pathways, their distance value")
    writeLines("would be small, and otherwise, they would be big.")
    writeLines("If you want to see the pathways themselves, please refer: RegulonPathwayAnnotation.rda.")
    writeLines(paste(rep("#", 100), collapse = ""))
  }
  
  ### save the distance matrices in RDA file
  save(list = c("distance_mats", "README"), file = paste0(params[[4]], "All_62_RegulonPathwayDistanceMats.rda"))
  
  ### get hubs that have the lowest distance with the threshold
  top_hubs <- lapply(distance_mats, function(x) {
    ### there are duplicates and useless (same hub) pairs, so filter them out
    x[lower.tri(x)] <- 1
    diag(x) <- 1
    
    ### get array indicies of the targets
    if(is.na(params[[2]])) {
      arrIndices <- which(x <= Inf, arr.ind = TRUE)
    } else {
      arrIndices <- which(x < params[[2]], arr.ind = TRUE)
    }
    
    ### make hub pair info
    if(nrow(arrIndices) > 0) {
      hubPair <- matrix(NA, nrow(arrIndices), 4)
      colnames(hubPair) <- c("Hub1", "Hub2", "Distance", "EmpiricalPVal")
      hubPair[,"Hub1"] <- as.numeric(rownames(x)[arrIndices[,1]])
      hubPair[,"Hub2"] <- as.numeric(colnames(x)[arrIndices[,2]])
      hubPair[,"Distance"] <- x[arrIndices]
      hubPair <- hubPair[order(hubPair[,"Distance"]),]
      totalPairNum <- (nrow(x) * ncol(x) - ncol(x)) / 2
      hubPair[,"EmpiricalPVal"] <- rank(hubPair[,"Distance"], ties.method = "min") / totalPairNum
      return(hubPair)
    } else {
      return(NULL)
    }
  })
  
  ### print the top hubs filtered by the threshold
  for(tissue in names(top_hubs)) {
    dir.create(paste0(params[[4]], tissue), showWarnings = FALSE)
    write.table(top_hubs[[tissue]], file = paste0(params[[4]], tissue, "/top_hub_pairs_", params[[2]], ".txt"), sep = "\t", row.names = FALSE)
  }
  
  ### the most appeared pathways in each tissue
  top_pathways <- vector("list", length(varGOnames))
  names(top_pathways) <- varGOnames
  for(tissue in varGOnames) {
    ### get the pathway analysis results for the given tissue
    pathRes <- get(tissue)
    
    ### create an empty pathway count object
    pathwayCntLen <- length(Reduce(union, lapply(pathRes, function(x) x[,"Description"])))
    pathwayCnt <- vector("integer", pathwayCntLen)
    names(pathwayCnt) <- Reduce(union, lapply(pathRes, function(x) x[,"Description"]))
    
    ### count pathways (how many times each appeared) for each tissue
    for(hub in names(pathRes)) {
      pathwayCnt[pathRes[[hub]][,"Description"]] <- pathwayCnt[pathRes[[hub]][,"Description"]] + 1
    }
    
    ### sort the pathway count in descending order
    pathwayCnt <- pathwayCnt[order(-pathwayCnt)]
    
    ### filter the pathways with the input threshold
    if(!is.na(params[[3]])) {
      pathwayCnt <- pathwayCnt[1:params[[3]]]
    }
    
    ### save the result to the list
    top_pathways[[tissue]] <- pathwayCnt
  }
  
  ### print the top pathways
  for(tissue in names(top_pathways)) {
    dir.create(paste0(params[[4]], tissue), showWarnings = FALSE)
    write.table(data.frame(Pathway=names(top_pathways[[tissue]]), Count=top_pathways[[tissue]]),
                file = paste0(params[[4]], tissue, "/top_pathways_", params[[3]], ".txt"),
                sep = "\t", row.names = FALSE)
  }
  
  ### set README function
  README <- function(){
    writeLines(paste(rep("#", 100), collapse = ""))
    writeLines("The \"top_hubs\" and the \"top_pathways\" are both list with 62 length.")
    writeLines("In each element of the \"top_hubs\", there is a table of hub pairs")
    writeLines(paste0("that have the lowest distance (< ", 0.1, ") values defined by Jaccard index."))
    writeLines("The first two columns indicate a hub pair (Entrez ID),")
    writeLines("the third column has the distance value of a given hub pair,")
    writeLines("and the fourth column means empirical p-values of the corresponding distances.")
    writeLines(paste0("In each element of the \"top_pathways\", there is an integer vector of length ", 50, "."))
    writeLines("They are top pathways which were most appeared among all the regulons in each tissue.")
    writeLines("The first column indicates GO pathway names, and the second column means")
    writeLines("how many times the given pathway was appeared among all the regulons in the tissue.")
    writeLines("If you want to see the pathways themselves, please refer: RegulonPathwayAnnotation.rda.")
    writeLines(paste(rep("#", 100), collapse = ""))
  }
  
  ### save the top_hubs & top_pathways objects
  save(list = c("top_hubs", "top_pathways", "README"),
       file = paste0(params[[4]], "All_62_RegulonPathwayInfo.rda"))
  
}

# ******************** cluster_regulons_with_pathways_heatmap *****************************
# After running cluster_regulons_with_pathways() function,
# we now have distances between every hub pairs of every tissue based on associated pathways.
# In the function, we also produced the top hubs with the designated threshold that is
# based on the closeness of two regulons, and also generated the most shown pathways based on
# their occurrences in each tissue. Now we want to make heatmap plots that which the top
# pathways are enriched with which top hubs in each tissue. 
#
# params[[1]]: The file path of the "All_62_RegulonPathwayInfo.rda" file
#              (a character vector of length 1)
# params[[2]]: The file path of the "RegulonPathwayAnnotation.rda" file
#              (a character vector of length 1)
# params[[3]]: The output directory that results will be printed out
#              (a character vector of length 1)
#
# e.g., params <- list("//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/All_62_RegulonPathwayInfo.rda",
#                      "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/RegulonPathwayAnnotation.rda",
#                      "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/regulon_pathway/")
# e.g., params <- list("./data/RDA_Files/All_62_RegulonPathwayInfo.rda", "./data/RDA_Files/RegulonPathwayAnnotation.rda", "./results/regulon_pathway/")

cluster_regulons_with_pathways_heatmap <- function(params) {
  
  ### argument checking
  assertString(params[[1]])
  assertString(params[[2]])
  assertString(params[[3]])
  
  ### load the necessary data
  load(params[[1]])
  load(params[[2]])
  
  ### get the hubs from the top hups info
  unique_hubs <- lapply(top_hubs, function(x) {
    y <- union(x[,"Hub1"], x[,"Hub2"])
    y <- y[order(as.numeric(y))]
    return(as.character(y))
  })
  
  ### get the pathways from the top pathways info
  unique_pathways <- lapply(top_pathways, function(x) {
    y <- names(x)
    y <- y[order(y)]
    return(y)
  })
  
  ### calculate the heatmap matrix (rows - top pathways, cols - top hubs)
  heatmap_pathway_hub <- vector("list", length(varGOnames))
  names(heatmap_pathway_hub) <- varGOnames
  for(tissue in varGOnames) {
    ### create an empty heatmap matrix
    heatmap_pathway_hub[[tissue]] <- matrix(0, length(unique_pathways[[tissue]]), length(unique_hubs[[tissue]]))
    rownames(heatmap_pathway_hub[[tissue]]) <- unique_pathways[[tissue]]
    colnames(heatmap_pathway_hub[[tissue]]) <- unique_hubs[[tissue]]
    
    ### get the regulon pathway annotation
    GOEnrich <- get(tissue)
    
    ### if the pathway is enriched by the given hub, 1, otherwise, 0
    for(pathway in unique_pathways[[tissue]]) {
      for(hub in unique_hubs[[tissue]]) {
        if(length(grep(pathway, GOEnrich[[hub]][,"Description"])) > 0) {
          heatmap_pathway_hub[[tissue]][pathway, hub] <- 1
        }
      }
    }
    
    ### change the colnames to gene symbols from entrez ids
    colnames(heatmap_pathway_hub[[tissue]]) <- entrezIDtoSymbol(colnames(heatmap_pathway_hub[[tissue]]))
    
    ### print out the heatmap
    png(paste0(params[[3]], tissue, "/heatmap_pathway_hub.png"), width = 2500, height = 1500, res = 120)
    par(oma=c(0,0,0,35))
    heatmap.3(heatmap_pathway_hub[[tissue]], main = paste0(tissue, "_Heatmap_Pathway_Hub"),
              xlab = "", ylab = "", col=c("white", "blue"),
              scale="none", key=F, keysize=0.5, dendrogram = "col", trace = 'none',
              labRow = rownames(heatmap_pathway_hub[[tissue]]), labCol = colnames(heatmap_pathway_hub[[tissue]]),
              Rowv = FALSE, Colv = TRUE,
              cexRow = 1.2, cexCol = 0.8, na.rm = TRUE)
    dev.off()
  }
  
}

# ******************** regulon_conservation_jaccard_vs_fet *****************************
# Distance between two regulons can be computed in several ways.
# We have done it with two ways.
# 1. P-value from Fisher's exact test based on target gene conservation 
# 2. Jaccard distance between associated pathways with the regulons
# In general, if two regulons share similar target genes, then their
# associated pathways should be also similar. Therefore, there would be
# no big deal between the results with the two measurements, but if there
# are any cases that a hub pair have small Jaccard distance but have 
# large FET p-value (or other way around), it would be interesting.
# 
# * You need 22GB memory to run this function.
#
# params[[1]]: The RDA file path of the Jaccard distanceS
#              All_62_RegulonPathwayDistanceMats.rda
#              (a character vector of length 1)
# params[[2]]: The RDA file path of the FET p-values
#              All_62_tfNetEnrich.rda
#              (a character vector of length 1)
# params[[3]]: The RDA file path of the Aracne networks
#              All_62_ARACNE.rda
#              (a character vector of length 1)
# params[[4]]: The RDA file path of the regulon pathway annotations
#              RegulonPathwayAnnotation.rda
#              (a character vector of length 1)
# params[[5]]: A cut-off for selecting top & bottom hub pairs
#              It means percentage, e.g., 10 = top & bottom 10% based on Jaccard or FET
#              It will be used to identify weird hub pairs such as very large Jaccard
#              distance but small FET P-value (or the other way around)
#              (a double between 0-100)
# params[[6]]: A cut-off for selecting hub pairs that have large differences between
#              Jaccard and FET. e.g., if 1000, then selecting top 1000 hub pairs
#              which have absolute large differences between Jaccard and FET.
#              (An integer)
# params[[7]]: The output directory that results will be printed out
#              (a character vector of length 1)
#
# e.g., params = list("//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/All_62_RegulonPathwayDistanceMats.rda",
#                     "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/All_62_tfNetEnrich.rda",
#                     "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/All_62_ARACNE.rda",
#                     "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RegulonPathwayAnnotation.rda",
#                     0.1, 1000,
#                     "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/regulon_pathway/")
# e.g., params = list("./data/RDA_Files/All_62_RegulonPathwayDistanceMats.rda",
#                     "./data/RDA_Files/All_62_tfNetEnrich.rda",
#                     "./data/RDA_Files/All_62_ARACNE.rda",
#                     "./data/RDA_Files/RegulonPathwayAnnotation.rda",
#                     0.1, 1000, "./results/regulon_pathway/")

regulon_conservation_jaccard_vs_fet <- function(params) {
  
  ### argument checking
  assertString(params[[1]])
  assertString(params[[2]])
  assertString(params[[3]])
  assertString(params[[4]])
  assertNumeric(params[[5]])
  assertIntegerish(params[[6]])
  assertString(params[[7]])
  
  ### load library
  if(!require(ggplot2, quietly = TRUE)) {
    install.packages("ggplot2")
    require(ggplot2, quietly = TRUE)
  }
  
  ### load the Jaccard distances and FET p-values
  load(params[[1]])
  load(params[[2]])
  load(params[[3]])
  load(params[[4]])
  
  ### adjust the list names for consistency between distance_mats & tfNetEnrich
  names(distance_mats) <- sapply(names(distance_mats), function(x) substr(x, 1, nchar(x)-2))
  
  for(tissue in names(distance_mats)) {
    ### progress print
    logLines(paste("\nProcessing tissue -> ", tissue))
    
    ### create a matrix for correlation plot
    cor_mat <- matrix(NA, nrow(distance_mats[[tissue]]) * (nrow(distance_mats[[tissue]])-1) / 2, 2)
    colnames(cor_mat) <- c("Jaccard", "FET")
    
    ### fill the empty cells
    cnt <- 1
    temp <- rep("", nrow(cor_mat))
    for(i in 1:(nrow(distance_mats[[tissue]])-1)) {
      hub1 <- rownames(distance_mats[[tissue]])[i]
      for(j in (i+1):ncol(distance_mats[[tissue]])) {
        hub2 <- colnames(distance_mats[[tissue]])[j]
        cor_mat[cnt,"Jaccard"] <- -log10(distance_mats[[tissue]][hub1, hub2])
        cor_mat[cnt,"FET"] <- tfNetEnrich[[tissue]][hub1, hub2]
        temp[cnt] <- paste0(hub1, "_", hub2)
        cnt <- cnt+1
      }
    }
    rownames(cor_mat) <- temp
    
    ### remove duplicates
    cor_mat <- cor_mat[which(!duplicated(cor_mat)),]
    
    ### print a correlation plot
    ggplot(data = as.data.frame(cor_mat), aes(x=Jaccard, y=FET)) +
      geom_point(color = "black", size = 1) +
      labs(title=paste0(tissue, "_Correlation_Jaccard_vs_FET")) +
      xlab("-log10(Jaccard distance)") +
      ylab("-log10(FET P Val)") +
      geom_smooth(method = lm, color="blue", se=FALSE) +
      theme_classic(base_size = 16)
    ggsave(filename = paste0(params[[7]], tissue, "GO/cor_jaccard_vs_fet.png"),
           width = 20, height = 20)
    
    ### change the Inf to a finite number
    cor_mat[which(cor_mat == Inf)] <- 10**(floor(log10(max(cor_mat[is.finite(cor_mat)])))+1)-1
    
    ### select top & bottom hub pairs based on Jaccard & FET
    cut_off <- floor(nrow(cor_mat) * params[[5]] / 100)
    top_jaccard <- rownames(cor_mat)[order(-cor_mat[,"Jaccard"])[1:cut_off]]
    bottom_jaccard <- rownames(cor_mat)[order(cor_mat[,"Jaccard"])[1:cut_off]]
    top_fet <- rownames(cor_mat)[order(-cor_mat[,"FET"])[1:cut_off]]
    bottom_fet <- rownames(cor_mat)[order(cor_mat[,"FET"])[1:cut_off]]
    
    ### get interesting hub pairs
    top_top <- intersect(top_jaccard, top_fet)
    bottom_bottom <- intersect(bottom_jaccard, bottom_fet)
    top_bottom <- intersect(top_jaccard, bottom_fet)
    bottom_top <- intersect(bottom_jaccard, top_fet)
    
    ### write out the results
    write.table(data.frame(paste(paste("TOP_TOP=", paste(top_top, collapse = " ")),
                                 paste("TOP_BOTTOM=", paste(top_bottom, collapse = " ")),
                                 paste("BOTTOM_TOP=", paste(bottom_top, collapse = " ")),
                                 paste("BOTTOM_BOTTOM=", paste(bottom_bottom, collapse = " ")),
                                 sep = "\n")),
                file = paste0(params[[7]], tissue, "GO/hub_pairs_jaccard_vs_fet.txt"),
                sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    ### get absolute difference between Jaccard and FET
    abs_diff <- abs(cor_mat[,"Jaccard"] - cor_mat[,"FET"])
    
    ### sort the abs_diff in descending order
    abs_diff <- abs_diff[order(-abs_diff)]
    
    ### filter with the cut-off
    abs_diff <- abs_diff[1:params[[6]]]
    
    ### additional info generatation
    hub1s <- sapply(names(abs_diff), function(x) strsplit(x, split = "_", fixed = TRUE)[[1]][1], USE.NAMES = FALSE)
    hub2s <- sapply(names(abs_diff), function(x) strsplit(x, split = "_", fixed = TRUE)[[1]][2], USE.NAMES = FALSE)
    pathways1 <- lapply(get(paste0(tissue, "GO"))[hub1s], rownames)
    pathways2 <- lapply(get(paste0(tissue, "GO"))[hub2s], rownames)
    shared_pathways <- vector("list", length(pathways1))
    for(i in 1:length(shared_pathways)) {
      shared_pathways[[i]] <- intersect(pathways1[[i]], pathways2[[i]])
    }
    targets1 <- lapply(get(tissue)[[2]][hub1s], rownames)
    targets2 <- lapply(get(tissue)[[2]][hub2s], rownames)
    shared_targets <- vector("list", length(targets1))
    for(i in 1:length(shared_targets)) {
      shared_targets[[i]] <- intersect(targets1[[i]], targets2[[i]])
    }
    
    ### write out the result
    write.table(data.frame(Hub_Pair=names(abs_diff),
                           Abs_Diff=abs_diff,
                           Jaccard=cor_mat[names(abs_diff),"Jaccard"],
                           FET=cor_mat[names(abs_diff),"FET"],
                           Hub1_Pathway_Num=sapply(get(paste0(tissue, "GO"))[hub1s], nrow),
                           Hub2_Pathway_Num=sapply(get(paste0(tissue, "GO"))[hub2s], nrow),
                           Shared_Pathway_Num=sapply(shared_pathways, length),
                           Hub1_Target_Num=get(tissue)[[1]][hub1s,3],
                           Hub2_Target_Num=get(tissue)[[1]][hub2s,3],
                           Shared_Targets=sapply(shared_targets, length)),
                file = paste0(params[[7]], tissue, "GO/abs_diff_jaccard_vs_fet.txt"),
                sep = "\t", row.names = FALSE)
    
    ### garbage collection
    rm(list = tissue)
    rm(list = paste0(tissue, "GO"))
    tfNetEnrich <- tfNetEnrich[-which(names(tfNetEnrich) == tissue)]
    distance_mats <- distance_mats[-which(names(distance_mats) == tissue)]
    gc()
  }
  
}

# ******************** count_regulon_pathways *****************************
# Before, We ran pathway analysis on every regulon of all the interactomes.
# Now we want to know interesting pathways that appeared in only one interactome.
# Therefore, we counted how many times they appeared in GTEx and in TCGA respectively.
#
# params[[1]]: The RDA file path of the regulon pathway annotations
#              RegulonPathwayAnnotation.rda
#              (a character vector of length 1)
#              
# params[[2]]: The number that indicate a separation between GTEx and TCGA for varGOnames
#              e.g., In varGOnames, 1:36 are GTEx tissues and 37:62 are TCGA tissues,
#                    then params[[2]] should be 36
#
# params[[3]]: The output directory that results will be printed out
#              (a character vector of length 1)
#
# e.g., params <- list("//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/RegulonPathwayAnnotation.rda",
#                      36, "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/regulon_pathway/")
# e.g., params <- list("./data/RDA_Files/RegulonPathwayAnnotation.rda", 36, "./results/regulon_pathway/")

count_regulon_pathways <- function(params) {
  
  ### argument checking
  assertString(params[[1]])
  assertIntegerish(params[[2]])
  assertString(params[[3]])
  
  ### load the pathway analysis results
  load(params[[1]])
  
  ### create a list of vectors of the union set of all the pathways for each tissue
  gtex_union_pathway <- lapply(varGOnames[1:params[[2]]], function(x) {
    pathRes <- get(x)
    return(Reduce(union, lapply(pathRes, function(y) y[,"ID"])))
  })
  names(gtex_union_pathway) <- sapply(varGOnames[1:params[[2]]], function(x) substr(x, 1, nchar(x)-2))
  tcga_union_pathway <- lapply(varGOnames[(params[[2]]+1):length(varGOnames)], function(x) {
    pathRes <- get(x)
    return(Reduce(union, lapply(pathRes, function(y) y[,"ID"])))
  })
  names(tcga_union_pathway) <- sapply(varGOnames[(params[[2]]+1):length(varGOnames)], function(x) substr(x, 1, nchar(x)-2))
  
  ### get the total pathway names in GTEx and in TCGA
  gtex_total_pathways <- Reduce(union, gtex_union_pathway)
  tcga_total_pathways <- Reduce(union, tcga_union_pathway)
  
  ### create empty matrices for pathway counts
  gtex_pathway_cnt <- matrix("", length(gtex_total_pathways), 4)
  tcga_pathway_cnt <- matrix("", length(tcga_total_pathways), 4)
  rownames(gtex_pathway_cnt) <- gtex_total_pathways
  rownames(tcga_pathway_cnt) <- tcga_total_pathways
  colnames(gtex_pathway_cnt) <- c("GOID", "Pathway", "Counts", "Tissue")
  colnames(tcga_pathway_cnt) <- c("GOID", "Pathway", "Counts", "Tissue")
  
  ### fill out the matrices
  gtex_pathway_cnt[,"GOID"] <- gtex_total_pathways
  tcga_pathway_cnt[,"GOID"] <- tcga_total_pathways
  gtex_pathway_cnt[,"Pathway"] <- goTermMap[gtex_pathway_cnt[,"GOID"]]
  tcga_pathway_cnt[,"Pathway"] <- goTermMap[tcga_pathway_cnt[,"GOID"]]
  for(i in 1:nrow(gtex_pathway_cnt)) {
    temp <- sapply(gtex_union_pathway, function(x) length(which(x == gtex_pathway_cnt[i,"GOID"])))
    gtex_pathway_cnt[i,"Counts"] <- sum(temp)
    gtex_pathway_cnt[i,"Tissue"] <- paste(names(temp[which(temp == 1)]), collapse = "/")
  }
  for(i in 1:nrow(tcga_pathway_cnt)) {
    temp <- sapply(tcga_union_pathway, function(x) length(which(x == tcga_pathway_cnt[i,"GOID"])))
    tcga_pathway_cnt[i,"Counts"] <- sum(temp)
    tcga_pathway_cnt[i,"Tissue"] <- paste(names(temp[which(temp == 1)]), collapse = "/")
  }
  
  ### order the matrices based on the counts in ascending order
  gtex_pathway_cnt <- gtex_pathway_cnt[order(as.numeric(gtex_pathway_cnt[,"Counts"])),]
  tcga_pathway_cnt <- tcga_pathway_cnt[order(as.numeric(tcga_pathway_cnt[,"Counts"])),]
  
  ### for the writing out, numerize the count column
  gtex_pathway_cnt <- data.frame(GOID=gtex_pathway_cnt[,"GOID"],
                                 Pathway=gtex_pathway_cnt[,"Pathway"],
                                 Counts=as.numeric(gtex_pathway_cnt[,"Counts"]),
                                 Tissue=gtex_pathway_cnt[,"Tissue"])
  tcga_pathway_cnt <- data.frame(GOID=tcga_pathway_cnt[,"GOID"],
                                 Pathway=tcga_pathway_cnt[,"Pathway"],
                                 Counts=as.numeric(tcga_pathway_cnt[,"Counts"]),
                                 Tissue=tcga_pathway_cnt[,"Tissue"])
  
  ### write out the matrices
  write.table(gtex_pathway_cnt, file = paste0(params[[3]], "gtex_pathway_counts.txt"),
              sep = "\t", row.names = FALSE)
  write.table(tcga_pathway_cnt, file = paste0(params[[3]], "tcga_pathway_counts.txt"),
              sep = "\t", row.names = FALSE)
  
}

# ******************** interesting_regulon_pathways *****************************
# We ran pathway analysis on every regulon of every tissue of GTEx and TCGA.
# The pathways indicate that they are associated with target genes of the hubs.
# Now we want to compare GTEx and TCGA based on the pathways that which pathways
# were exclusively showed up in GTEx only or in TCGA only.
# And it would be very interesting to also know if the target genes found in those
# pathways are differentially expressed between cancer and normal.
# Additionally, we should also check if the hubs carrying those pathways are 
# differentially activated in terms of Viper activity.
# 
# This function needs [RegulonPathwayAnnotation.rda] file which contains pathway analysis
# results of all the regulons of all the tissues of GTExa and TCGA, and [GTEx_TCGA_Map.rda]
# file that has the same tissue mapping info between GTEx and TCGA. Additionally,
# [All_62_raw_counts] file and [All_62_ViperMats.rda] file are needed for the further analyses.
#
# params[[1]]: The file path of the "RegulonPathwayAnnotation.rda" file
#              (a character vector of length 1)
# params[[2]]: The file path of the "GTEx_TCGA_Map.rda" file
#              (a character vector of length 1)
# params[[3]]: The file path of the "All_62_raw_counts" file
#              (a character vector of length 1)
# params[[4]]: The file path of the "All_62_ViperMats.rda" file
#              (a character vector of length 1)
# params[[5]]: The file path of the "All_62_Aracne_hubs.rda" file
#              (a character vector of length 1)
# params[[6]]: A cut-off of selecting the top pathways that have different
#              counts between cancer and normal samples
#              if 0.05, then it selects pathways that have regulons enriched with
#              the given pathway with the top 5% p-values.
#              First, the function seeks for pathways that are exclusively appeared
#              only in GTEx or in TCGA, then filter them further with this cut-off.
#              The function uses all the regulons in the tissue (both GTEx and TCGA)
#              to get the total p-value vector and use the 5% lowest p-value as
#              a threshold to filter the exclusive pathways.
#              (a number between 0 and 1)
# params[[7]]: The directory path for the results
#              (a character vector of length 1)
#
# e.g., params = list("//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/RegulonPathwayAnnotation.rda",
#                     "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/GTEx_TCGA_Map.rda",
#                     "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/All_62_raw_counts.rda",
#                     "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/All_62_ViperMats.rda",
#                     "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/All_62_Aracne_hubs.rda",
#                     0.05,
#                     "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/regulon_pathway/GTEx_vs_TCGA/")
# e.g., params = list("./data/RDA_Files/RegulonPathwayAnnotation.rda",
#                     "./data/RDA_Files/GTEx_TCGA_Map.rda",
#                     "./data/RDA_Files/All_62_raw_counts.rda",
#                     "./data/RDA_Files/All_62_ViperMats.rda",
#                     "./data/RDA_Files/All_62_Aracne_hubs.rda",
#                     0.05,
#                     "./results/regulon_pathway/GTEx_vs_TCGA/")

interesting_regulon_pathways <- function(params) {
  
  ### load library for a beeswarm plot
  if(!require(ggplot2, quietly = TRUE)) {
    install.packages("ggplot2")
    require(ggplot2, quietly = TRUE)
  }
  if(!require(ggbeeswarm, quietly = TRUE)) {
    install.packages("ggbeeswarm")
    require(ggbeeswarm, quietly = TRUE)
  }
  if(!require(ggpubr, quietly = TRUE)) {
    install.packages("ggpubr")
    require(ggpubr, quietly = TRUE)
  }
  if(!require("fgsea", quietly = TRUE)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("fgsea")
    require("fgsea", quietly = TRUE)
  }
  
  ### argument checking
  assertString(params[[1]])
  assertString(params[[2]])
  assertString(params[[3]])
  assertString(params[[4]])
  assertString(params[[5]])
  assertNumeric(params[[6]])
  assertString(params[[7]])
  if(is.na(params[[6]]) || params[[6]] <= 0 || params[[6]] > 1) {
    stop("params[[6]] should be between 0 and 1")
  }
  
  ### load the data
  load(params[[1]])
  load(params[[2]])
  load(params[[3]])
  load(params[[4]])
  load(params[[5]])
  
  ### create an empty list for saving pathway count info
  regulon_pathway_count_info <- vector("list", length = nrow(GTEx_TCGA_Map))
  
  ### for every comparison between GTEx and TCGA, find interesting pathways
  for(i in 1:nrow(GTEx_TCGA_Map)) {
    
    ### get regulon pathways
    gtex_regulon_pathways <- get(paste0(GTEx_TCGA_Map[i,"GTEx"], "GO"))
    tcga_regulon_pathways <- get(paste0("tcga_", GTEx_TCGA_Map[i,"TCGA"], "GO"))
    
    ### get gene expressions
    gtex_gexp <- get(paste0("rcntmat_", GTEx_TCGA_Map[i,"GTEx"]))
    tcga_gexp <- get(paste0("rcntmat_tcga_", GTEx_TCGA_Map[i,"TCGA"]))
    
    ### get viper activities
    gtex_viper <- get(paste0("vmat_gtex_", GTEx_TCGA_Map[i,"GTEx"]))
    tcga_viper <- get(paste0("vmat_tcga_", GTEx_TCGA_Map[i,"TCGA"]))
    
    ### get aracne hubs
    gtex_aracne_hubs <- Aracne_hubs[[GTEx_TCGA_Map[i,"GTEx"]]]
    tcga_aracne_hubs <- Aracne_hubs[[paste0("tcga_", GTEx_TCGA_Map[i,"TCGA"])]]
    
    ### create a sub directory for the results
    subdirPath <- paste0("GTEx_", GTEx_TCGA_Map[i,"GTEx"], "_vs_TCGA_", toupper(GTEx_TCGA_Map[i,"TCGA"]))
    dir.create(file.path(params[[7]], subdirPath), showWarnings = FALSE)
    
    ### get union pathways
    gtex_pathways <- Reduce(union, lapply(gtex_regulon_pathways, function(x) x[,"ID"]))
    tcga_pathways <- Reduce(union, lapply(tcga_regulon_pathways, function(x) x[,"ID"]))
    
    ### create empty pathway count matrix
    pathway_cnt <- matrix(0, length(union(gtex_pathways, tcga_pathways)), 3)
    rownames(pathway_cnt) <- union(gtex_pathways, tcga_pathways)
    colnames(pathway_cnt) <- c("GTEx", "TCGA", "Abs(Diff)")
    
    ### count the pathway occurrence
    for(j in 1:length(gtex_regulon_pathways)) {
      for(k in 1:nrow(gtex_regulon_pathways[[j]])) {
        pathway_cnt[gtex_regulon_pathways[[j]][k,"ID"], "GTEx"] <- pathway_cnt[gtex_regulon_pathways[[j]][k,"ID"], "GTEx"] + 1
      }
    }
    for(j in 1:length(tcga_regulon_pathways)) {
      for(k in 1:nrow(tcga_regulon_pathways[[j]])) {
        pathway_cnt[tcga_regulon_pathways[[j]][k,"ID"], "TCGA"] <- pathway_cnt[tcga_regulon_pathways[[j]][k,"ID"], "TCGA"] + 1
      }
    }
    
    ### calculate the absolute differences
    pathway_cnt[,"Abs(Diff)"] <- abs(pathway_cnt[,"GTEx"] - pathway_cnt[,"TCGA"])
    
    ### order the pathway count matrix based on the absolute difference in descending order
    pathway_cnt <- pathway_cnt[order(-pathway_cnt[,"Abs(Diff)"]),]
    
    ### extract the interesting pathways
    interesting_pathways <- rownames(pathway_cnt[union(which(pathway_cnt[,"GTEx"] == 0), which(pathway_cnt[,"TCGA"] == 0)),])
    
    ### get the counts of the interesting pathways
    interesting_pathway_cnt <- data.frame(pathway_cnt[interesting_pathways,], Hubs="",
                                          stringsAsFactors = FALSE, check.names = FALSE)
    
    ### determine the threshold p-value
    total_pvs <- unlist(lapply(c(gtex_regulon_pathways, tcga_regulon_pathways), function(x) x[,"p.adjust"]))
    total_pvs <- total_pvs[order(total_pvs)]
    pVal_threshold <- total_pvs[floor(params[[6]] * length(total_pvs))]
    
    ### filter the interesting pathways with the p-value threshold
    for(j in 1:nrow(interesting_pathway_cnt)) {
      if(interesting_pathway_cnt[j,"GTEx"] == 0) {
        for(k in 1:length(tcga_regulon_pathways)) {
          pv <- tcga_regulon_pathways[[k]][rownames(interesting_pathway_cnt)[j],"p.adjust"]
          if(!is.na(pv) && (pv < pVal_threshold) && (length(which(gtex_aracne_hubs == names(tcga_regulon_pathways)[k])) > 0)) {
            interesting_pathway_cnt[j,"Hubs"] <- paste(interesting_pathway_cnt[j,"Hubs"], names(tcga_regulon_pathways)[k], sep = "/")
          }
        }
      } else {
        for(k in 1:length(gtex_regulon_pathways)) {
          pv <- gtex_regulon_pathways[[k]][rownames(interesting_pathway_cnt)[j],"p.adjust"]
          if(!is.na(pv) && (pv < pVal_threshold) && (length(which(tcga_aracne_hubs == names(gtex_regulon_pathways)[k])) > 0)) {
            interesting_pathway_cnt[j,"Hubs"] <- paste(interesting_pathway_cnt[j,"Hubs"], names(gtex_regulon_pathways)[k], sep = "/")
          }
        }
      }
    }
    interesting_pathway_cnt[,"Hubs"] <- substring(interesting_pathway_cnt[,"Hubs"], 2)
    
    ### order the interesting pathway count info
    interesting_pathway_cnt <- interesting_pathway_cnt[order(-interesting_pathway_cnt[,"Abs(Diff)"]),]
    interesting_pathway_cnt <- interesting_pathway_cnt[c(which(interesting_pathway_cnt[,"Hubs"] != ""), which(interesting_pathway_cnt[,"Hubs"] == "")),]
    
    ### write out the interesting pathway count info
    write.table(data.frame(GO_ID=rownames(interesting_pathway_cnt),
                           Pathway=goTermMap[rownames(interesting_pathway_cnt)],
                           interesting_pathway_cnt,
                           stringsAsFactors = FALSE, check.names = FALSE, row.names = NULL),
                file = paste0(params[[7]], subdirPath, "/", subdirPath, "_interesting_pathway_counts_", pVal_threshold, ".txt"),
                sep = "\t", row.names = FALSE)
    
    ### save the info to the list
    regulon_pathway_count_info[[i]] <- data.frame(GO_ID=rownames(interesting_pathway_cnt),
                                                  Pathway=goTermMap[rownames(interesting_pathway_cnt)],
                                                  interesting_pathway_cnt,
                                                  stringsAsFactors = FALSE, check.names = FALSE, row.names = NULL)
    
    ### quality control test
    ### read count depth of each sample
    common_genes <- intersect(rownames(gtex_gexp), rownames(tcga_gexp))
    gtex_gexp_colSums <- apply(gtex_gexp[common_genes,], 2, sum)
    tcga_gexp_colSums <- apply(tcga_gexp[common_genes,], 2, sum)
    
    ### make a data frame for a beeswarm plot
    bsdf <- data.frame(matrix(0, length(gtex_gexp_colSums) + length(tcga_gexp_colSums), 2))
    colnames(bsdf) <- c("Group", "Read_Depth")
    bsdf[1:length(gtex_gexp_colSums),"Group"] <- paste0("GTEx_", GTEx_TCGA_Map[i,"GTEx"])
    bsdf[(length(gtex_gexp_colSums)+1):nrow(bsdf),"Group"] <- paste0("TCGA_", toupper(GTEx_TCGA_Map[i,"TCGA"]))
    bsdf[1:length(gtex_gexp_colSums),"Read_Depth"] <- gtex_gexp_colSums
    bsdf[(length(gtex_gexp_colSums)+1):nrow(bsdf),"Read_Depth"] <- tcga_gexp_colSums
    
    ### draw a beeswarm plot with the read count depth of each sample
    ggplot(bsdf, aes(x=Group, y=Read_Depth)) +
      ggtitle(paste("Read Depth of GTEx", GTEx_TCGA_Map[i,"GTEx"], "vs TCGA", toupper(GTEx_TCGA_Map[i,"TCGA"]))) +
      theme_classic(base_size = 16) +
      geom_boxplot() +
      geom_beeswarm(aes(color=Group)) +
      stat_compare_means()
    ggsave(filename = paste0(params[[7]], subdirPath, "/", subdirPath, "_read_depth.png"), width = 12, height = 10)
    
    ### DE analysis
    gexp <- cbind(gtex_gexp[common_genes,], tcga_gexp[common_genes,])
    group <- c(rep("GTEx", ncol(gtex_gexp)), rep("TCGA", ncol(tcga_gexp)))
    deresult <- deseqWithComparisons(rCnt = gexp, grp = group, exp_class = "GTEx", ctrl_class = "TCGA")
    
    ### make a signature vector for GSEA
    de_sig <- deresult[,"stat"]
    names(de_sig) <- rownames(deresult)
    
    ### make a target gene list (like a pathway gene list in GSEA)
    target_genes <- vector("list", length = length(which(interesting_pathway_cnt[,"Hubs"] != "")))
    names(target_genes) <- rownames(interesting_pathway_cnt[which(interesting_pathway_cnt[,"Hubs"] != ""),])
    for(j in 1:length(target_genes)) {
      ### get all the found target genes of the given pathway across all the filtered regulons in the tissue
      tIdx <- strsplit(interesting_pathway_cnt[j,"Hubs"], "/", fixed = TRUE)[[1]]
      if(interesting_pathway_cnt[j,"GTEx"] == 0) {
        target_genes[[j]] <- Reduce(union, lapply(tcga_regulon_pathways[tIdx], function(x) {
          y <- x[rownames(interesting_pathway_cnt)[j],"geneID"]
          return(strsplit(y, "/", fixed = TRUE)[[1]])
        }))
      } else {
        target_genes[[j]] <- Reduce(union, lapply(gtex_regulon_pathways[tIdx], function(x) {
          y <- x[rownames(interesting_pathway_cnt)[j],"geneID"]
          return(strsplit(y, "/", fixed = TRUE)[[1]])
        }))
      }
      
      ### change the gene symbols to entrez ids
      target_genes[[j]] <- as.character(geneSymbolToEntrezId(target_genes[[j]]))
    }
    
    ### run gene set enrichment test
    set.seed(1234)
    fgseaRes <- fgsea(pathways = target_genes, stats = de_sig, nperm = 10000)
    
    ### some target genes may not have any de_sig at all
    ### only retain target genes that have at least one in intersection
    exist_idx <- sapply(target_genes, function(x) {
      if(length(intersect(x, names(de_sig))) > 0) {
        return(TRUE)
      } else {
        return(FALSE)
      }
    })
    
    ### add pathway labels
    fgseaRes <- data.frame(GO_ID=fgseaRes$pathway[exist_idx],
                           Pathway=goTermMap[fgseaRes$pathway[exist_idx]],
                           fgseaRes,
                           stringsAsFactors = FALSE, check.names = FALSE, row.names = NULL)
    
    ### write out the gsea result
    write.table(fgseaRes[,-which(colnames(fgseaRes) == "leadingEdge")],
                file = paste0(params[[7]], subdirPath, "/", subdirPath, "_GSEA_targets_", pVal_threshold, ".txt"),
                sep = "\t", row.names = FALSE)
    
    ### for interesting cases (pval < 0.05), print enrichment plot
    dir.create(paste0(params[[7]], subdirPath, "/DEG_GSEA_plots"), showWarnings = FALSE)
    for(j in 1:nrow(fgseaRes)) {
      if(fgseaRes$pval[j] < 0.05) {
        png(filename = paste0(params[[7]], subdirPath, "/DEG_GSEA_plots/GO_", substring(fgseaRes$GO_ID[j], 4), "_", pVal_threshold, ".png"),
            width = 1500, height = 1000, res = 200)
        print(plotEnrichment(target_genes[[fgseaRes$GO_ID[j]]], de_sig) + labs(title = fgseaRes$Pathway[j]))
        dev.off()
      }
    }
    
    ### get t-statistic of VIPER NES difference between GTEx and TCGA
    common_hubs <- intersect(rownames(gtex_viper), rownames(tcga_viper))
    t_diff <- sapply(1:length(common_hubs), function(x) {
      t <- t.test(gtex_viper[x,], tcga_viper[x,])
      return(t$statistic)
    })
    names(t_diff) <- common_hubs
    
    ### make a hub set list (like a pathway gene list in GSEA)
    hub_sets <- vector("list", length = length(which(interesting_pathway_cnt[,"Hubs"] != "")))
    names(hub_sets) <-  rownames(interesting_pathway_cnt[which(interesting_pathway_cnt[,"Hubs"] != ""),])
    for(j in 1:length(hub_sets)) {
      hub_sets[[j]] <- strsplit(interesting_pathway_cnt[j,"Hubs"], "/", fixed = TRUE)[[1]]
    }
    
    ### run gene set enrichment test
    set.seed(1234)
    fgseaRes2 <- fgsea(pathways = hub_sets, stats = t_diff, nperm = 10000)
    
    ### some hub sets may not have any t_diff at all
    ### only retain hub sets that have at least one in intersection
    exist_idx <- sapply(hub_sets, function(x) {
      if(length(intersect(x, names(t_diff))) > 0) {
        return(TRUE)
      } else {
        return(FALSE)
      }
    })
    
    ### add pathway labels
    fgseaRes2 <- data.frame(GO_ID=fgseaRes2$pathway[exist_idx],
                            Pathway=goTermMap[fgseaRes2$pathway[exist_idx]],
                            fgseaRes2,
                            stringsAsFactors = FALSE, check.names = FALSE, row.names = NULL)
    
    ### write out the gsea result
    write.table(fgseaRes2[,-which(colnames(fgseaRes2) == "leadingEdge")],
                file = paste0(params[[7]], subdirPath, "/", subdirPath, "_GSEA_hubs_", pVal_threshold, ".txt"),
                sep = "\t", row.names = FALSE)
    
    ### for interesting cases (pval < 0.05), print enrichment plot
    dir.create(paste0(params[[7]], subdirPath, "/DAH_GSEA_plots"), showWarnings = FALSE)
    for(j in 1:nrow(fgseaRes2)) {
      if(fgseaRes2$pval[j] < 0.05) {
        png(filename = paste0(params[[7]], subdirPath, "/DAH_GSEA_plots/GO_", substring(fgseaRes2$GO_ID[j], 4), "_", pVal_threshold, ".png"),
            width = 1500, height = 1000, res = 200)
        print(plotEnrichment(hub_sets[[fgseaRes2$GO_ID[j]]], t_diff) + labs(title = fgseaRes2$Pathway[j]))
        dev.off()
      }
    }
    
  }
  
  ### set README function
  README <- function() {
    writeLines(paste(rep("#", 100), collapse = ""))
    writeLines("A RDA file contains regulon pathway count info between")
    writeLines("GTEx and TCGA of 29 tissue mappings.")
    writeLines("The \"regulon_pathway_count_info\" object is a data frame that has")
    writeLines("exclusively appeared pathway info of the corresponding GTEx-TCGA")
    writeLines("comparison of the TISSUE.")
    writeLines(paste(rep("#", 100), collapse = ""))
  }
  
  ### save the DE result matrices in a RDA file
  save(list = c("regulon_pathway_count_info", "README"), file = paste0(params[[7]], "All_29_GTEx_vs_TCGA_regulon_pathway_info.rda"))
  
}

# ******************** investigate_interesting_regulon_pathways ********************
# For a given significant hub and a given interesting regulon pathway of a given tissue mapping,
# this function helps to identify what is happening with regard to the given inputs.
# The output will be consist of 6 results and they will be all combined into one PDF file.
# 1. A visualized network (graph) of the given hub and target genes.
#    The target genes found in the given pathway are marked with yellow.
#    From this graph, we could know which target genes the given hub have in both
#    GTEx and TCGA, and how many of them are associated with the given pathway.
# 2. A Venn diagram that describes how many target genes are shared between two regulons
#    of the given hub in GTEx and TCGA.
# 3. A table of Aracne info + DE analysis result of the found targets (Found pathway genes in the regulon)
#    Rows are the found targets, columns should be (Symbol, MI, MoA, Likelihood, Pvalue_MI,
#    baseMean, log2FoldChange, lfcSE, stat, pvalue_DE, padj_DE, Rank_DE)
# 4. A heatmap of the gene expression of the found target genes in both GTEx and TCGA
#    There should be a column side bar that represents where the samples are from (GTEx or TCGA).
# 5. A table and a line graph of MI comparison between GTEx and TCGA.
#    The mutual information is measured between the given hub and the found target genes.
#    Rows are the found targets and the number of column is 2 meaning GTEx and TCGA.
# 
# params[[1]]: The GO ID of interest
#              (a character vector of length 1)
# params[[2]]: The hub name of interest in Entrez ID
#              (a character vector of length 1)
# params[[3]]: The comparison of interest - the index in the GTEx-TCGA mapping (should be 1:29)
#              See the "GTEx_TCGA_Map" matrix for the details
#              e.g., params[[3]] = 21 indicates GTEx_Lung vs TCGA_LUAD
#              (an integer on [1, 29])
# params[[4]]: The file path of the "GTEx_TCGA_Map.rda" file
#              (a character vector of length 1)
# params[[5]]: The file path of Aracne network RDA file (All_62_ARACNE.rda)
#              (a character vector of length 1)
# params[[6]]: The file path of regulon pathway result RDA file (RegulonPathwayAnnotation.rda)
#              (a character vector of length 1)
# params[[7]]: The file path of DE results of the 29 GTEx-TCGA mapping (All_29_GTEx_vs_TCGA_DE_Results.rda)
#              (a character vector of length 1)
# params[[8]]: The file path of Aracne-ready gene expression RDA file (ALL_62_ARACNE_READY_EXPMAT.rda)
#              (a character vector of length 1)
# params[[9]]: The output results directory path
#              (a character vector of length 1)
#
# e.g., params = list("GO:0140013", "7272", 6,
#                     "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/GTEx_TCGA_Map.rda",
#                     "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/All_62_ARACNE.rda",
#                     "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/RegulonPathwayAnnotation.rda",
#                     "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/All_29_GTEx_vs_TCGA_DE_Results.rda",
#                     "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/ALL_62_ARACNE_READY_EXPMAT.rda",
#                     "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/regulon_pathway/GTEx_vs_TCGA/GTEx_BrainCerHem_vs_TCGA_GBM/)
# e.g., params = list("GO:0140013", "7272", 6,
#                     "./data/RDA_Files/GTEx_TCGA_Map.rda",
#                     "./data/RDA_Files/All_62_ARACNE.rda",
#                     "./data/RDA_Files/RegulonPathwayAnnotation.rda",
#                     "./data/RDA_Files/All_29_GTEx_vs_TCGA_DE_Results.rda",
#                     "./data/RDA_Files/ALL_62_ARACNE_READY_EXPMAT.rda",
#                     "./results/regulon_pathway/GTEx_vs_TCGA/GTEx_BrainCerHem_vs_TCGA_GBM/")

investigate_interesting_regulon_pathways <- function(params) {
  
  ### argument checking
  assertString(params[[1]])
  assertString(params[[2]])
  assertIntegerish(params[[3]])
  assertString(params[[4]])
  assertString(params[[5]])
  assertString(params[[6]])
  assertString(params[[7]])
  assertString(params[[8]])
  assertString(params[[9]])
  if(is.na(params[[3]]) || params[[3]] < 1 || params[[3]] > 29) {
    stop("params[[3]] should be an integer in 1:29")
  }
  
  ### load the data - only load the data of the given tissue
  load(params[[4]])
  
  env <- new.env()
  load(params[[5]], env)
  gtex_aracne_net <- env[[GTEx_TCGA_Map[params[[3]],"GTEx"]]]
  tcga_aracne_net <- env[[paste0("tcga_",GTEx_TCGA_Map[params[[3]],"TCGA"])]]
  gc()
  
  env <- new.env()
  load(params[[6]], env)
  gtex_regulon_pathway_result <- env[[paste0(GTEx_TCGA_Map[params[[3]],"GTEx"], "GO")]][[params[[2]]]]
  tcga_regulon_pathway_result <- env[[paste0("tcga_",GTEx_TCGA_Map[params[[3]],"TCGA"],"GO")]][[params[[2]]]]
  gc()
  
  env <- new.env()
  load(params[[7]], env)
  deresult <- env[[paste0("DEG_GTEx_", GTEx_TCGA_Map[params[[3]],"GTEx"], "_vs_TCGA_", toupper(GTEx_TCGA_Map[params[[3]],"TCGA"]))]]
  gc()
  
  env <- new.env()
  load(params[[8]], env)
  gtex_aracne_ready_gexp <- env[[paste0("expmat_gtex_",GTEx_TCGA_Map[params[[3]],"GTEx"])]]
  tcga_aracne_ready_gexp <- env[[paste0("expmat_tcga_",GTEx_TCGA_Map[params[[3]],"TCGA"])]]
  gc()
  
  ### get target genes
  gtex_regulon <- rownames(gtex_aracne_net[[2]][[params[[2]]]])
  tcga_regulon <- rownames(tcga_aracne_net[[2]][[params[[2]]]])
  gtex_regulon_symbol <- entrezIDtoSymbol(gtex_regulon)
  tcga_regulon_symbol <- entrezIDtoSymbol(tcga_regulon)
  
  ### determine GTEx exclusive or TCGA exclusive
  if(is.na(gtex_regulon_pathway_result[params[[1]],"geneID"])) {
    exclusive <- "TCGA"
  } else {
    exclusive <- "GTEx"
  }
  
  ### found pathway genes
  if(exclusive == "GTEx") {
    found_pathway_genes <- gtex_regulon[which(gtex_regulon_symbol %in% strsplit(gtex_regulon_pathway_result[params[[1]],"geneID"], split = "/", fixed = TRUE)[[1]])]
  } else {
    found_pathway_genes <- tcga_regulon[which(tcga_regulon_symbol %in% strsplit(tcga_regulon_pathway_result[params[[1]],"geneID"], split = "/", fixed = TRUE)[[1]])]
  }
  
  ### add rank column to deresult
  deresult <- data.frame(deresult, rank=rank(deresult[,"stat"], ties.method = "min"),
                         stringsAsFactors = FALSE, check.names = FALSE)
  
  ### create a sub directory for the results
  subdirPath <- paste0("GO_", substring(params[[1]], 4), "_Hub_", params[[2]])
  dir.create(file.path(params[[9]], subdirPath), showWarnings = FALSE)
  
  ### 1. A visualized network (graph) of the given hub and target genes.
  
  ### load required library
  if(!require(igraph, quietly = TRUE)) {
    install.packages("igraph")
    require(igraph, quietly = TRUE)
  }
  
  ### make an edge list for graph
  edgeList <- matrix(NA, length(gtex_regulon)+length(tcga_regulon), 2)
  edgeList[1:length(gtex_regulon),1] <- paste0("GTEx_", params[[2]])
  edgeList[1:length(gtex_regulon),2] <- gtex_regulon
  edgeList[(length(gtex_regulon)+1):nrow(edgeList),1] <- paste0("TCGA_", params[[2]])
  edgeList[(length(gtex_regulon)+1):nrow(edgeList),2] <- tcga_regulon
  
  ### visualize the graph
  net_graph <- graph.edgelist(edgeList, directed=FALSE)
  
  node_color <- rep("black", length(V(net_graph)))
  node_color[which(startsWith(V(net_graph)$name, "GTEx"))] <- "blue"
  node_color[which(startsWith(V(net_graph)$name, "TCGA"))] <- "blue"
  node_color[which(V(net_graph)$name == intersect(gtex_regulon, tcga_regulon))] <- "orange"
  node_color[which(V(net_graph)$name %in% as.character(found_pathway_genes))] <- "red"
  
  ### graph in Entrez ID
  png(paste0(params[[9]], subdirPath, "/Hub_", params[[2]], "_target_genes_entrez.png"),
      width = 1800, height = 1300, res = 150)
  plot(net_graph, main=paste0("Hub ", params[[2]], "'s Target Genes"),
       vertex.label.color=node_color, vertex.shape="none", layout=layout_with_dh)
  legend("topleft", xpd = TRUE, title = "Gene Color",
         legend = c("Hubs", "Shared Targets", paste0(params[[1]], "_Genes"), "The Other Targets"),
         fill = c("blue", "orange", "red", "black"), cex = 1, box.lty = 1)
  dev.off()
  
  ### graph in gene symbol
  edgeList <- matrix(NA, length(gtex_regulon)+length(tcga_regulon), 2)
  edgeList[1:length(gtex_regulon),1] <- paste0("GTEx_", entrezIDtoSymbol(params[[2]]))
  edgeList[1:length(gtex_regulon),2] <- entrezIDtoSymbol(gtex_regulon)
  edgeList[(length(gtex_regulon)+1):nrow(edgeList),1] <- paste0("TCGA_", entrezIDtoSymbol(params[[2]]))
  edgeList[(length(gtex_regulon)+1):nrow(edgeList),2] <- entrezIDtoSymbol(tcga_regulon)
  net_graph <- graph.edgelist(edgeList, directed=FALSE)
  
  png(paste0(params[[9]], subdirPath, "/Hub_", entrezIDtoSymbol(params[[2]]), "_target_genes_symbol.png"),
      width = 1800, height = 1300, res = 150)
  plot(net_graph, main=paste0("Hub ", entrezIDtoSymbol(params[[2]]), "'s Target Genes"),
       vertex.label.color=node_color, vertex.shape="none", layout=layout_with_dh)
  legend("topleft", xpd = TRUE, title = "Gene Color",
         legend = c("Hubs", "Shared Targets", paste0(params[[1]], "_Genes"), "The Other Targets"),
         fill = c("blue", "orange", "red", "black"), cex = 1, box.lty = 1)
  dev.off()
  
  ### 2. A Venn diagram between GTEx & TCGA regulons
  
  ### load required library
  if(!require(VennDiagram)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("VennDiagram")
    require(VennDiagram)
  }
  if(!require(gridExtra, quietly = TRUE)) {
    install.packages("gridExtra")
    require(gridExtra, quietly = TRUE)
  }
  
  ### draw the venn diagram
  v <- venn.diagram(list(gtex_regulon, tcga_regulon),
                    category.names = c("GTEx", "TCGA"),
                    cat.cex = 1.0, cex = 1.5,
                    filename = NULL)
  png(paste0(params[[9]], subdirPath, "/Hub_", params[[2]], "_regulon_venn.png"),
      width = 1200, height = 1000, res = 180)
  grid.arrange(gTree(children=v),
               top=paste0("Hub ", params[[2]], "'s Regulons"),
               bottom="")
  dev.off()
  
  ### 3. A table of Aracne info + DE analysis result of the found targets
  
  ### get Aracne info of the found targets
  if(exclusive == "GTEx") {
    aracne_de_table <- gtex_aracne_net[[2]][[params[[2]]]][found_pathway_genes,]
  } else {
    aracne_de_table <- tcga_aracne_net[[2]][[params[[2]]]][found_pathway_genes,]
  }
  
  ### get DE info of the found targts
  aracne_de_table <- merge(aracne_de_table, deresult[found_pathway_genes,], by="row.names", all=FALSE)
  rownames(aracne_de_table) <- aracne_de_table[,1]
  aracne_de_table <- aracne_de_table[,-c(1,2)]
  colnames(aracne_de_table) <- c(paste0("MI_",params[[2]]), paste0("MoA_",params[[2]]),
                                 paste0("Likelihood_",params[[2]]), "MI_Pval",
                                 "DE_baseMean", "DE_log2FC", "DE_lfcSE",
                                 "DE_stat", "DE_Pval", "DE_FDR", "DE_Rank")
  aracne_de_table <- data.frame(Target_Entrez_ID=rownames(aracne_de_table),
                                Gene_Symbol=entrezIDtoSymbol(rownames(aracne_de_table)),
                                aracne_de_table,
                                stringsAsFactors = FALSE, check.names = FALSE)
  
  ### write out the Aracne info + DE analysis result table
  write.table(aracne_de_table, file = paste0(params[[9]], subdirPath,
                                             "/Hub_", params[[2]], "_", exclusive,
                                             "_target_genes_info.txt"),
              sep = "\t", row.names = FALSE)
  
  ### 4. A heatmap of the gene expression of the found target genes in both GTEx and TCGA
  
  ### load required library
  if(!require(gplots, quietly = TRUE)) {
    install.packages("gplots")
    require(gplots, quietly = TRUE)
  }
  
  ### create a matrix for heatmap
  heatmap_mat <- merge(gtex_aracne_ready_gexp[which(rownames(gtex_aracne_ready_gexp) %in% found_pathway_genes),],
                       tcga_aracne_ready_gexp[which(rownames(tcga_aracne_ready_gexp) %in% found_pathway_genes),],
                       by="row.names", all=FALSE)
  rownames(heatmap_mat) <- heatmap_mat[,1]
  heatmap_mat <- heatmap_mat[,-1]
  
  ### set colside colors
  col_colors <- c(rep("#F8766D", ncol(gtex_aracne_ready_gexp)), rep("#619CFF", ncol(tcga_aracne_ready_gexp)))
  names(col_colors) <- c(rep("GTEx", ncol(gtex_aracne_ready_gexp)), rep("TCGA", ncol(tcga_aracne_ready_gexp)))
  
  ### create a heatmap
  png(paste0(params[[9]], subdirPath, "/Hub_", params[[2]], "_targets_expression_heatmap.png"),
      width = 1550, height = 1500, res = 120)
  par(oma=c(0,0,0,10))
  heatmap.3(as.matrix(heatmap_mat), main = paste0("Hub_", params[[2]], "'s Target Expressions"),
            xlab = "", ylab = "", col=greenred(100),
            scale="none", key=T, keysize=0.5, dendrogram = 'none', trace = 'none',
            labRow = rownames(heatmap_mat), labCol = "",
            Rowv = FALSE, Colv = FALSE,
            ColSideColors = cbind(as.vector(col_colors)),
            cexRow = 2.3, cexCol = 1.3, na.rm = TRUE)
  legend("topright", inset = -0.05, xpd = TRUE, title = "Sample Group",
         legend = unique(names(col_colors)), fill = unique(col_colors), cex = 1.3, box.lty = 0)
  dev.off()
  
  ### 5. A table and a line graph of MI comparison between GTEx and TCGA
  
  ### load required library
  if(!require(entropy, quietly = TRUE)) {
    install.packages("entropy")
    require(entropy, quietly = TRUE)
  }
  if(!require(ggplot2, quietly = TRUE)) {
    install.packages("ggplot2")
    require(ggplot2, quietly = TRUE)
  }
  
  ### create an empty mi table
  mi_table <- matrix(NA, length(found_pathway_genes), 8)
  rownames(mi_table) <- found_pathway_genes
  colnames(mi_table) <- c("GTEx_MI", "TCGA_MI", "GTEx_MI_PVal", "TCGA_MI_PVal", "GTEx_Pearson", "TCGA_Pearson", "GTEx_Spearman", "TCGA_Spearman")
  
  ### calculate MIs
  gtex_mi <- rep(NA, nrow(gtex_aracne_ready_gexp))
  names(gtex_mi) <- rownames(gtex_aracne_ready_gexp)
  tcga_mi <- rep(NA, nrow(tcga_aracne_ready_gexp))
  names(tcga_mi) <- rownames(tcga_aracne_ready_gexp)
  set.seed(1234)
  for(gene in rownames(gtex_aracne_ready_gexp)) {
    gtex_mi[gene] <- mi.plugin(rbind(gtex_aracne_ready_gexp[params[[2]],],
                                     gtex_aracne_ready_gexp[gene,]))
  }
  for(gene in rownames(tcga_aracne_ready_gexp)) {
    tcga_mi[gene] <- mi.plugin(rbind(tcga_aracne_ready_gexp[params[[2]],],
                                     tcga_aracne_ready_gexp[gene,]))
  }
  gtex_mi <- gtex_mi[order(gtex_mi)]
  tcga_mi <- tcga_mi[order(tcga_mi)]
  
  ### calculate Correlations
  for(target in found_pathway_genes) {
    if(length(which(rownames(gtex_aracne_ready_gexp) == target)) > 0) {
      ### mi calculation
      mi_table[target,"GTEx_MI"] <- gtex_mi[target]
      
      ### mi p-value
      mi_table[target,"GTEx_MI_PVal"] <- which(names(gtex_mi) == target) / length(gtex_mi)
      
      ### pearson
      mi_table[target,"GTEx_Pearson"] <- cor(gtex_aracne_ready_gexp[params[[2]],],
                                             gtex_aracne_ready_gexp[target,],
                                             method = "pearson")
      ### spearman
      mi_table[target,"GTEx_Spearman"] <- cor(gtex_aracne_ready_gexp[params[[2]],],
                                              gtex_aracne_ready_gexp[target,],
                                              method = "spearman")
    }
    if(length(which(rownames(tcga_aracne_ready_gexp) == target)) > 0) {
      ### mi calculation
      mi_table[target,"TCGA_MI"] <- tcga_mi[target]
      
      ### mi p-value
      mi_table[target,"TCGA_MI_PVal"] <- which(names(tcga_mi) == target) / length(tcga_mi)
      
      ### pearson
      mi_table[target,"TCGA_Pearson"] <- cor(tcga_aracne_ready_gexp[params[[2]],],
                                             tcga_aracne_ready_gexp[target,],
                                             method = "pearson")
      ### spearman
      mi_table[target,"TCGA_Spearman"] <- cor(tcga_aracne_ready_gexp[params[[2]],],
                                              tcga_aracne_ready_gexp[target,],
                                              method = "spearman")
    }
  }
  
  ### write out the mi table
  write.table(data.frame(Target_Entrez_ID=rownames(mi_table),
                         Gene_Symbol=entrezIDtoSymbol(rownames(mi_table)),
                         mi_table,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(params[[9]], subdirPath,
                            "/Hub_", params[[2]], "_targets_mi_comparison.txt"),
              sep = "\t", row.names = FALSE)
  
  ### make a data frame for a line graph
  line_df <- data.frame(matrix(NA, 2*nrow(mi_table), 6))
  colnames(line_df) <- c("Gene", "MI", "MI_PVal", "Pearson", "Spearman", "Group")
  line_df[1:nrow(mi_table),"Gene"] <- rownames(mi_table)
  line_df[1:nrow(mi_table),"MI"] <- mi_table[,"GTEx_MI"]
  line_df[1:nrow(mi_table),"MI_PVal"] <- mi_table[,"GTEx_MI_PVal"]
  line_df[1:nrow(mi_table),"Pearson"] <- mi_table[,"GTEx_Pearson"]
  line_df[1:nrow(mi_table),"Spearman"] <- mi_table[,"GTEx_Spearman"]
  line_df[1:nrow(mi_table),"Group"] <- "GTEx"
  line_df[(nrow(mi_table)+1):nrow(line_df),"Gene"] <- rownames(mi_table)
  line_df[(nrow(mi_table)+1):nrow(line_df),"MI"] <- mi_table[,"TCGA_MI"]
  line_df[(nrow(mi_table)+1):nrow(line_df),"MI_PVal"] <- mi_table[,"TCGA_MI_PVal"]
  line_df[(nrow(mi_table)+1):nrow(line_df),"Pearson"] <- mi_table[,"TCGA_Pearson"]
  line_df[(nrow(mi_table)+1):nrow(line_df),"Spearman"] <- mi_table[,"TCGA_Spearman"]
  line_df[(nrow(mi_table)+1):nrow(line_df),"Group"] <- "TCGA"
  
  ### line graphs
  p <- vector("list", length = 4)
  p[[1]] <- ggplot(data = line_df, aes(x=Gene, y=MI, group=Group)) +
    geom_line(aes(color=Group)) +
    theme_classic(base_size = 16)
  p[[2]] <- ggplot(data = line_df, aes(x=Gene, y=MI_PVal, group=Group)) +
    geom_line(aes(color=Group)) +
    theme_classic(base_size = 16)
  p[[3]] <- ggplot(data = line_df, aes(x=Gene, y=Pearson, group=Group)) +
    geom_line(aes(color=Group)) +
    theme_classic(base_size = 16)
  p[[4]] <- ggplot(data = line_df, aes(x=Gene, y=Spearman, group=Group)) +
    geom_line(aes(color=Group)) +
    theme_classic(base_size = 16)
  
  ### draw out the graphs
  png(paste0(params[[9]], subdirPath, "/Hub_", params[[2]], "_targets_mi_comparison.png"),
      width = 2000, height = 1000)
  multiplot(plotlist = p, cols = 1, title = paste0("Hub_", params[[2]], "_Targets_Correlation_Comparisons"))
  dev.off()
  
}
