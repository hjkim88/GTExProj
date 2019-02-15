### this is a sanity check code after running /ifs/scratch/hk2990/Pathway/scripts/arrayRunMaster.sh

### get aracne networks files
f <- list.files(path = "X:/hk2990/Pathway/data/") 
f <- f[which(endsWith(f, ".rda"))]

### extract tissue names
tissue_names <- sapply(f, function(x) substr(x, 1, nchar(x)-4), USE.NAMES = FALSE)

### retreive number of result files for each tissue
outfile_num <- sapply(tissue_names, function(x) {
  f2 <- list.files(path = paste0("X:/hk2990/Pathway/runs/", x))
  f2 <- f2[which(endsWith(f2, ".rda"))]
  return (length(f2)-1)
})

### retreive number of jobs originally submitted to the cluster
job_num <- sapply(tissue_names, function(x) {
  f3 <- list.files(path = paste0("X:/hk2990/Pathway/runs/", x, "/out/"), full.names = TRUE)
  f3 <- f3[sapply(f3, file.size) > 0]
  return(length(f3))
})

### check the numbers of input jobs and the numbers of output results are the same
isOK = all.equal(as.integer(outfile_num), as.integer(job_num))

error_tissue <- NULL
for(i in 1:length(outfile_num)) {
  if(outfile_num[i] != job_num[i]) {
    error_tissue <- c(error_tissue,  tissue_names[i])
    writeLines(paste(tissue_names[i], " outfile_num[", i, "] = ", outfile_num[i], ", job_num[", i, "] = ", job_num[i]))
  }
}

### a function to load RDA and returns as a variable
loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# ******************************************************************************************
# Pathway Analysis with clusterProfiler package
# Input: geneList     = a vector of gene Entrez IDs for pathway analysis [numeric or character]
#        org          = organism that will be used in the analysis ["human" or "mouse"]
#                       should be either "human" or "mouse"
#        database     = pathway analysis database (KEGG or GO) ["KEGG" or "GO"]
#        title        = title of the pathway figure [character]
#        pv_threshold = pathway analysis p-value threshold (not DE analysis threshold) [numeric]
#        displayNum   = the number of pathways that will be displayed [numeric]
#                       (If there are many significant pathways show the few top pathways)
#        imgPrint     = print a plot of pathway analysis [TRUE/FALSE]
#        dir          = file directory path of the output pathway figure [character]
#
# Output: Pathway analysis results in figure - using KEGG and GO pathways
#         The x-axis represents the number of DE genes in the pathway
#         The y-axis represents pathway names
#         The color of a bar indicates adjusted p-value from the pathway analysis
#         For Pathview Result, all colored genes are found DE genes in the pathway,
#         and the color indicates log2(fold change) of the DE gene from DE analysis
# ******************************************************************************************
pathwayAnalysis_CP <- function(geneList,
                               org = "human",
                               database,
                               title="Pathway_Results",
                               pv_threshold=0.05,
                               displayNum=Inf,
                               imgPrint=FALSE,
                               dir="./") {
  
  ### load library
  if(!require(clusterProfiler)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("clusterProfiler")
    library(clusterProfiler)
  }
  if(!require(ggplot2)) {
    install.packages("ggplot2")
    library(ggplot2)
  }
  
  
  ### colect gene list (Entrez IDs)
  geneList <- geneList[which(!is.na(geneList))]
  
  if(!is.null(geneList)) {
    ### make an empty list
    p <- list()
    
    if(database == "KEGG") {
      ### KEGG Pathway
      kegg_enrich <- enrichKEGG(gene = geneList, organism = org, pvalueCutoff = pv_threshold)
      
      if(is.null(kegg_enrich)) {
        writeLines("KEGG Result does not exist")
        return(NULL)
      } else {
        kegg_enrich@result <- kegg_enrich@result[which(kegg_enrich@result$p.adjust < pv_threshold),]
        
        if(imgPrint == TRUE) {
          if((displayNum == Inf) || (nrow(kegg_enrich@result) <= displayNum)) {
            result <- kegg_enrich@result
            description <- kegg_enrich@result$Description
          } else {
            result <- kegg_enrich@result[1:displayNum,]
            description <- kegg_enrich@result$Description[1:displayNum]
          }
          
          if(nrow(kegg_enrich) > 0) {
            p[[1]] <- ggplot(result, aes(x=Description, y=Count)) + labs(x="", y="Gene Counts") + 
              theme_classic(base_size = 16) + geom_bar(aes(fill = p.adjust), stat="identity") + coord_flip() +
              scale_x_discrete(limits = rev(description)) +
              guides(fill = guide_colorbar(ticks=FALSE, title="P.Val", barheight=10)) +
              ggtitle(paste0("KEGG ", title))
            
            png(paste0(dir, "kegg_", title, "_CB.png"), width = 2000, height = 1000)
            print(p[[1]])
            dev.off()
          } else {
            writeLines("KEGG Result does not exist")
          }
        }
        
        return(kegg_enrich@result)
      }
    } else if(database == "GO") {
      ### GO Pathway
      if(org == "human") {
        go_enrich <- enrichGO(gene = geneList, OrgDb = 'org.Hs.eg.db', readable = T, ont = "BP", pvalueCutoff = pv_threshold)
      } else if(org == "mouse") {
        go_enrich <- enrichGO(gene = geneList, OrgDb = 'org.Mm.eg.db', readable = T, ont = "BP", pvalueCutoff = pv_threshold)
      } else {
        go_enrich <- NULL
        writeLines(paste("Unknown org variable:", org))
      }
      
      if(is.null(go_enrich)) {
        writeLines("GO Result does not exist")
        return(NULL)
      } else {
        go_enrich@result <- go_enrich@result[which(go_enrich@result$p.adjust < pv_threshold),]
        
        if(imgPrint == TRUE) {
          if((displayNum == Inf) || (nrow(go_enrich@result) <= displayNum)) {
            result <- go_enrich@result
            description <- go_enrich@result$Description
          } else {
            result <- go_enrich@result[1:displayNum,]
            description <- go_enrich@result$Description[1:displayNum]
          }
          
          if(nrow(go_enrich) > 0) {
            p[[2]] <- ggplot(result, aes(x=Description, y=Count)) + labs(x="", y="Gene Counts") + 
              theme_classic(base_size = 16) + geom_bar(aes(fill = p.adjust), stat="identity") + coord_flip() +
              scale_x_discrete(limits = rev(description)) +
              guides(fill = guide_colorbar(ticks=FALSE, title="P.Val", barheight=10)) +
              ggtitle(paste0("GO ", title))
            
            png(paste0(dir, "go_", title, "_CB.png"), width = 2000, height = 1000)
            print(p[[2]])
            dev.off()
          } else {
            writeLines("GO Result does not exist")
          }
        }
        
        return(go_enrich@result)
      }
    } else {
      stop("database prameter should be \"GO\" or \"KEGG\"")
    }
  } else {
    writeLines("geneList = NULL")
  }
}

# *****************************************************************************
# Run the GO enrichment analysis. Keeps track of the execution time, for
# logging purposes
# 
# ARGUMENTS:
# * net:		the slice of the interactome on which to run the analysis.
# * goOnt:		string specifying which GO ontology to use (BP, MF, CC).
# * logFile:	file for storing logging info.
# *****************************************************************************
doEnrichment <- function(net, org = "human", Ont = "GO", imgPrint = FALSE, logFile = NULL){
  start = Sys.time()
  res = lapply(net, function(x, org, Ont, imgPrint){
    return(pathwayAnalysis_CP(abs(x[,1]), org = org, database = Ont, imgPrint = imgPrint))
  }, org, Ont, imgPrint)
  end = Sys.time()
  names(res) = names(net)
  if (!is.null(logFile))
    cat(paste("\n\n\tElapsed time -> ", paste(toString(end-start), attr( end-start, "units"))), file = logFile, append=TRUE, sep="\n")
  return(res)
}

### run the remaining jobs
if(isOK != TRUE) {
  chunk_size = 30
  ontology = "GO"
  
  for(i in 1:length(error_tissue)) {
    
    ### get result files
    f4 <- list.files(path = paste0("X:/hk2990/Pathway/runs/", error_tissue[i]))
    f4 <- f4[which(startsWith(f4, paste0(error_tissue[i], "Res")))]
    
    ### get chunk ids from the result files
    f4 <- sapply(f4, function(x) substr(x, 1, nchar(x)-4), USE.NAMES = FALSE)
    f4 <- substring(f4, nchar(error_tissue[i])+5)
    f4 <- as.integer(f4)
    f4 <- f4[order(f4)]
    
    ### load the aracne network for the tissue
    net <- loadRData(fileName = paste0("X:/hk2990/Pathway/data/", error_tissue[i], ".rda"))
    
    ### the number of total hubs in the tissue
    L <- length(net[[2]])	
    
    ### the ids that should be in the result
    answer <- seq(1, L, chunk_size)
    
    ### the chunks that do not have results
    chunk_needed <- setdiff(answer, f4)
    
    ### do the remaining jobs
    for(j in 1:length(chunk_needed)) {
      res <- doEnrichment(net[[2]][chunk_needed[j]:(chunk_needed[j]+chunk_size)], Ont = ontology)
      resVarName <- paste0(error_tissue[i], "Res_", chunk_needed[j])
      fileName <- paste0(error_tissue[i], "Res_", chunk_needed[j], ".rda")
      assign(resVarName, res, envir = globalenv())
      save(list = c(resVarName), file=paste0("X:/hk2990/Pathway/runs/", error_tissue[i], "/", fileName))
    }
    
  }
  
}




