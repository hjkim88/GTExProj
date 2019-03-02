# *****************************************************************************
# This file contains common functions used by different pieces of code
# *****************************************************************************

### load library
if(!require("xlsx", quietly = TRUE)) {
  install.packages("xlsx")
  require(xlsx, quietly = TRUE)
}
if(!require("checkmate", quietly = TRUE)) {
  install.packages("checkmate")
  require("checkmate", quietly = TRUE)
}
if(!require("dendextend", quietly = TRUE)) {
  install.packages("dendextend")
  require("dendextend", quietly = TRUE)
}
if(!require("viper", quietly = TRUE)) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("viper")
  require("viper", quietly = TRUE)
}
if(!require("annotate", quietly = TRUE)) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("annotate")
  require("annotate", quietly = TRUE)
}
if(!require("org.Hs.eg.db", quietly = TRUE)) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("org.Hs.eg.db")
  require("org.Hs.eg.db", quietly = TRUE)
}
if(!require("GO.db", quietly = TRUE)) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("GO.db")
  require("GO.db", quietly = TRUE)
}
if(!require("topGO", quietly = TRUE)) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("topGO")
  require("topGO", quietly = TRUE)
}

# Load the list that maps gene ids to GO terms
# FIXME: this is only for human genes, should name the variable more clearly
if (!exists("geneToGoMap"))	
  load("Data/geneToGoMap.rda")

# Load the vector that maps gene names to their type (TF, co_TF, Signaling)
if (!exists("geneTypeMap"))	
  load("Data/geneTypeMap.rda")

# Convenience array listing of GO evidence types
goEvidenceTypes = c("inferred from mutant phenotype", "inferred from genetic interaction", 
		"inferred from physical interaction", "inferred from sequence similarity", "inferred from direct assay", 
		"inferred from expression pattern", "inferred from electronic annotation", "traceable author statement", 
		"non-traceable author statement", "no biological data available", "inferred by curator")
names(goEvidenceTypes) = c("IMP", "IGI", "IPI", "ISS", "IDA", "IEP", "IEA", "TAS", "NAS", "ND", "IC")

# Convenience variable: a named vector mapping GO terms (the vector names) to 
# their corresponding term description (the vector values):
if (!exists("goTermMap"))
	goTermMap = Term(GOTERM)

# *****************************************************************************
#
# Write the results of a MARINa analysis to the output "file_name" by appending 
# a worksheet with the MARINa results. If synergy analysis was run as well then
# a second worksheet is also attached listing the detailed results of the 
# synergy calculation.
# 
# * mrs:		An msviper object such as that returned by the methods msviper(),
#		msviperCombinatorial(), msviperSynergy()
# * file_name:	The file name where to write the results.
# * worksheet:	The name for the 2 worksheets (the second will just be named
#		"<worksheet> - Synergy"
# * format:		Specifies if the results will be stored in an Excel or a R
#		binary data file.
# *****************************************************************************
writeMarinaResults = function(mrs, file_name, worksheet, regul, format = c("xlsx", "rda")){
   # get the symbol names, the nes, the p-values, the regulon sizes, and the ledge genes.
   names = names(mrs$es$nes)
   nes   = mrs$es$nes
   pvals = mrs$es$p.value
   sizes = mrs$es$size
   ledge = mrs$ledge
   
   # Variables used to store results
	marina_hubs = marina_synergy = NULL
   # generate FDR values
   fdr = signif(p.adjust(mrs$es$p.value, "fdr"), 3)
   
   # convert 'ledge" into a character vector, to prepare for output
   ledge_genes = vector(mode="character", length=length(names))
   ledge_sizes = vector(mode="integer", length=length(names))
   for (i in 1:length(ledge)){
      ledge_sizes[i] = length(ledge[[i]])
      ledge_genes[i] = toString(ledge[[i]])
   }
   
   if (length(ledge) < length(names))
      for (i in (length(ledge)+1):length(names)){
         ledge_sizes[i] = 0
         ledge_genes[i] = ""
      }
   
   # order everything by increasing FDR and write to file
	pos = order(fdr)
	marina_hubs = data.frame(Gene = names[pos], NES = nes[pos], FDR = fdr[pos], P_value = pvals[pos], 
			Regulon_size = sizes[pos], Ledge_size = ledge_sizes[pos], Ledge = ledge_genes[pos])
	if (format == "xlsx")
		write.xlsx(marina_hubs, file = file_name, row.names = FALSE, col.names = TRUE, sheetName = worksheet, append=TRUE)

   # prepare to write out the second worksheet, containing the synergy results   
   L = length(ledge)
   N = length(names)
   if (L < N){
      names = names[(L+1):N]
      nes   = nes[(L+1):N]
      pvals = pvals[(L+1):N]
      sizes = sizes[(L+1):N]
      fdr   = fdr[(L+1):N]

      ledge_genes = vector(mode="character", N-L+1)
      ledge_sizes = vector(mode="integer", N-L+1)
      for (i in 1:length(names)){
         gene_pair_names =  strsplit(c(names[i]), c("--"))[[1]]
         reg1 = names(regul[[gene_pair_names[1]]]$tfmode)
         reg2 = names(regul[[gene_pair_names[2]]]$tfmode)
         ledge_sizes[i] = length(intersect(reg1, reg2))
         ledge_genes[i] = toString(intersect(reg1, reg2))
      }
      
	  pos = order(fdr)
	  marina_synergy = data.frame(Gene = names[pos], NES = nes[pos], FDR = fdr[pos], P_value = pvals[pos], 
			  Common_enriched = sizes[pos], Shared_regulon = ledge_sizes[pos], Shared_genes = ledge_genes[pos])
	  if (format == "xlsx")
		  write.xlsx(marina_synergy, file = file_name, row.names = FALSE, col.names = TRUE, 
				  sheetName = paste(worksheet, " - Synergy"), append=TRUE)
  
   }
   
   if (format == "rda")
	   save(list = c("marina_hubs", "marina_synergy"), file = file_name)
}


# *****************************************************************************
#
# Plot a dendrogram with colored labels.
#
# hc:		The results of an hclust() call.
# groupings:	A vector where names(groupings) contains the labels for the hc
#		leaf nodes and groupings[i] is an integer indicatings which 
#		group the leaf names(groupings)[i] belongs to. The labels od 
#		all leaves in the same group will be printed with the same font
#		color. Often we will want to call this methods as:
#			plotColoredCluster(hc, cutree(hc, <N>))
#		where <N> is an interger that instructs cutree to aumatically
#		create a "groupings" vector by splicing the cluster into <N>
#		groups, according to the structure of the tree
# hangVal:	A graphics parameter, to instruct how the labels should be hang
#		from the leaf branches (most of the time the defualt value 
#		should be fine).
# title:	Text to use as the title of the plot.
# *****************************************************************************
plotColoredCluster <- function(hc, groupings, hangVal = 0.1, title = "") {
  
  # The different colors to use for the various label groups. If there are 
  # more groups than distinct colors, cycle through
  # Red, blue, green, medium-green, black, magenta, orange, dark green, brown
  labelColors = c("#FF0000", "#0000FF", "#00FF00", "#14A890", "#000000", "#FF00FF", "#FFCC99", "#196307", "#996600")
  L = length(labelColors)
  
  labels = hc$labels
  if (length(setdiff(names(groupings), labels)) > 0)
     stop("Some labels in groupings are not in the hc ojbect")
     
  # function to get color labels
  colLab <- function(n) {
    if (is.leaf(n)) {
      a <- attributes(n)
      ind = groupings[which(names(groupings) == a$label)] %% L
      if (ind == 0)
        ind = L
      attr(n, "nodePar") <- c(a$nodePar, lab.col = labelColors[ind])
    }
    n
  }
  
  hcd = as.dendrogram(hc)
  clusDendro = dendrapply(hcd, colLab)
  #plot(clusDendro)
  plot(hang.dendrogram(clusDendro, hang=hangVal), main = title)
}


# *****************************************************************************
#
# Map Entrez Gene IDs to Gene Symbols
#
# geneID:	A single integer, or a vector of integers, or a single strings, or
#			a vector of strings. In each case, integers or strings are 
#			Entrez gene ids.
# 
# Returns a vector of the same size as "geneID" where the i-th entry is the 
# gene symbol corresponsing to the i-th Entrez id in "geneID". The return 
# vector entries also are named using the Entrez Ids. 
#
# -----------------------------------------------------------------------------
# -------FIXME: Modify the code to use "select", as in method geneSymbolToEntrezId
# -----------------------------------------------------------------------------
# *****************************************************************************
entrezIDtoSymbol <- function(geneID){
  if (is.numeric(geneID))
    geneID = as.character(geneID)
  return(getSYMBOL(geneID, data="org.Hs.eg"))
}

# *****************************************************************************
#
# Map Entrez Gene IDs to Gene Descriptions
#
# geneID:	A single integer, or a vector of integers, or a single strings, or
#			a vector of strings. In each case, integers or strings are 
#			Entrez gene ids.
# 
# Returns a vector of the same size as "geneID" where the i-th entry is the 
# descriptipm corresponsing to the i-th Entrez id in "geneID". The return 
# vector entries also are named using the Entrez Ids.
# -----------------------------------------------------------------------------
# -------FIXME: Modify the code to use "select", as in method geneSymbolToEntrezId
# -----------------------------------------------------------------------------
# *****************************************************************************
entrezIDtoDescription <- function(geneID){
	map <- org.Hs.egGENENAME
	if (is.numeric(geneID))
		geneID = as.character(geneID)
	res = sapply(geneID, function(e){return(map[[e]])})
	res = sapply(res, function(x){
				if (is.null(x))
					return("")
				return(x)
			})
	return(res)
}

# *****************************************************************************
#
# Map Gene Symbols to Gene IDs
#
# geneSymbol:	A single string or a vector of strings representing gene symbol(s).
# 
# Returns a vector of the same size as "geneSymbol" where the i-th entry is the 
# gene ID (as an integer) corresponsing to the i-th gene symbol in "geneSymbol". 
# The return vector entries are named using the gene symbols. For gene symbols mapped
# to more than one entrez ids, only the first id is returned. If the i-th entry is
# a gene alias, then entrez Id of its corresponding official gene symbol is returned.
#
# ATTENTION: 
# Before used geneSymbol is first stripped of symbols which are not mapped to 
# at least one entrez ID. If no gene symbols remain after this stripping, a NULL 
# value is returned.
# *****************************************************************************
geneSymbolToEntrezId <- function(geneSymbol){
  # res = as.integer(unlist(mget(geneSymbol, org.Hs.egSYMBOL2EG)))
  
  keys = keys(org.Hs.eg.db,  keytype="ALIAS")
  geneSymbol = geneSymbol[geneSymbol %in% keys]
  if (length(geneSymbol) == 0)
    return(NA)
  
  tmp = select(org.Hs.eg.db, keys=geneSymbol, keytype="ALIAS", columns=c("ENTREZID"))
  # Remove duplicate mappings, if any
  t = which(duplicated(tmp[,1]))
  if (length(t) > 0)
    tmp = tmp[-t, ]
  res = as.integer(tmp[,2])
  names(res) = tmp[, 1]
  return(res)
}


getChromosome <- function(genes, numeric = TRUE){
	g = as.entrezId(genes)
	map <- org.Hs.egCHR
	res = sapply(g, function(e){
				c = map[[e]] 
				# Handle sex chromosomes, taking into consideration the possiblity
				# that there could be no valid mapping for a gene ID. This can happen, 
				# e.g., if a gene ID has been depricated.
				if (is.null(c) || is.na(c))
					c = NA
				else if((c=="X" || c=="Y") && numeric)
					c = "23"
				else if ((c=="X" || c=="Y"))
					c = toString(c)
				return(c)})
	if (numeric)
		res = strtoi(res)
	names(res) = genes
	return(res)
}

# *****************************************************************************
# Check if one or more genes are TFs.
#
# gids:	A single integer, or a vector of integers, or a single string, or
#		a vector of strings. If 'gids' contains strings then these can be either 
#		*all* gene symbols or *all* Entrez gene ids. If 'gids' contains integers
#		then these must be Entrez ids.
# 
# Returns a vector V of the same size as "gids" where the i-th entry is TRUE or 
# FASLE depending on if the gene gids[i] is a transction factor or not. V is 
# named after the gene IDs of the query genes(s), specifically, names(V) = gids.
# *****************************************************************************
is.tf <- function(gids){
	gids = strtoi(as.entrezId(gids))
	return(check.gene.type(gids, 1))
}

# *****************************************************************************
# Check if one or more genes are coTFs.
#
# gids:	A single integer, or a vector of integers, or a single string, or
#		a vector of strings. If 'gids' contains strings then these can be either 
#		*all* gene symbols or *all* Entrez gene ids. If 'gids' contains integers
#		then these must be Entrez ids.
# 
# Returns a vector V of the same size as "gids" where the i-th entry is TRUE or 
# FALSE depending on if the gene gids[i] is a co-transction factor or not. V is 
# named after the gene IDs of the query genes(s), specifically, names(V) = gids.
# *****************************************************************************
is.cotf <- function(gids){
	gids = strtoi(as.entrezId(gids))
	return(check.gene.type(gids, 2))
}

# *****************************************************************************
# Check if one or more genes code signaling proteins.
#
# gids:	A single integer, or a vector of integers, or a single string, or
#		a vector of strings. If 'gids' contains strings then these can be either 
#		*all* gene symbols or *all* Entrez gene ids. If 'gids' contains integers
#		then these must be Entrez ids.
# 
# Returns a vector V of the same size as "gids" where the i-th entry is TRUE or 
# FALSE depending on if the gene gids[i] is a signaling protein or not. V is 
# named after the gene IDs of the query genes(s), specifically, names(V) = gids.
# *****************************************************************************
is.sign <- function(gids){
	gids = strtoi(as.entrezId(gids))
	return(check.gene.type(gids, 3))
}

# *****************************************************************************
# Annotate input genes as "TF", "co-TF", "SIGN", or OTHER (coded as NA).
#
# gids:	A single integer, or a vector of integers, or a single string, or
#		a vector of strings. If 'gids' contains strings then these can be either 
#		*all* gene symbols or *all* Entrez gene ids. If 'gids' contains integers
#		then these must be Entrez ids.
# 
# Returns a vector V of the same size as "gids" where the i-th entry is either
# "TF", "coTF", "SIGN" or NA depending on if the gene gids[i] is a transription
# factor, a co-transcription factor, a signaling protein, or none of these. The
# vectir is named after the query genes(s), specifically, names(V) = gids.
# *****************************************************************************
gene.type <-function(gids){
	gids = strtoi(as.entrezId(gids))
	res = rep(NA, length(gids))
	res[check.gene.type(gids, 1)] = "TF"
	res[check.gene.type(gids, 2)] = "co-TF"
	res[check.gene.type(gids, 3)] = "SIGN"
	names(res) = gids
	return(res)
}

# *****************************************************************************
# Helper function, does the dirty work for the funcitons is.tf, is.cotf, is.sign.
#
# gids:	A character vector containing Entrez IDs 
# val:	The integer value that codes the gene types in the matrix geneTypeMap:
#		* for TFs, val = 1.
#		* for co-TFs, val = 2
#		* for singaling proteins, val = 3
#
# Returns a boolean vector V with of size length(gids), where V[i] == TRUE iff
# gids[i] od of type == 'val'.
# *****************************************************************************

check.gene.type <- function (gids, val){
	res = rep(FALSE, length(gids))
	names(res) = gids
	for (i in 1:length(gids)){
		ind = findInterval(gids[i], geneTypeMap[,1])
		if((ind > 0) && (geneTypeMap[ind,1] == gids[i]))
			res[i] = (geneTypeMap[ind,2] == val)
	}
	return(res)
}


# *****************************************************************************
# Returns the most consistently highly expressed genes in an expression matrix.
# Specifically, genes that are found among the 'top' most highly expressed across
# mulitple samples.
#
# * expmat:	Gene expression matrix with columns corresponsing to samples and 
#			rows corresponding to genes. Rows are expected to be named with
#			gene Entrez ids.
# * top:	Positive integer
# 
# Retunrs a table-vector V where entries are named with entrez ids and for an
# id G the entry V["G"] is the number of samples where the gene G is among
# the 'top' most highly expressed genes. Entries in V are ordered according to
# that number (decreasing order). 
# *****************************************************************************
getTopExpressed <- function(expmat, top = 50){
	
	# First "reverse" the values in the expession matrix, so that the "rank"
	# method will work the way we want
	max_val = max(expmat)
	expmat = max_val - expmat
	
	# Rank the column values
	for (i in 1:ncol(expmat)) 
		expmat[,i] = rank(expmat[,i])
	
	return(sort(table(apply(expmat, 2, function(x){
										return(strtoi(names(sort(x))[1:top]))
									})), decreasing = TRUE))
}




# *****************************************************************************
# Return the GO terms associated with the query gene(s)
#
# ARGUMENTS:
# * genes:	the query genes, as a vector of Entrez ids (character or integer) or
# 			gene symbols.
# * ontology: a vector containing one or more of the strings "BP", "MF", "CC".
#			Designates the ontologies for which to report terms. If NULL then
#			terms of all ontologies are returned.
# * evidenceType: a character vector containing one or evidence types, i.e.,
#			one or more strings from 'names(goEvidenceTypes)'. Only terms with
#			evidence types among thost listed in this vector will be returned.
#			If NULL then terms for all evidence types are returned.
# * def:	If TRUE, the definitions of each GO term in also included in the 
#			results object.
#
# RETURNS
# A named list with length(genes) elements, one for each query gene. The i-th 
# list entry is named genes[i] and is a data frame with one row per GO term 
# that passes the filters specified by the funciton arguments 'ontology' and 
# 'evidenceType'. Each data frame row lists the GO term id, the GO term 
# mame, the evidence type; and, if 'def' = TRUE, the full (long...) GO term
# definition.
#
# EXAMPLE USE:
# The examples below use gene symbols as queries. These could be as well replaced 
# with Entrez IDs.
#           getGoTerms(c("FOXM1", "TPT1"))
#			getGoTerms(c("FOXM1", "TPT1"), "BP")
#			getGoTerms(c("FOXM1", "TPT1"), c("BP", "MF"))
#			getGoTerms(c("FOXM1", "TPT1"), c("BP", "MF"), c("TAS", "IC"))
#			getGoTerms(c("FOXM1", "TPT1"), c("BP", "MF"), c("TAS", "IC"), TRUE)
# *****************************************************************************
getGoTerms <- function(genes, ontology = NULL, evidenceType = NULL, def = FALSE){
	g = as.entrezId(genes)
	res = geneToGoMap[g]
	if (is.null(evidenceType))
		evidenceType = names(goEvidenceTypes)
	if (is.null(ontology))
		ontology = c("BP", "MF", "CC")
	res = lapply(res, function(goMap){
				x = goMap
				ind = (x$ONTOLOGY %in% ontology) & (x$EVIDENCE %in% evidenceType)
				cols = c("GO", "TERM", "ONTOLOGY", "EVIDENCE")
				if (def)
					cols = c(cols, "DEFINITION")
				return(x[ind, cols])
			})
	names(res) = genes
	return(res)
}


# *****************************************************************************
# Generates list object with one entry per gene. Each entry contains the GO
# terms for the corresponding gene.
#
# The goal of this method is to speed up the retrieval of GO terms for query 
# genes by improvinb on the slow performance of querying the "org.Hs.eg.db"
# annotation object via the standard "select()" interface. This goal is achieved
# by generating a list 'geneToGoMap' which is saved on a binary R file and is
# loaded everytime the "Utils.R" code is sourced. This object is then utilized
# by the getGoTerms() method, to run gene-based queries.
#
# The "standard" approach is to invoke this method without arguments. This will
# generate the 'geneToGoMap' object and store it in the file "geneToGoMap.rda",
# in the current working directory. After the operation is finished (the run 
# will take a good amount of time) the file "geneToGoMap.rda" must be copied
# the subdirectory "Data", under the directory where the "Utils.R" code resides.
# Notice that the method can also be invoked with the following arguments:
#
# ARGUMENTS
# * genes:	a vector of Entrez IDs or gene symbols. If provided the code is run
#			only for these genes, instead of all the genes.
# * save:	if FALSE, the resulting list object is not saved in a binary R file,
#			it is instead returned by the function.
# * fileName: the name of file where the binary R object will be stored. The 
#			default value, "geneToGoMap.rda", is the one expected by this code
#			when loading the binary object. So, don't change if planning to 
#			copy the file into the "Data" directory.
# *****************************************************************************
makeGeneToGoMap <- function(genes = NULL, save = TRUE, fileName = "geneToGoMap.rda"){
	if (is.null(genes))
		g = mappedkeys(org.Hs.egGO)

	else
		g = as.entrezId(genes)
	
	gTgMap = geneToGo(g)
	
	if (save){
		geneToGoMap = gTgMap
		save(list=c("geneToGoMap"), file=fileName)
		return()
	}
	
	return(gTgMap)
}


# *****************************************************************************
# Runs GO Term enrichment analysis for a set of query genes, using the topGO
# package in Bioconductor.
#
# ARGUMENTS:
# * genes:	the set of query genes. This can be:
# 		1. a vector of character strings representing Entrez IDs, e,g. c("1", "100", ..)
# 		2. a vector of intergers representing Entrez IDs, e.g., c(1, 100, ...)
# 		3. a vector of character strings representing gene symbols, e.g., c("FOXM1",
#	 		"TP53", ...)
# * ont:	the GO ontology against which to run the enrichment analysis, it
#			can be on of "BP" (biological process), "MF" (molecular function),
#			"CC" (cellular component).
# * universe:	the organism where the query genes come from.
# * method:	the enrichment method to use, see topGO documentation. The default
#			value, "weight01", corrects for the fact the gene membership in
#			various GO terms is not an independent event but rather depends on
#			the structure of the ontology.
# * test:	the statistical test of enrichment to use. See topGO documentation
#			for full details.
# * minSize:	The minimum allowed value length(genes). If the set of query
#			genes is less that this then no enrichment sore is calculated. See
#			topGO documentation for details.
# * correct:	Indicates if all enriched terms should be returned (correct = FALSE)
#			or only terms whose enrichment p-value is above an FDR theshold 
#			(correct = TRUE).
# * fdr:	FDR threshold to use when correct == TRUE.
#
# RETURN VALUE:
# A list res[[]] with 3 members:
# * res[[1]]:	a data frame listing the GO terms found to be enriched in query
#			genes, at the specified FDR level. Contains one row per enriched GO
#			term. Each row lists the GO term id, the enrichment p-value, the 
#			number of genes expected to be annotated to the term in a random
#			query gene set of size length(genes), the total number of genes from
#			the 'universe' organism annotated to that term, the number of genes 
#			from 'genes' annotated to that term, and the full descrption of the
#			term. NOTE: the 'universe' of genes against which the enrichment is
#			computed comprises all the genes from the 'universe' ogranimss that 
#			are annotation to at least one 'ont' ontology GO term.
# * res[[2]]:	the topGOData object from the topGO analysis - see package 
#			documentation for more details
# * res[[3]]:	the topGOresult object from the topGO analysis - see package 
#			documentation for more details
#
# We expect that res[[1]] will be the list element that is most useful as a 
# result value. However, the topGO package contains methods that can post-
# process the topGOdata and topGOresults objects, e.g., to generate images of
# the enriched terms against the ontology tree. Since re-generating these 
# objects from scratch is time consuming, we return than as elements of the 
# results list, in case the user would like to leverage them for that purpose.
# *****************************************************************************
goEnrichment <- function(genes, ont = "BP", universe = "human",  method = "weight01", 
		test = "fisher", minSize = 5, correct = TRUE, fdr = 0.05){
	
	# Return gracefully if no genes are provided
	if (is.null(genes) || length(genes) == 0)
		return(NULL)
	
	if (universe == "human"){
		all.genes = unique(unlist(annFUN.org(ont, mapping = "org.Hs.eg.db", ID = "entrez")))
		query.genes = as.entrezId(genes)
		geneList <- factor(as.integer(all.genes %in% query.genes))
		names(geneList) <- all.genes
		# If none of the query genes is annotated to a GO term, return an empty table
		# with the same column names lke a regular result, for cosnistency
		if (length(levels(geneList)) < 2){
			res = matrix(nrow = 0, ncol = 6)
			colnames(res) = c("GO.ID", "p-value", "Expected", "Annotated", "Significant", "Term")
			return(data.frame(res))
		}
		
		GOdata <- new("topGOdata", ontology = ont,	allGenes = geneList, nodeSize = minSize,
				annot = annFUN.org, mapping = "org.Hs.eg.db",ID = "entrez")	
	}
	
	GOres = runTest(GOdata, algorithm = method, statistic = test)
	summary <- GenTable(GOdata, GOres, topNodes = length(GOres@score), numChar = 100)
	# reduce results to GO terms passing Benjamini-Hochberg multiple hypothesis 
	# corrected pval <= 0.05, FDR <= 5%
	if (correct){
		summary[, "result1"] = gsub("<", "", summary[, "result1"])  # remove non-numeric characters
		summary <- summary[which(p.adjust(as.numeric(summary[,"result1"]),method="BH")<=fdr),]
	}
	# Fix the truncated term descriptions
#	if (nrow(summary) > 0)
#		for (i in 1:nrow(summary))
#			summary$Term[i] = goTermMap[summary$GO.ID[i]]
	# Rearrange columns, for convenience
	newOrder = c("GO.ID", "result1", "Expected", "Annotated", "Significant", "Term")
	summary = summary[, newOrder]
	colnames(summary)[2] = "p-value"
	return(list(summary, GOdata, GOres))
}


# *****************************************************************************
# Utility function used by 'makeGeneToGoMap()'. It does the actual mapping of
# genes to GO terms by extracting the key-value pairs from 'org.Hs.eg.db' and 
# joining with 'GO.db' to retrieve the detailed GO terms descriptions.
# *****************************************************************************
geneToGo <- function(genes){
	g = as.entrezId(genes)
	
	res = lapply(g, function(gid){
				k = c(gid)
				m = select(org.Hs.eg.db, keys=k, keytype="ENTREZID", columns=c("GO", "EVIDENCE", "ONTOLOGY"))[,-1]
				m1 = sapply(m[,"GO"], function(goid){
							k = c(goid)
							select(GO.db, keys=k, keytype="GOID", columns=c("TERM", "DEFINITION"))[, -1]
						})
				return(cbind(m, t(m1)))
			})
	names(res) = genes
	return(res)
}


# *****************************************************************************
# Convenience method. Takes as input any of the following vectors X:
# 1. a vector of character strings representing Entrez IDs, e,g. c("1", "100", ..)
# 2. a vector of intergers representing Entrez IDs, e.g., c(1, 100, ...)
# 3. a vector of character strings representing gene symbols, e.g., c("FOXM1",
#	 "TP53", ...)
#
# It returns a vector V of character strings representing Entrez IDs. For each
# of the 3 vectors X describe above, the vector V is as follows:
# 1. V[i] = X[i]
# 2. V[i] = toStrings(X[i])
# 3. V[i] = entrez ID correponding to gene symbol X[i]
# *****************************************************************************
as.entrezId <- function(genes){
	if (is.numeric(genes[1]))
		return(as.character(genes))
	else if (is.na(strtoi(genes[1])))
		return(as.character(geneSymbolToEntrezId(genes)))
	else
		return(genes)
}


# *********************************************************************************
# I did not write this code. The code is written by Davis McCarthy.
# It helps users to organize many ggplots in one figure. 
# You could add ggplots themselves consecutively or make a list and add it.
# You could put options like layout, colum numbers, file name, and title.
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
# *********************************************************************************
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL, title="") {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (nchar(title)>0){
    layout<-rbind(rep(0, ncol(layout)), layout)
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout), heights =if(nchar(title)>0){unit(c(0.5, rep(5,nrow(layout)-1)), "null")}else{unit(c(rep(5, nrow(layout))), "null")} )))
    
    # Make each plot, in the correct location
    if (nchar(title)>0){
      grid.text(title, vp = viewport(layout.pos.row = 1, layout.pos.col = 1:ncol(layout)))
    }
    
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


###
#   This function was downloaded from: https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R
#   Make a heatmap using gplots
#   Same as heatmap.2() of gplots, but this allows multiple ColSideColors and rowSideColors
###
heatmap.3 <- function(x,
                      Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
                      distfun = dist,
                      hclustfun = hclust,
                      dendrogram = c("both","row", "column", "none"),
                      symm = FALSE,
                      scale = c("none","row", "column"),
                      na.rm = TRUE,
                      revC = identical(Colv,"Rowv"),
                      add.expr,
                      breaks,
                      symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",
                      col = "heat.colors",
                      colsep,
                      rowsep,
                      sepcolor = "white",
                      sepwidth = c(0.05, 0.05),
                      cellnote,
                      notecex = 1,
                      notecol = "cyan",
                      na.color = par("bg"),
                      trace = c("none", "column","row", "both"),
                      tracecol = "cyan",
                      hline = median(breaks),
                      vline = median(breaks),
                      linecol = tracecol,
                      margins = c(5,5),
                      ColSideColors,
                      RowSideColors,
                      side.height.fraction=0.3,
                      cexRow = 0.2 + 1/log10(nr),
                      cexCol = 0.2 + 1/log10(nc),
                      labRow = NULL,
                      labCol = NULL,
                      key = TRUE,
                      keysize = 1.5,
                      density.info = c("none", "histogram", "density"),
                      denscol = tracecol,
                      symkey = max(x < 0, na.rm = TRUE) || symbreaks,
                      densadj = 0.25,
                      main = NULL,
                      xlab = NULL,
                      ylab = NULL,
                      lmat = NULL,
                      lhei = NULL,
                      lwid = NULL,
                      ColSideColorsSize = 1,
                      RowSideColorsSize = 1,
                      KeyValueName="Value",...){
  
  invalid <- function (x) {
    if (missing(x) || is.null(x) || length(x) == 0)
      return(TRUE)
    if (is.list(x))
      return(all(sapply(x, invalid)))
    else if (is.vector(x))
      return(all(is.na(x)))
    else return(FALSE)
  }
  
  x <- as.matrix(x)
  scale01 <- function(x, low = min(x), high = max(x)) {
    x <- (x - low)/(high - low)
    x
  }
  retval <- list()
  scale <- if (symm && missing(scale))
    "none"
  else match.arg(scale)
  dendrogram <- match.arg(dendrogram)
  trace <- match.arg(trace)
  density.info <- match.arg(density.info)
  if (length(col) == 1 && is.character(col))
    col <- get(col, mode = "function")
  if (!missing(breaks) && (scale != "none"))
    warning("Using scale=\"row\" or scale=\"column\" when breaks are",
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
  if (is.null(Rowv) || is.na(Rowv))
    Rowv <- FALSE
  if (is.null(Colv) || is.na(Colv))
    Colv <- FALSE
  else if (Colv == "Rowv" && !isTRUE(Rowv))
    Colv <- FALSE
  if (length(di <- dim(x)) != 2 || !is.numeric(x))
    stop("`x' must be a numeric matrix")
  nr <- di[1]
  nc <- di[2]
  if (nr <= 1 || nc <= 1)
    stop("`x' must have at least 2 rows and 2 columns")
  if (!is.numeric(margins) || length(margins) != 2)
    stop("`margins' must be a numeric vector of length 2")
  if (missing(cellnote))
    cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
  if (!inherits(Rowv, "dendrogram")) {
    if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
                                                 c("both", "row"))) {
      if (is.logical(Colv) && (Colv))
        dendrogram <- "column"
      else dedrogram <- "none"
      warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting row dendogram.")
    }
  }
  if (!inherits(Colv, "dendrogram")) {
    if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
                                                 c("both", "column"))) {
      if (is.logical(Rowv) && (Rowv))
        dendrogram <- "row"
      else dendrogram <- "none"
      warning("Discrepancy: Colv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting column dendogram.")
    }
  }
  if (inherits(Rowv, "dendrogram")) {
    ddr <- Rowv
    rowInd <- order.dendrogram(ddr)
  }
  else if (is.integer(Rowv)) {
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Rowv)) {
    Rowv <- rowMeans(x, na.rm = na.rm)
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else {
    rowInd <- nr:1
  }
  if (inherits(Colv, "dendrogram")) {
    ddc <- Colv
    colInd <- order.dendrogram(ddc)
  }
  else if (identical(Colv, "Rowv")) {
    if (nr != nc)
      stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
    if (exists("ddr")) {
      ddc <- ddr
      colInd <- order.dendrogram(ddc)
    }
    else colInd <- rowInd
  }
  else if (is.integer(Colv)) {
    hcc <- hclustfun(distfun(if (symm)
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Colv)) {
    Colv <- colMeans(x, na.rm = na.rm)
    hcc <- hclustfun(distfun(if (symm)
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else {
    colInd <- 1:nc
  }
  retval$rowInd <- rowInd
  retval$colInd <- colInd
  retval$call <- match.call()
  x <- x[rowInd, colInd]
  x.unscaled <- x
  cellnote <- cellnote[rowInd, colInd]
  if (is.null(labRow))
    labRow <- if (is.null(rownames(x)))
      (1:nr)[rowInd]
  else rownames(x)
  else labRow <- labRow[rowInd]
  if (is.null(labCol))
    labCol <- if (is.null(colnames(x)))
      (1:nc)[colInd]
  else colnames(x)
  else labCol <- labCol[colInd]
  if (scale == "row") {
    retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
    x <- sweep(x, 1, rm)
    retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
    x <- sweep(x, 1, sx, "/")
  }
  else if (scale == "column") {
    retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
    x <- sweep(x, 2, rm)
    retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
    x <- sweep(x, 2, sx, "/")
  }
  if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
    if (missing(col) || is.function(col))
      breaks <- 16
    else breaks <- length(col) + 1
  }
  if (length(breaks) == 1) {
    if (!symbreaks)
      breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                    length = breaks)
    else {
      extreme <- max(abs(x), na.rm = TRUE)
      breaks <- seq(-extreme, extreme, length = breaks)
    }
  }
  nbr <- length(breaks)
  ncol <- length(breaks) - 1
  if (class(col) == "function")
    col <- col(ncol)
  min.breaks <- min(breaks)
  max.breaks <- max(breaks)
  x[x < min.breaks] <- min.breaks
  x[x > max.breaks] <- max.breaks
  if (missing(lhei) || is.null(lhei))
    lhei <- c(keysize, 4)
  if (missing(lwid) || is.null(lwid))
    lwid <- c(keysize, 4)
  if (missing(lmat) || is.null(lmat)) {
    lmat <- rbind(4:3, 2:1)
    
    if (!missing(ColSideColors)) {
      #if (!is.matrix(ColSideColors))
      #stop("'ColSideColors' must be a matrix")
      if (!is.character(ColSideColors) || nrow(ColSideColors) != nc)
        stop("'ColSideColors' must be a matrix of nrow(x) rows")
      lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
      #lhei <- c(lhei[1], 0.2, lhei[2])
      lhei=c(lhei[1], side.height.fraction*ColSideColorsSize/2, lhei[2])
    }
    
    if (!missing(RowSideColors)) {
      #if (!is.matrix(RowSideColors))
      #stop("'RowSideColors' must be a matrix")
      if (!is.character(RowSideColors) || ncol(RowSideColors) != nr)
        stop("'RowSideColors' must be a matrix of ncol(x) columns")
      lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[,2] + 1)
      #lwid <- c(lwid[1], 0.2, lwid[2])
      lwid <- c(lwid[1], side.height.fraction*RowSideColorsSize/2, lwid[2])
    }
    lmat[is.na(lmat)] <- 0
  }
  
  if (length(lhei) != nrow(lmat))
    stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
  if (length(lwid) != ncol(lmat))
    stop("lwid must have length = ncol(lmat) =", ncol(lmat))
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  
  layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
  
  if (!missing(RowSideColors)) {
    if (!is.matrix(RowSideColors)){
      par(mar = c(margins[1], 0, 0, 0.5))
      image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
    } else {
      par(mar = c(margins[1], 0, 0, 0.5))
      rsc = t(RowSideColors[,rowInd, drop=F])
      rsc.colors = matrix()
      rsc.names = names(table(rsc))
      rsc.i = 1
      for (rsc.name in rsc.names) {
        rsc.colors[rsc.i] = rsc.name
        rsc[rsc == rsc.name] = rsc.i
        rsc.i = rsc.i + 1
      }
      rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
      image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
      if (length(rownames(RowSideColors)) > 0) {
        axis(1, 0:(dim(rsc)[2] - 1)/max(1,(dim(rsc)[2] - 1)), rownames(RowSideColors), las = 2, tick = FALSE)
      }
    }
  }
  
  if (!missing(ColSideColors)) {
    
    if (!is.matrix(ColSideColors)){
      par(mar = c(0.5, 0, 0, margins[2]))
      image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
    } else {
      par(mar = c(0.5, 0, 0, margins[2]))
      csc = ColSideColors[colInd, , drop=F]
      csc.colors = matrix()
      csc.names = names(table(csc))
      csc.i = 1
      for (csc.name in csc.names) {
        csc.colors[csc.i] = csc.name
        csc[csc == csc.name] = csc.i
        csc.i = csc.i + 1
      }
      csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
      image(csc, col = as.vector(csc.colors), axes = FALSE)
      if (length(colnames(ColSideColors)) > 0) {
        axis(2, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), colnames(ColSideColors), las = 2, tick = FALSE)
      }
    }
  }
  
  par(mar = c(margins[1], 0, 0, margins[2]))
  x <- t(x)
  cellnote <- t(cellnote)
  if (revC) {
    iy <- nr:1
    if (exists("ddr"))
      ddr <- rev(ddr)
    x <- x[, iy]
    cellnote <- cellnote[, iy]
  }
  else iy <- 1:nr
  image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
  retval$carpet <- x
  if (exists("ddr"))
    retval$rowDendrogram <- ddr
  if (exists("ddc"))
    retval$colDendrogram <- ddc
  retval$breaks <- breaks
  retval$col <- col
  if (!invalid(na.color) & any(is.na(x))) { # load library(gplots)
    mmat <- ifelse(is.na(x), 1, NA)
    image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
          col = na.color, add = TRUE)
  }
  axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
       cex.axis = cexCol)
  if (!is.null(xlab))
    mtext(xlab, side = 1, line = margins[1] - 1.25)
  axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
       cex.axis = cexRow)
  if (!is.null(ylab))
    mtext(ylab, side = 4, line = margins[2] - 1.25)
  if (!missing(add.expr))
    eval(substitute(add.expr))
  if (!missing(colsep))
    for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, length(csep)), xright = csep + 0.5 + sepwidth[1], ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
  if (!missing(rowsep))
    for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
  min.scale <- min(breaks)
  max.scale <- max(breaks)
  x.scaled <- scale01(t(x), min.scale, max.scale)
  if (trace %in% c("both", "column")) {
    retval$vline <- vline
    vline.vals <- scale01(vline, min.scale, max.scale)
    for (i in colInd) {
      if (!is.null(vline)) {
        abline(v = i - 0.5 + vline.vals, col = linecol,
               lty = 2)
      }
      xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
      xv <- c(xv[1], xv)
      yv <- 1:length(xv) - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (trace %in% c("both", "row")) {
    retval$hline <- hline
    hline.vals <- scale01(hline, min.scale, max.scale)
    for (i in rowInd) {
      if (!is.null(hline)) {
        abline(h = i + hline, col = linecol, lty = 2)
      }
      yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
      yv <- rev(c(yv[1], yv))
      xv <- length(yv):1 - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (!missing(cellnote))
    text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
         col = notecol, cex = notecex)
  par(mar = c(margins[1], 0, 0, 0))
  if (dendrogram %in% c("both", "row")) {
    plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
  }
  else plot.new()
  par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
  if (dendrogram %in% c("both", "column")) {
    plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
  }
  else plot.new()
  if (!is.null(main))
    title(main, cex.main = 1.5 * op[["cex.main"]])
  if (key) {
    par(mar = c(5, 4, 2, 1), cex = 0.75)
    tmpbreaks <- breaks
    if (symkey) {
      max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
      min.raw <- -max.raw
      tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
      tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
    }
    else {
      min.raw <- min(x, na.rm = TRUE)
      max.raw <- max(x, na.rm = TRUE)
    }
    
    z <- seq(min.raw, max.raw, length = length(col))
    image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
          xaxt = "n", yaxt = "n")
    par(usr = c(0, 1, 0, 1))
    lv <- pretty(breaks)
    xv <- scale01(as.numeric(lv), min.raw, max.raw)
    axis(1, at = xv, labels = lv)
    if (scale == "row")
      mtext(side = 1, "Row Z-Score", line = 2)
    else if (scale == "column")
      mtext(side = 1, "Column Z-Score", line = 2)
    else mtext(side = 1, KeyValueName, line = 2)
    if (density.info == "density") {
      dens <- density(x, adjust = densadj, na.rm = TRUE)
      omit <- dens$x < min(breaks) | dens$x > max(breaks)
      dens$x <- dens$x[-omit]
      dens$y <- dens$y[-omit]
      dens$x <- scale01(dens$x, min.raw, max.raw)
      lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
            lwd = 1)
      axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
      title("Color Key\nand Density Plot")
      par(cex = 0.5)
      mtext(side = 2, "Density", line = 2)
    }
    else if (density.info == "histogram") {
      h <- hist(x, plot = FALSE, breaks = breaks)
      hx <- scale01(breaks, min.raw, max.raw)
      hy <- c(h$counts, h$counts[length(h$counts)])
      lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
            col = denscol)
      axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
      title("Color Key\nand Histogram")
      par(cex = 0.5)
      mtext(side = 2, "Count", line = 2)
    }
    else title("Color Key")
  }
  else plot.new()
  retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
                                  high = retval$breaks[-1], color = retval$col)
  invisible(retval)
}

# ***************************************************************************************
# Pathway Analysis with TCGAbiolinks package
#
# Input: geneList       = a vector of gene symbols for pathway analysis [character vector]
#        title          = title of the pathway figure [character]
#        FDR_threshold  = pathway analysis FDR threshold (not DE analysis threshold) [numeric]
#        imgPrint       = print a plot of pathway analysis [TRUE/FALSE]
#        displayNum     = the number of pathways that will be displayed [numeric]
#                         (If there are many significant pathways show the few top pathways)
#        fName          = file path and file name of the output pathway figure [character]
#
# Output: Pathway analysis results in figure-using GO biological process and GO pathways
#         The return object is a list containing two data frames
#         [[1]]: GO biological process result
#         [[2]]: GO pathway result
#         if imgPrint=TRUE, then there will be also pathway result plot.
#         The left figure is from GO biological processes and the right figure is from
#         GO pathways. The x-axis represents -log10(FDR) of the pathways and the red line
#         indicates the ratio between the found hubs in a pathway and the total number of
#         genes in the pathway.
# ***************************************************************************************
pathwayAnalysis_TB <- function(geneList,
                               title="Pathway_Results",
                               FDR_threshold=0.05,
                               imgPrint=TRUE,
                               displayNum=Inf,
                               fName="pathway_results.pdf") {
  
  ### load library
  if(!require(TCGAbiolinks)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("TCGAbiolinks")
    library(TCGAbiolinks)
  }
  
  ### get pathway results from the gene list
  ansEA <- TCGAanalyze_EAcomplete(TFname=title, as.character(geneList))
  
  ### if no pathways were found - return "Result does not exist"
  if(is.null(ncol(ansEA$ResBP)) && is.null(ncol(ansEA$ResPat))) {
    writeLines("Result does not exist")
  } else {
    
    ### FDR with zero or below than 1.00E-10 will be changed to 1.00E-10
    ### FDR with zero will make an error in the print plot function
    ### I personally fixed the bug in TCGAbiolinks package
    fdrZeroToMin <- function(ansEA_result) {
      cnt <- 0
      
      for(i in 1:ncol(ansEA_result)) {
        Go_temp <- strsplit(ansEA_result[1,i], "; ", fixed = TRUE)[[1]]
        Go_temp2 <- strsplit(Go_temp[2], "e", fixed = TRUE)[[1]]
        if(grepl("FDR= 0.00e", ansEA_result[1,i])) {
          cnt <- cnt+1
        } else if((!is.na(Go_temp2[2])) && as.numeric(Go_temp2[2]) <= -10){
          cnt <- cnt+1
        } else {
          break
        }
      }
      
      if(cnt > 0) {
        for(i in 1:cnt) {
          Go_temp <- strsplit(ansEA_result[1,i], "; ", fixed = TRUE)[[1]]
          Go_temp[2] <- "FDR= 1.00e-10"
          ansEA_result[1,i] <- paste(Go_temp, collapse = "; ")
        }
      }
      
      return(ansEA_result)
    }
    
    filterWithFDR <- function(ansEA_result) {
      idx <- NULL
      
      for(i in 1:ncol(ansEA_result)) {
        Go_temp <- strsplit(ansEA_result[1,i], "; ", fixed = TRUE)[[1]]
        fdr <- as.numeric(substr(Go_temp[2], 6, nchar(Go_temp[2])))
        
        if(fdr >= FDR_threshold) {
          idx <- c(idx, i)
        }
      }
      
      if(length(idx) > 0) {
        return(ansEA_result[,-idx])
      } else {
        return(ansEA_result)
      }
    }
    
    if(imgPrint == TRUE) {
      ### Normally, I would like to show both results from GO biological processes and GO pathways
      ### When if it is impossible, only show one
      if(is.null(ncol(ansEA$ResBP))) {
        TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResPat), PathTab = fdrZeroToMin(filterWithFDR(ansEA$ResPat)), nRGTab = geneList, nBar = displayNum, mfrow = c(1,1), filename = fName, text.size = 2)
      } else if(is.null(ncol(ansEA$ResPat))) {
        TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP), GOBPTab = fdrZeroToMin(filterWithFDR(ansEA$ResBP)), nRGTab = geneList, nBar = displayNum, mfrow = c(1,1), filename = fName, text.size = 2)
      } else {
        TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP), GOBPTab = fdrZeroToMin(filterWithFDR(ansEA$ResBP)), PathTab = fdrZeroToMin(filterWithFDR(ansEA$ResPat)), nRGTab = geneList, nBar = displayNum, mfrow = c(1,2), filename = fName, text.size = 2, xlim = c(0,10))
      }
    }
    
    return(list(filterWithFDR(ansEA$ResBP), filterWithFDR(ansEA$ResPat)))
  }
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
                               org,
                               database,
                               title="Pathway_Results",
                               pv_threshold=0.05,
                               displayNum=Inf,
                               imgPrint=TRUE,
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


# ******************************************************************************************
# A function to get genomic location information from a vector of gene names
#
# The function takes input as a vector of genes and returns a matrix with
# length(genes) rows and 4 columns:
# "Chr": an integer in the range 1-23, representing the chromosome where the gene is located
#        (you can map both X and Y chromosomes to 23)
# "First": an integer representing the genomic location of the gene's start site
#        on the chromosome. In case a gene has multiple start sites, use the earliest one
#        (for gens on the FWD strand this is the one with the smallest absolute coordinate
#        while for genes on the REV strand this is the one with the highest absolute coordinate)
# "Last": an integer representing the genomic location of the gene's end site
#        on the chromosome. In case a gene has multiple end sites, use the earliest one
#        (for gens on the FWD strand this is the one with the smallest absolute coordinate
#        while for genes on the REV strand this is the one with the highest absolute coordinate)
# "Strand": ether 1 or 2, representing if the gene is on the forward or reverse strand
#
# The input vector can be either gene symbol or Entrez ID
# The rows in the results matrix should be named after the values in "genes" and the columns
# specified above
# ******************************************************************************************
getGenomicLocation <- function(genes) {
  if(is.na(genes) || is.null(genes)) {
    return(NA)
  } else {
    ### get character vector of the genes in Entrez ID
    g <- as.entrezId(genes)
    
    ### make a map of Entrez_ID - Genomic location info
    map <- org.Hs.egCHRLOC
    mapped_genes <- mappedkeys(map)
    map <- as.list(map[mapped_genes])
    map2 <- org.Hs.egCHRLOCEND
    mapped_genes2 <- mappedkeys(map2)
    map2 <- as.list(map2[mapped_genes2])
    
    ### make the result
    result <- sapply(g, function(x) {
      ### get location info
      info <- map[[x]]
      info2 <- map2[[x]]
      
      if(is.na(info) || is.null(info)) {
        return(c(NA, NA, NA))
      } else {
        ### remove alternative haplotype chromosomes
        info <- info[which(nchar(names(info)) < 3)]
        info2 <- info2[which(nchar(names(info2)) < 3)]
        
        ### For genes on the FWD strand = one with the smallest absolute coordinate
        ### For genes on the REV strand = one with the highest absolute coordinate
        info <- info[which(info == min(info))]
        info2 <- info2[which(info2 == min(info2))]
        
        ### Chr
        c <- names(info[1])
        if((c == "X") || (c == "Y")) {
          r <- 23
        } else {
          r <- c
        }
        
        ### First
        r <- c(r, abs(info[1]))
        
        ### Last
        r <- c(r, abs(info2[1]))
        
        ### Strand
        if(info[1] > 0) {
          r <- c(r, 1)
        } else {
          r <- c(r, 2)
        }
        
        
        return(as.numeric(r))
      }
    })
    
    ### organize the result
    result <- t(result)
    rownames(result) <- as.character(genes)
    colnames(result) <- c("Chr", "First", "Last", "Strand")
    
    return(result)
  }
}


#******************************************************************************************
# Replace NA values in a matrix.
#
# ARGUMENTS:
# * mat:	the matrix to be processed
# * val:	the value that will replace the NAs
#
# RETURNS
# A matrix "res" of the same dimensions as "mat" where:
#		res[i,j] = mat[i,j]			if mat[i,j] is not NA
#		res[i,j] = val				if mat[i,j] is NA
# ******************************************************************************************
replace.NA <- function(mat, val){
	x = as.numeric(mat)
	ind = which(is.na(x))
	x[ind] = val
	res = matrix(x, nrow(mat), ncol(mat))
	rownames(res) = rownames(mat)
	colnames(res) = colnames(mat)
	return(res)
}


#******************************************************************************************
# Generate MDS plot for matrix columns.
#
# * mat:		matrix whose columns are to be clustered.
# * plot_names:	if TRUE, name the points on the MDS plot.
# * alt_names:	String vector. By default, if plot_names == TRUE, the plotted points are 
#		named using colnames(mat). However, if alt_names in not NULL, then this is used
#		instead. In that case, length(alt_names) should be equal to ncol(mat) and the point
#		mat[, i] will be named after alt_names[i].
# * groups:		if not NULL, it must be a string vector such that names(groups) == colnames(mat)
#		and groups[i] is the name of the group where the i-th column of "mat" belongs. This
#		info will be used to color members of the same group using the same color. If the
#		value of this argument is NULL, then no group-based coloring is performed.
# * dist_fun:	distance function to use for computing distances between the column vectors in
#		mat. This should take as input a matrix object and return an object of class "dist". If
#		the value of this argument is NULL, then the value of the argument dist_options below
# 		is used for determining how to compute the distances.
# * dist_options:	if dist_fun is NULL the the standard "dist" function in R is used for 
#		computing distances. In that case, the value of dist_options specifies which distance 
#		option to use when calling "dist". The default option is "euclidan".
# * save:		if TRUE, save plot to file. Otherwise just plot on screen.
# * f_name:		if save == TRUE, this is the full pathname of the file were the plot will be saved.
# * width, height, res:	values of graphical parameters to use when generating the plot.
# * xlab, ylab, main:	titles for x-axis, y-axis, and entire plot, respectively 
mdsPlot <- function(mat, plot_names = FALSE, alt_names = NULL, groups = NULL, dist_fun = NULL,
		dist_options = c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman", "kendall"), 
		save = FALSE, f_name = "./mds_plot.png", width = 1000, height = 1000, res = 130,
		xlab = "", ylab="", main=""){
	
  ### load library
  if(!require(ArgumentCheck)) {
    install.packages("ArgumentCheck")
    library(ArgumentCheck)
  }
  
  ### argument checking
  check <- ArgumentCheck::newArgCheck()
  if(is.null(mat)) {
    ArgumentCheck::addError(
      msg = "[mat]: matrix whose columns are to be clustered should exists",
      argcheck = check
    )
  }
  if((!is.null(alt_names)) && (!length(alt_names) == ncol(mat))) {
    ArgumentCheck::addError(
      msg = "[alt_names]: length(alt_names) should be equal to ncol(mat)",
      argcheck = check
    )
  }
  if((!is.null(groups)) && (!length(groups) == ncol(mat))) {
    ArgumentCheck::addError(
      msg = "[groups]: length(groups) should be equal to ncol(mat)",
      argcheck = check
    )
  }
  if((!is.null(dist_fun)) && (!class(dist_fun(mat)) == "dist")) {
    ArgumentCheck::addError(
      msg = "[dist_fun]: it should be a function which generates \"dist\" object",
      argcheck = check
    )
  }
  if((is.null(dist_fun)) && (is.null(dist_options))) {
    ArgumentCheck::addError(
      msg = "[dist_options]: if the \"[dist_fun]\" was not provided, \"[dist_options]\" should be provided and should be one of \"euclidean\", \"maximum\", \"manhattan\", \"canberra\", \"binary\", \"minkowski\", \"pearson\", \"spearman\", or \"kendall\"",
      argcheck = check
    )
  } else if((is.null(dist_fun)) && (!dist_options %in% c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman", "kendall"))) {
    ArgumentCheck::addError(
      msg = "[dist_options]: if the \"[dist_fun]\" was not provided, \"[dist_options]\" should be provided and should be one of \"euclidean\", \"maximum\", \"manhattan\", \"canberra\", \"binary\", \"minkowski\", \"pearson\", \"spearman\", or \"kendall\"",
      argcheck = check
    )
  }
  if((save == TRUE) && is.null(f_name)) {
    ArgumentCheck::addError(
      msg = "[f_name]: if \"[save]\" == TRUE, then \"[f_name]\" should be provided",
      argcheck = check
    )
  }
  ArgumentCheck::finishArgCheck(check)
  
  
  ### save the mds plot as a png format
	if(save){
		png(f_name, width = width, height = height, res = res)
	}
	
  ### make a distance matrix
  if(!is.null(dist_fun)) {
    d <- dist_fun(mat)
  } else if(dist_options[1] == "euclidean") {
    d <- dist(t(mat), method = "euclidean")
  } else if(dist_options[1] == "maximum") {
    d <- dist(t(mat), method = "maximum")
  } else if(dist_options[1] == "manhattan") {
    d <- dist(t(mat), method = "manhattan")
  } else if(dist_options[1] == "canberra") {
    d <- dist(t(mat), method = "canberra")
  } else if(dist_options[1] == "binary") {
    d <- dist(t(mat), method = "binary")
  } else if(dist_options[1] == "minkowski") {
    d <- dist(t(mat), method = "minkowski")
  } else if(dist_options[1] == "pearson") {
    d <- as.dist(1-cor(mat, method = "pearson"))
  } else if(dist_options[1] == "spearman") {
    d <- as.dist(1-cor(mat, method = "spearman"))
  } else if(dist_options[1] == "kendall") {
    d <- as.dist(1-cor(mat, method = "kendall"))
  }
	
  ### get MDS points
  fit <- cmdscale(d,eig=TRUE, k=2)
  Dimension1 <- fit$points[,1]
  Dimension2 <- fit$points[,2]
  
  ### group coloring
  if(is.null(groups)) {
    ### make a MDS plot
    plot(Dimension1, Dimension2, main=main, xlab=xlab, ylab=ylab)
    
    ### print point names
    if(plot_names == TRUE) {
      if(is.null(alt_names)) {
        text(Dimension1, Dimension2, labels = labels(d), cex=.7, pos=3)
      } else {
        text(Dimension1, Dimension2, labels = alt_names, cex=.7, pos=3)
      }
    }  
  } else {
    ###
	if (length(unique(as.character(groups))) < 6)
		colors = c("black", "red", "blue", "magenda", "green")[1:length(unique(as.character(groups)))]
	else
		colors = rainbow(length(unique(as.character(groups))))
    names(colors) = unique(as.character(groups))
    
    ### make a MDS plot
    plot(Dimension1, Dimension2, main=main, xlab=xlab, ylab=ylab,
         col = colors[as.character(groups)])
    legend("topright", legend = unique(as.character(groups)),
           col = colors[unique(as.character(groups))], pch = 15,
           title = "Sample Groups", cex = 0.7)
    
    ### print point names
    if(plot_names == TRUE) {
      if(is.null(alt_names)) {
        text(Dimension1, Dimension2, labels = labels(d), cex=0.7, pos=3, col = colors[as.character(groups)])
      } else {
        text(Dimension1, Dimension2, labels = alt_names, cex=0.7, pos=3, col = colors[as.character(groups)])
      }
    }
  }
  
  ### print out the plot
	if(save)
		dev.off()
  
}


# ******************************************************************************
# Method to read RNA-seq read counts
#
# ARGUMENTS
# * file_names:		character vector; each entry is the full pathname to a 
#		tab-delimited text file containing read count data, in the usual format: 
#		one column for each of the N samples in the file, the first column 
#		contains gene symbols or gene ids, column headers are sample names. The
#		header line may contains N or (N+1) entries. If the former, each entry
#		is a sample name. If the latter, the first header entry is assummed to
#		be the column name for the gene column (no need to retain it, but the 
#		code must be prepared to handle this dual possibility for the header line).
# * entrez_ids:		if the variable is FALSE, the gene identifiers found in the 
#		gene column are used without change. Otherwise, these identifiers can be
#		assumed to be gene symbols and must be converted to entrez Ids by using 
#		the method Utils::geneSymbolToEntrezId(). If a gene symbol cannot be 
# 	mapped to an entrez Id, then the data matrix row correponding to that gene
#		is removed. If mulitple symbols map to the same entrez id, then the value of
#		the argument "combine" (see below) specifies how to combine the rows 
#		corresponsing to these symbols:
# * combine:		specifes how to handle a situation where multiple gene identifiers
#		g1, g2, ..., gk map to the same gene (here we assume that gi are listed in the
#		order in which they are encountered in the first column of the data file). 
#		Acceptable values and their meaning are:
#		- "first":	keep only the data row g1 (the first gene) and remove all other rows.
#		- "sum":	create one aggregate row by summing the reads for all gi in each sample.
#		- "min":	keep the minimum read count among all gi within each sample.
#		- "max":	keep the maximum read count among all gi within each sample.
#		- "average":average the read counts for all gi within each sample and round 
#					to the closest integer.
#		- "abort":	write out on the console an error message stating that there are duplicate
#					gene identifiers and return NULL.
# * merge:	boolean value specifying if data matrices should be combined - relevant
#		only when "file_names" contains more than one entry. If TRUE, all matrices are
#		merged into one matrix. Otherwise, the function returns a separate matrix for each
#		file in "file_names". If merge == TRUE and also combine == TRUE, the combination of
#		gene identifiers is performed first for each file in "file_names" separately, followed
#		by the data matrix merging operation. If the same sample name appears in more than one
#		input file, then write out on the console an error message stating that the are duplicate
#		sample names and return NULL.
# * merge_mode:		a character vector, specifying how data matrices should be merged (used
#		when "merge" == TRUE). Possible values are:
#		- "shared": only retains gene identifiers that are common among all individual data
#					matrices.
#		- "all":	retains all identifiers encountered across all data matrices. If an
#					identifier G is not present in matrix D, then, in the merged matrix, the
#					value of G in all samples in D is set to the value of argument 
#					"missing_values".
# * missing_values:	an integer or NA. Used when merge == TRUE and merge_mode == "all", as 
#		described above. 
#
# RETURN VALUE
# If "file_names" contains exactly one file name, the return value is a data matrix M x N
# containing the read counts for each (gene, sample) combination, with one column per sample 
# and one row per gene. Columns are named with the sample names found in the header of the input 
# data file. If entrez_ids == FALSE, the row names are the gene identifiers found in the first column of
# the data file. If entez_ids == TRUE, row names are entrez ids, derived as described above.
# 
# If "file_names" contains more than one entry, then the value of the argument "merge" determines
# what the function returns. If "merge" == FALSE the function returns a list with one entry
# for each file in "file_names". The i-th entry correposponds to the i-th file name and contains
# a read count data matrix generated by parsing that file, as described in the previous paragraph. 
# The list is named and the name of the i-th entry is the file name stem of the i-th entry in 
# file_names, after the file extension has been removed. E.g., if:
#		file_names[i] = "/ifs/scratch/af_lab/data/ranseqfile.txt"
# then
#		names(results_object)[i] = "ranseqfile"
# If merge == TRUE, then the function returns a data matrix with columns correponding to the samples
# in all files in "file_names" and rows corresponding to gene identifiers generated according to the 
# values of the arguments "entrez_ids", "combine", and "merge_mode", as describe above. Column names 
# are sample names and row names are gene identifiers, the same as in the case where length(file_names) == 1.

readRNAseqData <- function(file_names,
                           entrez_ids = FALSE,
                           combine = c("first", "sum", "min", "max", "average", "abort"),
	                         merge = FALSE,
	                         merge_mode = c("shared", "all"),
	                         missing_values = NA) {
	
  ### load library
  if(!require(checkmate)) {
    install.packages("checkmate")
    library(checkmate)
  }
  
  ### argument checking
  assertCharacter(file_names)
  assertFlag(entrez_ids)
  assertChoice(combine[1], c("first", "sum", "min", "max", "average", "abort"))
  assertFlag(merge)
  assertChoice(merge_mode[1], c("shared", "all"))
  assertIntegerish(missing_values)
  
  ### make an empty list for saving the files
  file_list <- vector("list", length(file_names))
  
  ### load the files
  for(i in 1:length(file_names)) {
    ### read each file from the file_names
    file_list[[i]] <- read.table(file = file_names[i], header = TRUE, sep = "\t",
                                 stringsAsFactors = FALSE, check.names = FALSE)
    
    ### if there is no header line for the gene column, then it automatically
    ### set the first column as rownames. Copy the rownames to the first column
    ### so that the further analysis can be done accurately.
    if(!identical(rownames(file_list[[i]]), as.character(1:nrow(file_list[[i]])))) {
      file_list[[i]] <- data.frame(Gene=rownames(file_list[[i]]), file_list[[i]],
                                   stringsAsFactors = FALSE, check.names = FALSE)
    }
    
    ### if the found gene identifier is gene symbol and we want to convert it to entrez ids
    if(entrez_ids) {
      ### get Entrez IDs for the corresponding gene symbols
      eIDs <- geneSymbolToEntrezId(file_list[[i]][,1])
      
      ### if the column is not gene symbol, there would be no entrez ids returned
      if(is.null(eIDs) || is.na(eIDs)) {
        stop("ERROR: The first column is not GENE SYMBOL")
      }
      
      ### remove genes that do not have entrez ids
      file_list[[i]] <- file_list[[i]][which(file_list[[i]][,1] %in% names(eIDs)),]
      
      ### change the the first column to entrez ids
      file_list[[i]][,1] <- eIDs[file_list[[i]][,1]]
    }
    
    ### if there are other character columns, remove them
    ### keep only the raw counts except the first (entrez id) column
    chr_idx <- which(sapply(file_list[[i]][-1], is.character))
    if(length(chr_idx) > 0) {
      file_list[[i]] <- file_list[[i]][,-(chr_idx+1)]
    }
    
    ### if there are duplicated gene identifiers
    if(length(file_list[[i]][,1]) != length(unique(file_list[[i]][,1]))) {
      ### get duplicated indices (this does not retain the first ones)
      dups <- which(duplicated(file_list[[i]][,1]))
      
      ### get the indicies of first duplicated elements
      first_dups <- setdiff(which(duplicated(file_list[[i]][,1], fromLast = TRUE)), dups)
      
      ### handle the case based on the combine parameter
      if(combine[1] == "first") {
        ### there is nothing to do here if you want to keep only the first thing
      } else if(combine[1] == "sum") {
        ### sum up the duplicated rows and save them to the first dups
        for(j in 1:length(first_dups)) {
          file_list[[i]][first_dups[j],2:ncol(file_list[[i]])] <-
            apply(file_list[[i]][which(file_list[[i]][,1] == file_list[[i]][first_dups[j],1]),
                                 2:ncol(file_list[[i]])], 2, sum)
        }
      } else if(combine[1] == "min") {
        ### get min values of the duplicated rows by each sample and save them to the first dups
        for(j in 1:length(first_dups)) {
          file_list[[i]][first_dups[j],2:ncol(file_list[[i]])] <-
            apply(file_list[[i]][which(file_list[[i]][,1] == file_list[[i]][first_dups[j],1]),
                                 2:ncol(file_list[[i]])], 2, min)
        }
      } else if(combine[1] == "max") {
        ### get max values of the duplicated rows by each sample and save them to the first dups
        for(j in 1:length(first_dups)) {
          file_list[[i]][first_dups[j],2:ncol(file_list[[i]])] <-
            apply(file_list[[i]][which(file_list[[i]][,1] == file_list[[i]][first_dups[j],1]),
                                 2:ncol(file_list[[i]])], 2, max)
        }
      } else if(combine[1] == "average") {
        ### get rounded average values of the duplicated rows by each sample and save them to the first dups
        for(j in 1:length(first_dups)) {
          file_list[[i]][first_dups[j],2:ncol(file_list[[i]])] <-
            apply(file_list[[i]][which(file_list[[i]][,1] == file_list[[i]][first_dups[j],1]),
                                 2:ncol(file_list[[i]])], 2, function(x) round(mean(x)))
        }
      } else if(combine[1] == "abort") {
        stop("ERROR: There are duplicated GENE SYMBOLs")
      }
      
      ### remove the duplicates and only keep the first ones
      file_list[[i]] <- file_list[[i]][-dups,]
    }
    
    ### set the row names with the first column and remove it
    rownames(file_list[[i]]) <- file_list[[i]][,1]
    file_list[[i]] <- file_list[[i]][,-1]
    
    ### change data.frame to matrix
    file_list[[i]] <- as.matrix(file_list[[i]])
  }
  
  ### handle the merge case
  if(merge) {
    ### if there are duplicated sample names, print out an error message and stop the process
    if(length(Reduce(union, lapply(file_list, colnames))) != length(unlist(lapply(file_list, colnames)))) {
      stop("ERROR: There are duplicated sample names and you are trying to merge them")
    }
    
    if(merge_mode[1] == "shared") {
      file_list <- Reduce(
        function(df1, df2) {
          temp <- merge(df1, df2, by = "row.names", all = FALSE, sort = FALSE)
          rownames(temp) <- temp[,1]
          return(temp[,-1])
        }, file_list)
    } else if(merge_mode[1] == "all") {
      file_list <- Reduce(
        function(df1, df2) {
          temp <- merge(df1, df2, by = "row.names", all = TRUE, sort = FALSE)
          rownames(temp) <- temp[,1]
          return(temp[,-1])
        }, file_list)
      
      ### handle missing values
      file_list[is.na(file_list)] <- missing_values
    }
  }
  
  ### organize the result for returning
  if(is(file_list, 'list') && length(file_list) == 1) {
    file_list <- file_list[[1]]
  } else if(!merge) {
    names(file_list) <- sapply(basename(file_names), function(x) {
      temp <- strsplit(x, ".", fixed = TRUE)[[1]]
      return(paste(temp[1:length(temp)-1], collapse = "."))
    })
  }
  
  return(file_list)
}

####################################################
### A function to perform DE analysis with limma ###
####################################################
#' @title limmaWithComparisons
#' @param normCnt normalized count matrix
#' @param grp a character vector of class info of the samples
#' @param exp_class a string of the experiment group's name
#' @param ctrl_class a string of the control group's name
#' @param bat_eff a character vector of batch effect info of the samples
#' @return data.frame
#' @export
#' @author Hyunjin Kim
####################################################
limmaWithComparisons <- function(normCnt, grp, exp_class, ctrl_class, bat_eff=NULL) {
  
  ### load library
  if(!require(limma)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("limma")
    library(limma)
  }
  
  ### sometimes, there are some variables which can not be transformed into R variable names
  ### so, just to be safe, change all the variables to R-usuable ones
  grp <- make.names(grp)
  exp_class <- make.names(exp_class)
  ctrl_class <- make.names(ctrl_class)
  
  ### make a design matrix for DE analysis
  sampleType <- relevel(as.factor(grp), ref = ctrl_class)
  if(is.null(bat_eff)) {
    design <- model.matrix(~0+sampleType)
    colnames(design) <- levels(sampleType)
  } else {
    bat_eff <- make.names(bat_eff)
    bat_eff <- as.factor(bat_eff)
    design <- model.matrix(~0+sampleType+bat_eff)
    colnames(design) <- c(levels(sampleType), levels(bat_eff)[-1])
  }
  
  ### fir the linear model
  fit <- lmFit(normCnt, design)
  
  ### extract specific comparison of interest
  contrastMat <- makeContrasts(contrasts=paste(exp_class,ctrl_class,sep="-"), levels=design)
  
  ### fit the contrasts
  fit2 <- contrasts.fit(fit, contrastMat)
  fit2 <- eBayes(fit2)
  
  ### get the differentially expressed genes
  result <- topTable(fit2, adjust.method="BH", number=Inf)
  
  ### order based on adj.p.val
  result <- result[order(result$adj.P.Val),]
  
  return(result)
}

