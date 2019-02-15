# *****************************************************************************
# Run GO enrichment analysis on the hubs of an interactome, in chunks. The code 
# is invoked with the following command line parameters and is designed for 
# cluster array execution (see the master submission script):
#	<tcga_acro> <GOontology> <ChunkSize> <Increment>
# where:
# * tcga_acro:	tcga acronym (e.g., brca, gbm, ov) indicating which interactome
#		to run the analysis on. NOTE: the code expects that a binary R file named
#				<tcga_acro>.rda
#		will be available in the directory where the sipts is invoked form and will
#		contain a variable nasmed <tcga_acro> containing the query interactome.
# * GOontology:	either of the strings BP, MF, CC. Indicates which GO ontology to
#		use.
# * ChunkSize:	An integer used for incremental processing, see next argument.
# * Increment: An integer used to specify the range of hub genes to process.
#		E.g., say that tcga_acro = brca. Then for a Increment value of N, the 
#		code will process the regulons for the hubs brca[[2]][start:end] where: 
#				start = (N-1)*ChunkSize+1 
#				end   =  N*ChunkSize.
# 
# The analysis generateds a named list with entries named after processed hub genes.
# The entry for hub gene G is the data frame Utils::goEnrichment(regulon(G))[[1]].
# This list is assigned to a variable named:
#	<tcga_acro>Res_<start>           (e.g., brcaRes_1) 
# and is saved for post-processing to a binary R file with the same name.
# *****************************************************************************


UTILSLIB = "~/cvs/R/Common/Utils.R"

# *****************************************************************************
# Run the GO enrichment analysis. Keeps track of the execution time, for
# logging purposes
# 
# ARGUMENTS:
# * net:		the slice of the interactome on which to run the analysis.
# * goOnt:		string specifying which GO ontology to use (BP, MF, CC).
# * logFile:	file for storing logging info.
# *****************************************************************************
doEnrichment <- function(net, goOnt = "BP", logFile = NULL){
	start = Sys.time()
	source(UTILSLIB, chdir = TRUE)
	res = lapply(net, function(x, ontology){
				return(goEnrichment(abs(x[,1]), ont = ontology)[[1]])
			}, ontology = goOnt)
	end = Sys.time()
	names(res) = names(net)
	if (!is.null(logFile))
		cat(paste("\n\n\tElapsed time -> ", paste(toString(end-start), attr( end-start, "units"))), file = logFile, append=TRUE, sep="\n")
	return(res)
}

# Parse command-line arguments
args = commandArgs(TRUE)
if (length(args) != 4)
	stop("use:\t<command_name> tcga_acro GOontology ChunkSize Increment")
load(paste(args[1], ".rda", sep=""))
net = get(args[1])
ontology = args[2]
INC = as.integer(args[3])
start = (as.integer(args[4]) - 1)*INC+1
L = length(net[[2]])	

# Invoke the analysis on the specified slice of the interactome and store the 
# result list in a binary R file.
if (start <= L){
	end = min(L, start+INC-1)
	res = doEnrichment(net[[2]][start:end], goOnt = ontology, paste(args[1], "Log.txt", sep=""))
	resVarName = paste(args[1], "Res_", start, sep="")
	fileName = paste(args[1], "Res_", start, ".rda", sep="")
	assign(resVarName, res, envir = globalenv())
	save(list = c(resVarName), file=fileName)
}

