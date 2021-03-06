#!/usr/bin/env Rscript

# Takes as input one or more 4-col network files generated by running ARACNe 
# on a gene expression data set. Using the *same* data set, it converts them 
# to the CPTAC 6-col format. Produces two files, the first is the text version
# of the 6-col network and the second is a binary R file containing the 
# corresponding 'regulon' object, the format needed for downstream analyses by
# the viper package.
#
# ARGUMENTS
# * expmat.file:	location of original gene expression file used for generating  
#					the ARACNe networks.
# * netFiles:		location of the 4-column formatted ARACNe network text files
# * resultDir:		directory where to store the output files
# * resultName:		text to use for building the file names of the result files
consolidateNets <- function(expmat.file = "./gep.txt", 
		netFiles = c("./tf_run/finalNetwork_4col.tsv", "./cotf_run/finalNetwork_4col.tsv", 
				"./signal_run/finalNetwork_4col.tsv"),
		resultDir = "./", resultName = "finalNetwork"){
	
	library(viper)
  resultDir = paste(resultDir, "/", sep = "") # to be safe...
		
# Turn the 4-col files network files into one merged 3-col file, as needed by the method
# aracne2regulon, and store it in a temporary directory
	temp_dir = paste(resultDir, "temp_", as.integer(Sys.time()), "/", sep="")
	dir.create(temp_dir)
	temp_combined = paste(temp_dir, "consolidatedNet.tsv", sep="") # Put all interacions here
	master_mat = NULL
	for (i in 1:length(netFiles)){
		mat = as.matrix(read.table(netFiles[i]))
		master_mat = rbind(master_mat, mat)
	}
	
	write.table(master_mat[, 1:3], file = temp_combined, row.names = FALSE, col.names = FALSE, sep="\t")
	
# Run the method that will create the 'regulon' object
	regulon <- aracne2regulon(temp_combined, expmat.file, format="3col")
	
# order the regulon targets in increasing order of their Entrez id
	for (i in 1:length(regulon)){
		x = as.integer(names(regulon[[i]]$tfmode))
		regulon[[i]]$tfmode = regulon[[i]]$tfmode[order(x)]
		regulon[[i]]$likelihood = regulon[[i]]$likelihood[order(x)]
	}
	regulon = regulon[order(as.integer(names(regulon)))]
	
# Prepare the final 6-column matrix and order the rows first by the hub Entrez gene id
# and then by the Entrez target gene id.
	int_no = sum(sapply(regulon, function(x){return(length(x$tfmode))})) # total num of interactions in regulon
	final_mat = matrix(nrow = int_no, ncol=6)
	colnames(final_mat) = c("Hub", "Target", "MI", "MoA", "likelihood", "pvalue")
	mi_mat = matrix(nrow = length(unique(master_mat[,1])), ncol = length(unique(master_mat[,2])))
	rownames(mi_mat) = unique(master_mat[,1])
	colnames(mi_mat) = unique(master_mat[,2])
	p_value_mat = matrix(nrow = length(unique(master_mat[,1])), ncol = length(unique(master_mat[,2])))
	rownames(p_value_mat) = unique(master_mat[,1])
	colnames(p_value_mat) = unique(master_mat[,2])
	
	for (i in 1:nrow(master_mat)){
		mi_mat[toString(master_mat[i,1]), toString(master_mat[i,2])] = master_mat[i,3]
		p_value_mat[toString(master_mat[i,1]), toString(master_mat[i,2])] = master_mat[i,4]
	}
	
	pos = 1
	for (i in 1:length(names(regulon))){
		for (j in 1:length(regulon[[i]]$tfmode)){
			final_mat[pos, 1] = names(regulon)[i]
			final_mat[pos, 2] = names(regulon[[i]]$tfmode)[j]
			final_mat[pos, 3] = mi_mat[names(regulon)[i], names(regulon[[i]]$tfmode)[j]]
			final_mat[pos, 4] = regulon[[i]]$tfmode[j]
			final_mat[pos, 5] = regulon[[i]]$likelihood[j]
			final_mat[pos, 6] = p_value_mat[names(regulon)[i], names(regulon[[i]]$tfmode)[j]]
			pos = pos + 1
		}
	}
	
# Save the regulon files
	colnames(final_mat) = c("Hub", "Target", "MI", "MoA", "likelihood", "pvalue")
	final_6col_file = paste(resultDir, resultName, ".txt", sep="")
	write.table(final_mat, file = final_6col_file, row.names = FALSE, col.names = TRUE, sep="\t", quote = FALSE)
	final_reg_file = paste(resultDir, resultName, ".rda", sep="")
	save(regulon, file=final_reg_file)
	
# Cleanup
	unlink(paste0(dirname(temp_dir), "/", basename(temp_dir)), recursive=TRUE)
	
}