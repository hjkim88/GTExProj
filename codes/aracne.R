
# Source the Utils.R file, if not already loaded.
if (Sys.info()["nodename"] %in% c("C2B2AFPL8", "C2B2AFPL9", "C2B2AFPL10", "C2B2AFPD11", "C2B2AFPL7")){
  UTILS_FILE = "C:/Users/floratos/workspace/R/labProjects/Common/Utils.R"
  ARACNE_GO_FILE = "C:/Users/floratos/workspace/R/labProjects/GTEx/R/aracne_go.R"
} else if (Sys.info()["nodename"] == "MDDSB-L004"){
  UTILS_FILE = "C:/Users/af2202/workspace/R/labProjects/Common/Utils.R"
  ARACNE_GO_FILE = "C:/Users/af2202/workspace/R/labProjects/GTEx/R/aracne_go.R"
} else if (Sys.info()["nodename"] == "C2B2AFPD9" || Sys.info()["nodename"] == "DESKTOP-F24420B"){
  UTILS_FILE = "C:/Research/CUMC/GTExProj/codes/Utils.R"
  ARACNE_GO_FILE = "C:/Research/CUMC/GTExProj/codes/aracne_go.R"
} else if (Sys.info()["nodename"] == "afdev5.c2b2.columbia.edu"){
  UTILS_FILE = "/ifs/home/c2b2/af_lab/floratos/cvs/labProjects/Common/Utils.R"
  ARACNE_GO_FILE = "/ifs/home/c2b2/af_lab/floratos/cvs/labProjects/GTEx/R/aracne_go.R"
} else if  (Sys.info()["user"]=='jb3401' || Sys.info()["user"]=='joshuabroyde'){
  source ('~/jb3401/scripts/labProjects/Common/Broyde/rFunctions.R')
  UTILS_FILE='~/jb3401/scripts/Utils/FloratosUtils.R'
} else if (Sys.info()["nodename"] == "C2B2AFPL6"){
  UTILS_FILE = "C:/Users/floratos/eclipse-workspace/R/labProjects/Common/Utils.R"
  ARACNE_GO_FILE = "C:/Users/floratos/eclipse-workspace/R/labProjects/GTEx/R/aracne_go.R"
} else if (any(grepl(Sys.info()["nodename"], c("C2B2AFPD6", "C2B2AFPD10", "C2B2AFPD12"), ignore.case = TRUE))) {
  if (Sys.info()["sysname"] == "Windows") {
    UTILS_FILE = "C:/repository/floratosLabCVS/labProjects/Common/Utils.R"
    ARACNE_GO_FILE = "C:/repository/floratosLabCVS/labProjects/GTEx/R/aracne_go.R"
  } else if (Sys.info()["sysname"] == "Linux"){
    UTILS_FILE = "/mnt/c/repository/floratosLabCVS/labProjects/Common/Utils.R"
    ARACNE_GO_FILE = "/mnt/c/repository/floratosLabCVS/labProjects/GTEx/R/aracne_go.R"
  }
} else {
  stop("No UTILS_FILE defined")
}

source(UTILS_FILE, chdir = TRUE)
source(ARACNE_GO_FILE)

# Map of TCGA abbreviations to full tumor name
tcga_abbr = c("Bladder Urothelial Carcinoma", "Breast invasive carcinoma", "Colon adenocarcinoma", 
		"Glioblastoma multiforme", "Head and Neck squamous cell carcinoma", "Kidney renal clear cell carcinoma", 
		"Kidney renal papillary cell carcinoma", "Acute Myeloid Leukemia", "Brain Lower Grade Glioma", 
		"Liver hepatocellular carcinoma", "Lung adenocarcinoma", "Lung squamous cell carcinoma", 
		"Ovarian serous cystadenocarcinoma", "Prostate adenocarcinoma", "Rectum adenocarcinoma", "Sarcoma", 
		"Skin Cutaneous Melanoma", "Stomach adenocarcinoma", "Thyroid carcinoma", "Uterine Corpus Endometrial Carcinoma")
names(tcga_abbr) = c("blca", "brca", "coad", "gbm", "hnsc", "kirc", "kirp", "laml", "lgg", "lihc", "luad", 
		"lusc", "ov", "prad", "read", "sarc", "skcm", "stad", "thca", "ucec")

# Names of ARACNe files on disk.
# fileNames = c("Adipose-Subcutaneous_vst/Adipose-Subcutaneous_vst_6cols.txt", 
#               "Adipose-Visceral-Omentum-_vst/Adipose-Visceral-Omentum-_vst_6cols.txt", 
#               "Adrenal_Gland_vst/Adrenal_Gland_vst_6cols.txt", "Artery-Aorta_vst/Artery-Aorta_vst_6cols.txt", 
#               "Artery-Coronary_vst/Artery-Coronary_vst_6cols.txt", "Artery-Tibial_vst/Artery-Tibial_vst_6cols.txt", 
#               "Brain-Caudate-basalganglia-_vst/Brain-Caudate-basalganglia-_vst_6cols.txt", 
#               "Brain-CerebellarHemisphere_vst/Brain-CerebellarHemisphere_vst_6cols.txt", 
#               "Brain-Cerebellum_vst/Brain-Cerebellum_vst_6cols.txt", "Brain-Cortex_vst/Brain-Cortex_vst_6cols.txt", 
#               "Brain-FrontalCortex-BA9-_vst/Brain-FrontalCortex-BA9-_vst_6cols.txt", 
#               "Brain-Nucleusaccumbens-basalganglia-_vst/Brain-Nucleusaccumbens-basalganglia-_vst_6cols.txt", 
#               "Breast_vst/Breast_vst_6cols.txt", 
#               "Cells-EBV-transformedlymphocytes_vst/Cells-EBV-transformedlymphocytes_vst_6cols.txt", 
#               "Cells-Transformedfibroblasts_vst/Cells-Transformedfibroblasts_vst_6cols.txt", 
#               "Colon-Sigmoid_vst/Colon-Sigmoid_vst_6cols.txt", "Colon-Transverse_vst/Colon-Transverse_vst_6cols.txt", 
#               "Esophagus-GastroesophagealJunction_vst/Esophagus-GastroesophagealJunction_vst_6cols.txt", 
#               "Esophagus-Mucosa_vst/Esophagus-Mucosa_vst_6cols.txt", 
#               "Esophagus-Muscularis_vst/Esophagus-Muscularis_vst_6cols.txt", 
#               "Heart-AtrialAppendage_vst/Heart-AtrialAppendage_vst_6cols.txt", 
#               "Heart-LeftVentricle_vst/Heart-LeftVentricle_vst_6cols.txt", "Liver_vst/Liver_vst_6cols.txt", 
#               "Lung_vst/Lung_vst_6cols.txt", "Muscle_vst/Muscle_vst_6cols.txt", "Nerve_vst/Nerve_vst_6cols.txt", 
#               "Pancreas_vst/Pancreas_vst_6cols.txt", "Pituitary_vst/Pituitary_vst_6cols.txt", 
#               "Prostate_vst/Prostate_vst_6cols.txt", 
#               "Skin-NotSunExposed-Suprapubic-_vst/Skin-NotSunExposed-Suprapubic-_vst_6cols.txt", 
#               "Skin-SunExposed-Lowerleg-_vst/Skin-SunExposed-Lowerleg-_vst_6cols.txt", "Spleen_vst/Spleen_vst_6cols.txt", 
#               "Stomach_vst/Stomach_vst_6cols.txt", "Testis_vst/Testis_vst_6cols.txt", "Thyroid_vst/Thyroid_vst_6cols.txt", 
#               "WholeBlood_vst/WholeBlood_vst_6cols.txt")

# Names of ARACNe files on disk.
fileNames = c("tcga_blca/tcga_blca_6cols2.txt","tcga_brca/tcga_brca_6cols2.txt",
		"tcga_cesc/tcga_cesc_6cols2.txt","tcga_coad/tcga_coad_6cols2.txt",
		"tcga_esca/tcga_esca_6cols2.txt","tcga_gbm/tcga_gbm_6cols2.txt",
		"tcga_hnsc/tcga_hnsc_6cols2.txt","tcga_kirc/tcga_kirc_6cols2.txt",
		"tcga_kirp/tcga_kirp_6cols2.txt","tcga_laml/tcga_laml_6cols2.txt",
		"tcga_lgg/tcga_lgg_6cols2.txt","tcga_lihc/tcga_lihc_6cols2.txt",
		"tcga_luad/tcga_luad_6cols2.txt","tcga_lusc/tcga_lusc_6cols2.txt",
		"tcga_ov/tcga_ov_6cols2.txt","tcga_paad/tcga_paad_6cols2.txt",
		"tcga_pcpg/tcga_pcpg_6cols2.txt","tcga_prad/tcga_prad_6cols2.txt",
		"tcga_read/tcga_read_6cols2.txt","tcga_sarc/tcga_sarc_6cols2.txt",
		"tcga_skcm/tcga_skcm_6cols2.txt","tcga_stad/tcga_stad_6cols2.txt",
		"tcga_tgct/tcga_tgct_6cols2.txt","tcga_thca/tcga_thca_6cols2.txt",
		"tcga_thym/tcga_thym_6cols2.txt","tcga_ucec/tcga_ucec_6cols2.txt")

# Names of DIGGIT files of disk
# fileNamesDG = c("associations-blca.dat",  "associations-laml.dat",  "associations-read.dat",
# 		"associations-brca.dat",  "associations-lgg.dat",   "associations-sarc.dat",
# 		"associations-coad.dat",  "associations-lihc.dat",  "associations-skcm.dat",
# 		"associations-gbm.dat",   "associations-luad.dat",  "associations-stad.dat",
# 		"associations-hnsc.dat",  "associations-lusc.dat",  "associations-thca.dat",
# 		"associations-kirc.dat",  "associations-ov.dat",    "associations-ucec.dat",
# 		"associations-kirp.dat",  "associations-prad.dat")

# Names of VIPER files of disk
# fileNamesVP = c("viperValues-blca.dat", "viperValues-brca.dat", "viperValues-coad.dat", 
# 		"viperValues-gbm.dat", "viperValues-hnsc.dat", "viperValues-kirc.dat", 
# 		"viperValues-kirp.dat", "viperValues-laml.dat", "viperValues-lgg.dat", 
# 		"viperValues-lihc.dat", "viperValues-luad.dat", "viperValues-lusc.dat", 
# 		"viperValues-ov.dat", "viperValues-prad.dat", "viperValues-read.dat", 
# 		"viperValues-sarc.dat", "viperValues-skcm.dat", "viperValues-stad.dat", 
# 		"viperValues-thca.dat", "viperValues-ucec.dat")

# fileNamesVP = c("viperValues-brca.dat", "viperValues-gbm.dat")

# Variable names to be used for storing the data for each network.
# if (!exists("varNames"))
# 	varNames = c("AdiposeSub", "AdiposeVis", "Adrenal", "ArteryAor", "ArteryCor", "ArteryTib", "BrainCau", 
# 	             "BrainCerHem", "BrainCer", "BrainCor", "BrainFroCor", "BrainNuc", "Breast", "CellsEBV", 
# 	             "CellsTra", "ColonSig", "ColonTra", "EsophagusGasJun", "EsophagusMuc", "EsophagusMus", 
# 	             "HeartAtrApp", "HeartLefVen", "Liver", "Lung", "Muscle", "Nerve", "Pancreas", "Pituitary", 
# 	             "Prostate", "SkinNotSun", "SkinSunExp", "Spleen", "Stomach", "Testis", "Thyroid", "WholeBlood")

# Variable names to be used for storing the data for each network.
if (!exists("varNames"))
	varNames = c("tcga_blca","tcga_brca","tcga_cesc","tcga_coad","tcga_esca","tcga_gbm","tcga_hnsc",
			"tcga_kirc","tcga_kirp","tcga_laml","tcga_lgg","tcga_lihc","tcga_luad","tcga_lusc",
			"tcga_ov","tcga_paad","tcga_pcpg","tcga_prad","tcga_read","tcga_sarc","tcga_skcm",
			"tcga_stad","tcga_tgct","tcga_thca","tcga_thym","tcga_ucec")

# Variable names to be used for storing the DIGGIT data for each cancer.
# if (!exists("varNamesDG"))
# 	varNamesDG = c("blcaDG",  "lamlDG",  "readDG", "brcaDG",  "lggDG",   "sarcDG", "coadDG",  "lihcDG",  "skcmDG",
# 		"gbmDG",   "luadDG",  "stadDG", "hnscDG",  "luscDG",  "thcaDG", "kircDG",  "ovDG",    "ucecDG",
# 		"kirpDG",  "pradDG")

# Variable names to be used for storing the VIPER data for each cancer.
# if (!exists("varNamesVP"))
# 	varNamesVP = c("blcaVP", "brcaVP", "coadVP", "gbmVP", "hnscVP", "kircVP", "kirpVP", "lamlVP", 
# 		"lggVP", "lihcVP", "luadVP", "luscVP", "ovVP", "pradVP", "readVP", "sarcVP", "skcmVP", 
# 		"stadVP", "thcaVP", "ucecVP")

# This is a number larger that the max gene ID for a TF among all networks. This 
# should be pre-computed and reset if a different set of networks is ever used.
max_geneId = 120000000

# Total number of genes used to construct Aracne networks
total_geneNum <- 31674

# Set this to control logging of messages
LOGGING_ON = TRUE


# **********************************************************************
# Reads the ARACNe network files generated for the CPTAC project and
# organizes the data into data structures convenient for exploration. 
# The intention is to run this commnand once, to generate the data 
# structures and then store them to an .rda file for subsequent use.
#
# One variable is created for each ARACNe network.Variables are names 
# with the standard TCGA tumor type acronyms (blca, brca, etc.). Each 
# variable A is a list with two elements:
# *  A[[1]] is a Nx3 matrix, where N is the number of TF is in the original ARACNe network, i.e., the 
#    number of unique genes B such that (B, C) is a network edge. Each row R in that matrix corresponds
#    to a single TF B and contains the following data:
#    ** A[[1]][R,1] is the Gene ID of gene B.
#    ** A[[1]][R,2] is the index corresponding to B in the list A[[2]] -- this will be explained below.
#    ** A[[1]][R,3] is the size of the regulon of B in the network, i.e., the number of its predicted targets.
# *  A[[2]] is a list with N items, again each corresponding to a TF in the ARACNe network. The R-th element
#    of the list contains information related to the TF B for which A[[1]][,2] = R. The R-th list element 
#    A[[2]][[R]] is a matrix with dimensions M x 5 where M is the size if the regulon of the TF B. Each 
#    Row S in that matrix corresponds to a predicted target C of B and contains the following data:
#    ** A[[2]][[R]][S, 1] is the gene ID of gene C.
#    ** A[[2]][[R]][S, 2] is the mutual information of the edge (B,C), as computed by ARACNe.          
#    ** A[[2]][[R]][S, 3] is the Mechanism of Action of the edge (B,C), as computed by ARACNe.
#    ** A[[2]][[R]][S, 4] is the Likelihood of the edge (B,C), as computed by ARACNe.
#    ** A[[2]][[R]][S, 3] is the P-value of the edge (B,C), as computed by ARACNe\n
# According to the above scheme, the targets of the TF stored at A[[1]][i,] can be found at the list
# element A[[2]][[A[[1]][i,2]]]. A few more things to note about variable A:
# * The matrix A[[1]] is ordered according to TF regulon size, with the largest TFs first. 
# * The targets contained in the matrix A[[2]][[R]] are listed in decreasing order of significance, based on
#   P-value, with the lowest P-value listed first. 
#
#
# The (basic) logic of the code below is as follows (for each processed ARACNe file):
# * Create a vector with a unique list of TF gene ids, say of length N.
# * Create a matrix mat with dimensions Nx3. Store the following:
#   ** mat[,1] = the N unique TF ids.
#   ** mat[,2] = indices 1:N
#   ** mat[,3] = initialize all values to 0.
# * Go over the ARACNe network and for each row corresponding to TF A, increase by one the contents of mat[index(A),3].
# * Create a list L of length N where L[[i]] contains an empty numerical vector of lenght mat[i,3].
# * Reset mat[,3] to all 0.
# * Go over the ARACNe network again. For each row corresponding to TF A, describing and edge (A,B).
#   ** increase by one mat[index(A),3].
#   ** set L[[index(A)]][mat[index(A),3]] = B.
# **********************************************************************
readAracneNetworks <- function(rootDir = "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/CPTAC/ARACNe/"){
	
	rootDir = paste(rootDir, "/", sep = "") 	# to be safe...
	
	# We use this vector to track duplicate interactions. About 10% of interactions
	# in each network involve pairs of TFs A, B and appear both as (A, B) and (B, A).
	# The vector below will help us identify those. They will still be both stored
	# but the second one will be stored as (B, -A). By making the gene ID of A negative, 
	# we can still maintain the info we need (abs(A) will still give us the ID of A) 
	# and at the same time have a way to identify duplicates easily.
	trackDuplicates = vector(mode = "integer", length = max_geneId)
	
	# Read ARACNe networks one by one
	for(index in 1:length(fileNames)){
		logLines(paste("\nProcessing network -> ", fileNames[index]))
		network = as.matrix(read.table(paste(rootDir, fileNames[index], sep=""), header = TRUE))
		
		# Reconstitute the TF->TF interactions (y, x) where both x and y are TFs, (x, y)
		# is reported in the network by (y, x) is not
		missing = getMissingInteractions(network)
		colnames(missing) = colnames(network)
		logLines(paste("\t# of missing TF -> TF interactions =", length(missing[,1])))
		
		# Create the final network by adding the missing interactions.
		network = rbind(network, missing)
		
		# Order first by p-value. This will guarantee that later on, when building the list of 
		# targets for a TF, the targets will be added in order of significance, starting with
		# the one that has the smallest p-value.
		network = network[order(network[,6]),]   
		
		# Order then by TF gene id. This step is important for our code that computes
		# duplicates. Specifically, it guarantees that if the gene Id of A is smaller
		# than the gene Id of B then then interaction (A,B) will be encountered before
		# the interaction (B,A)
		network = network[order(network[,1]),] 
		
		# Initialize the duplicates tracking array. Ne do not need to do this the
		# first time around, as the values are already set correctly.
		if (index != 1)
			for (i in 1:max_geneId)
				trackDuplicates[i] = 0
		
		uniqueTFs = sort(unique(network[,1]))
		N = length(uniqueTFs)
		mat = matrix(nrow = N, ncol = 3)
		for (i in 1:N){
			mat[i,1] = uniqueTFs[i]
			mat[i,2] = i
			mat[i,3] = 0
			trackDuplicates[uniqueTFs[i]] = i
		}
		
		# count the number of targets for each TF
		j = 1
		mat[1,3] = 1
		M = length(network[,1])
		for (i in 2:M){
			if(network[i,1] != network[i-1,1])
				j = j+1
			mat[j, 3] =  mat[j, 3]+1
		}
		
		
		listOfTargets = list()
		# for eact TF, create a matrix Tx5 where T is the number of targets of the TF.
		# each matrix row will represent an edge in the network and will be named by the
		# gene ID of the target gene
		for (i in 1:N){
			listOfTargets[[i]] = matrix(nrow = mat[i,3], ncol=5)
			colnames(listOfTargets[[i]]) = c("Target", "MI", "MoA", "Likelihood", "Pvalue")
		}
		
		# populate the matrices with the actual targets
		j = 1
		k = 1
		# Fill in the first entry
		for (r in 1:5)
			listOfTargets[[1]][1,r] = network[1,r+1]
		# Fill in the rest of the entries
		for (i in 2:M){
			if(network[i,1] != network[i-1,1]){
				rownames(listOfTargets[[j]]) = abs(listOfTargets[[j]][,1])
				j = j+1
				k = 0;
			}
			k = k+1
			for (r in 1:5)
				listOfTargets[[j]][k, r] = network[i, r+1]
			
			# Check if this a duplicate interaction
			if (network[i,2] <= max_geneId && trackDuplicates[network[i,2]] != 0  &&
					network[i,1] %in% listOfTargets[[trackDuplicates[network[i,2]]]][,1]){
				listOfTargets[[j]][k, 1] = -network[i,2]
			}
		}
		rownames(listOfTargets[[j]]) = abs(listOfTargets[[j]][,1])
		
		results = list()
		# Use the gene IDs of the TFs to name the rows of the 'mat' matrix and 
		# the elements of 'listOfTargets' list, to facilitate searching 
		rownames(mat) = mat[,1]
		names(listOfTargets) = mat[,1]
		
		# Order mat so that that TFs with the largest regulons are listed first.
		mat = mat[order(mat[,3], decreasing=TRUE),]
		results[[1]] = mat
		results[[2]] = listOfTargets
		assign(varNames[index], results, envir = globalenv())
	}
	
	assign("netSizes", interactionCounts(), envir = globalenv())
	assign("pairWise", interactionCounts("pairwise"), envir = globalenv())
}


# *****************************************************************************
# Print out a description of what the variables generated from processing the
# ARACNe network files.
# *****************************************************************************
README = function(){
	writeLines("This workspace contains a number of variable related to the ARACNe networks generated")
	writeLines("run the command ls() to see the full listing.")
	
	writeLines("---> The number of ARACNe networks: K")
	
	writeLines("---> One variable for each ARACNe network was generated. These")
	writeLines("variables are named with tissue names")
	writeLines("Specifically, each such variable A is a list with two elements:")
	writeLines("* A[[1]] is a Nx3 matrix, where N is the number of TF is in the original ARACNe network, i.e., the ")
	writeLines("   number of unique genes B such that (B, C) is a network edge. Each row R in that matrix corresponds")
	writeLines("   to a single TF B and contains the following data:")
	writeLines("   ** A[[1]][R,1] is the Gene ID of gene B.")
	writeLines("   ** A[[1]][R,2] is the index corresponding to B in the list A[[2]] -- this will be explained below.")
	writeLines("   ** A[[1]][R,3] is the size of the regulon of B in the network, i.e., the number of its predicted targets.")
	writeLines("* A[[2]] is a list with N items, again each corresponding to a TF in the ARACNe network. The R-th element")
	writeLines("   of the list contains information related to the TF B for which A[[1]][,2] = R. The R-th list element ")
	writeLines("   A[[2]][[R]] is a matrix with dimensions M x 5 where M is the size if the regulon of the TF B. Each ")
	writeLines("   Row S in that matrix corresponds to a predicted target C of B and contains the following data:")
	writeLines("   ** A[[2]][[R]][S, 1] is the gene ID of gene C. If C is TF and the interaction (C,B) has been seen")
	writeLines("      before the value of this cell is -(gene ID of gene C).")
	writeLines("   ** A[[2]][[R]][S, 2] is the mutual information of the edge (B,C), as computed by ARACNe.          ")
	writeLines("   ** A[[2]][[R]][S, 3] is the Mechanism of Action of the edge (B,C), as computed by ARACNe.")
	writeLines("   ** A[[2]][[R]][S, 4] is the Likelihood of the edge (B,C), as computed by ARACNe.")
	writeLines("   ** A[[2]][[R]][S, 3] is the P-value of the edge (B,C), as computed by ARACNe")
	writeLines("According to the above scheme, the targets of the TF stored at A[[1]][i,] can be found at the list")
	writeLines("element A[[2]][[A[[1]][i,2]]]. A few more things to note about variable A:")
	writeLines("* The matrix A[[1]] is ordered according to TF regulon size, with the largest TFs first. ")
	writeLines("* The targets contained in the matrix A[[2]][[R]] are listed in decreasing order of significance, based on")
	writeLines("   P-value, with the lowest P-value listed first.\n")
	
	writeLines("--->  A matrix called 'pairWise' of dimensions K x K where, for i <> j, pairWise[i, j] is the")
	writeLines("number of interactions shared between the i-th and the j-th network. The network names are also listed")
	writeLines("as column and row headings. When i = h then pairWise[i, i] is the number of interactions in the -th")
	writeLines("network. In both cases, only unique interactions are counted. I.e., if both (A,B) and (B,A) are ")
	writeLines("reported for a given network then only one of these is taken into consideration in the interaction")
	writeLines("counts.\n")
	
	writeLines("--->  A vector called 'netSizes' with K entries one for each networks. The i-th entry contains the")
	writeLines("number of unique interactions in network varNames[i], i.e., duplicates have been discarded. The")
	writeLines("entries in the vector are named	after 'varNames'.\n")
	
	writeLines("--->  A variable called 'tfPairEnrich' describing the TF-specific regulon conservation across pairs of")
	writeLines("interactomes. It comprises the following elements:")
	writeLines("* tfPairNet[[1]] is a K x K symmetric matrix where the entry [i,i] is zero and the entry [i, j] is")
	writeLines("  the index within tfPairEnrich[[2]] where the results of the pairwise comparison between the i-th")
	writeLines("  and the j-th interactome are stored (i and j are indices within the varNames variable).")
	writeLines("* tfPairEnrich[[2]] is a list with choose(K, 2) elements, each corresponding to a pairwise")
	writeLines("  interactome comparison. Each element is a N x 5 matrix M where N is the # number of TFs that appear")
	writeLines("  as hubs both in the i-th and the j-th interactome. Each row corresponds to a TF A and contains the")
	writeLines("  following 5 columns:")
	writeLines("  ** The gene id of the TF A.")
	writeLines("  ** The value round(log(P), 0) where P is the p-value of the the Fisher exact test to assess the size")
	writeLines("     of the intersection of the regulons of A in the i-th and the j-th interactome.")
	writeLines("  ** The size of the intersection of the 2 regulons.")
	writeLines("  ** The size of the regulon of A in the i-th interactome.")
	writeLines("  ** The size of the regulon of A in the j-th interactome.")
	writeLines("  Within the matrix M, rows are ordered in increasing value of the 2nd column, i.e., the most enriched")
	writeLines("  TFs (those with the smallest log(P)) are listed first. When computing Fisher's exact test we use the")
	writeLines("  following 2 x 2 contingency matrix:")
	writeLines("               X	Y")
	writeLines("               Z	W")
	writeLines("where:")
	writeLines("* X is the size of the intersection of the 2 regulons.")
	writeLines("* Y is the size of regulon(A) in the i-th interactome minus X.")
	writeLines("* Z is the size of regulon(A) in the j-th interactome minus X.")
	writeLines("* W is the total number of interactions in the j-th interactome minus the size of regulon(A) in the")
	writeLines("  j-th interactome.")
	writeLines("Notice that Y, Z, W are normalized based on the size of the intersections of the 2 interactomes.\n")
	
	writeLines("--->  A variable called 'tfPairProb' describing the another TF-specific regulon conservation across")
	writeLines("pairs of interactomes. It comprises the following elements:")
	writeLines("* tfPairProb[[1]] is a K x K symmetric matrix where the entry [i,i] is zero and the entry [i, j] is")
	writeLines("  the index within tfPairProb[[2]] where the results of the pairwise comparison between the i-th")
	writeLines("  and the j-th interactome are stored (i and j are indices within the varNames variable).")
	writeLines("* tfPairProb[[2]] is a list with choose(K, 2) elements, each corresponding to a pairwise")
	writeLines("  interactome comparison. Each element is a N x 5 matrix M where N is the # number of TFs that appear")
	writeLines("  as hubs both in the i-th and the j-th interactome. Each row correponds to a TF A and contains the")
	writeLines("  following 5 columns:")
	writeLines("  ** The gene id of the TF A.")
	writeLines("  ** The value log(P) where P is the probability that a regulon pair share same target genes")
	writeLines("  ** The size of the intersection of the 2 regulons.")
	writeLines("  ** The size of the regulon of A in the i-th interactome.")
	writeLines("  ** The size of the regulon of A in the j-th interactome.")
	writeLines("  Within the matrix M, rows are ordered in increasing value of the 2nd column.\n")
	
	writeLines("--->  A variable called 'tfNetEnrich' has regulon conservation info between every possible")
	writeLines("hub pair in each interactome (tissue). This is different from tfPairEnrich,")
	writeLines("since this occurred in each interactome of different hubs while the tfPairEnrich")
	writeLines("took place in different interactomes of regulons of one same hub gene.")
	writeLines("It has a list with length of \"varNames\", which means the length of")
	writeLines("existing interactomes. And each element in the list, there is a matrix")
	writeLines("that contains the regulon conservation info.")
	writeLines("Detailed descriptions of the result are:")
	writeLines("* tfNetEnrich is a list object such that:")
	writeLines("- length(tfNetEnrich) = length(varNames).")
	writeLines("- names(tfNetEnrich) = varNames")
	writeLines("* Let net_name be a value from varNames and let X = get(net_name)")
	writeLines("Then tfNetEnrich[[net]] is a a symmetric NxN matrix M, where N is the number of hubs")
	writeLines("in the interactome \"net\" and M[hub1, hub2] = M[hub2, hub1] = -log10(p),")
	writeLines("where p is the p-value of the FET for the intersection of the regulons of hub1 and hub2.")
	writeLines("Also, rownames(M) = colnames(M) = rownames(get(net)[[1]]) and M[i,i] = Inf.")
	writeLines("")
	writeLines("When computing Fisher's exact test we use the following 2 x 2 contingency matrix:")
	writeLines("")
	writeLines("                 Regulon2  No-Regulon2")
	writeLines("                -----------------------")
	writeLines("   Regulon1    |     X	          Y")
	writeLines("   No-Regulon1 |     Z	          W")
	writeLines("")
	writeLines("   where:")
	writeLines("    - X is the size of the intersection of the 2 regulons : the number of shared genes")
	writeLines("    - Y is the size of the first regulon minus X.")
	writeLines("    - Z is the size of the other (second) regulon minus X.")
	writeLines("    - W is the total number of genes in the interactome minus (X plus Y plus Z).")
	writeLines("The p-value of the FET will be one-sided p-value with alternative = \"greater\" option,")
	writeLines("which means it is a test of the odds ratio being bigger than 1.\n")
}



# *****************************************************************************
# Count number of unique intreractions.
#
# ARGUMENTS:
# * selector:	What is returned depends on the value of the variable "selector":
# 	-- "single":	Returns a named vector with 20 entries, one for each networks. The
# 		i-th entry contains the number of unique interactions in network varNames[i]
#		i.e., duplicates have been discarded. The entris in the vector are named
#		after 'varNames'
# 	-- "pairwise":		Returns a N x N matrix with one row and one column for each
#		interactome. The [i,j] entry contains the number of interactions shared
#		by the i-th and j-th interactome. The matrix has column and row names,
#		named after the interactomes
# 	-- "allUnique":	Returns a vector that contains 2 values
#		- The number of unique interactions across all networks.
#		- The number of all interactions across all networks.
# * filterFUN: Ff nor NULL, 'filterFUN' is expected to be one of the functions 
# 		is.tf(), is.cotf(), or is.sign(). This function will then be used to filter
# 		interactions involving only TF, co-TF, or signaling hubs respectively. If
# 		filterFUN == NULL, then all interactions are counted.
# * thresh:	only interactions with p-values <= 'thresh' will be used in the 
#		calculations. If 'thresh' is NULL then all interactions will be used.
#
# RETURN VALUE:
# Depends on the value of the argument selector, as describe above.
# *****************************************************************************
interactionCounts <- function(selector = "single", filterFUN = NULL, thresh = NULL){
	
	# Check that argument is OK
	if (selector != "single" & selector != "pairwise" & selector != "allUnique"){
		print("The value of the argument should be either 'single' or 'pairwise' or 'allUnique'")
		stop()
	}
	
	finalCounts = vector(mode="integer", length=length(varNames))
	fcInd = 1
	
	if (is.null(filterFUN))
		filterFUN <- function(e) {return(TRUE)}
	if (is.null(thresh))
		thresh = 1
	
	if (selector == "single"){
		# Go over each interactome and count non-duplicate interactions, i.e., 
		# interactions (A, B) where B is a positive number (A is assummed to
		# be the TF we are currently traversing).
		for (i in 1:length(varNames)){
			count = 0;
			nList = get(varNames[i])[[2]]
			for (j in 1:length(nList)){
				if (filterFUN(names(nList)[j])){
					x = nList[[j]][,"Target"]
					y = nList[[j]][x>0, "Pvalue"]
					count = count + length(y[y <= thresh])
				}
			}
			finalCounts[fcInd] = count
			fcInd = fcInd + 1
		}
		names(finalCounts) = varNames
		return(finalCounts)
	}
	else if (selector == "pairwise"){
		# Construct the results array and name rows and columns
		pairCounts = matrix(nrow = length(varNames), ncol = length(varNames))
		row.names(pairCounts) = varNames
		colnames(pairCounts) = varNames
		
		for (i in 1:length(varNames)){
			logLines(paste("net1 = ", i))
			net1 = get(varNames[i])[[1]]
			net1 = net1[order(net1[,1]),]
			net1_det = get(varNames[i])[[2]]
			
			for (j in (i+1):length(varNames)){
				if (j > length(varNames))
					break
				logLines(paste("\tnet2 = ", j))
				count = 0
				net2 = get(varNames[j])[[1]]
				net2 = net2[order(net2[,1]),]
				net2_det = get(varNames[j])[[2]]
				n_1 = 1
				n_2 = 1
				# Step down the list of TFs in the two networks that are being compared
				# and look for identical TFs. When such a pair is found we count how many 
				# interactions they have in common, excluding duplicate interactions
				while (n_1 <= length(net1[,1]) & n_2 <= length(net2[,1])){
					if (net1[n_1,1] < net2[n_2, 1])
						n_1 = n_1 + 1
					else if (net1[n_1,1] > net2[n_2, 1])
						n_2 = n_2 + 1
					else{
						if (filterFUN(net1[n_1,1])){
							ind1 = net1[n_1, 2]
							ind2 = net2[n_2, 2]
							t1 = net1_det[[ind1]][, "Pvalue"]
							t2 = net2_det[[ind2]][,"Pvalue"]
							t1 = t1 <= thresh
							t2 = t2 <= thresh
							x = net1_det[[ind1]][t1, "Target"]
							y = net2_det[[ind2]][t2, "Target"]
							count = count + length(intersect(x[x>0], y[y>0]))
						}
						n_1 = n_1 + 1
						n_2 = n_2 + 1
					}
				}
				logLines(paste("count = ", count))
				pairCounts[i,j] = count
				pairCounts[j,i] = count
			}
		}
		# Finally, count the number of unique interactions in each network and 
		# populate the diagonal of the results matrix. 
		x = interactionCounts(thresh = thresh)
		for (i in 1:length(x))
			pairCounts[i, i] = x[i]
		return(pairCounts)
	}
	else if (selector == "allUnique"){
		tfCounts = vector(mode = "integer", length = max_geneId)
		
		# Count the number of interactions for each TF, excluding duplicates
		for (ind in 1:length(varNames)){
			net = get(varNames[ind])[[1]]
			netInts = get(varNames[ind])[[2]]
			L = length(net[,1])
			for (i in 1:L){
				if (filterFUN(net[i, 1])){
					tf_id = net[i,1]
					tf_ind = net[i,2]
					t1 = netInts[[tf_ind]][, "Pvalue"]
					t1 = t1 <= thresh
					tf_ints = netInts[[tf_ind]][t1,"Target"]
					count = length(tf_ints[tf_ints>0])
					tfCounts[tf_id] = tfCounts[tf_id] + count
				}
			}
		}
		logLines("Done counting the number of targets for all TFs across all nets")
		
		# This is the total number of interactions summed up across all networks
		totCount = sum(tfCounts)
		
		# Create the record-keeping apparatus:
		# - One vector of length N for each TF, when N are the total number of
		#   targets of TF across all networks.
		# - One index for each TF, to indicate the next available position in 
		#   the TFs vector
		N = length(tfCounts[tfCounts > 0])
		listOfTargets = list()
		indVector = vector(mode="integer", length = N)
		indVector[1:N] = 1
		j = 1
		for (i in 1:max_geneId){
			if (tfCounts[i] > 0){
				listOfTargets[[j]] = vector(mode="integer", length = tfCounts[i])
				tfCounts[i] = j
				j = j + 1
			}
		}
		logLines("Done allocating record-keeping space.")
		
		# Now go over all networks again and  populate the TF target vectors
		for (ind in 1:length(varNames)){
			net = get(varNames[ind])[[1]]
			netInts = get(varNames[ind])[[2]]
			L = length(net[,1])
			for (i in 1:L){
				if (filterFUN(net[i, 1])){
					tf_id = net[i,1]
					tf_ind = net[i,2]
					t1 = netInts[[tf_ind]][, "Pvalue"]
					t1 = t1 <= thresh
					tf_ints = netInts[[tf_ind]][t1,"Target"]
					tf_ints = tf_ints[tf_ints>0]
					# logLines(paste(length(tf_ints)))
					if (length(tf_ints) > 0){
						beginIndex = indVector[tfCounts[tf_id]]
						endIndex = beginIndex + length(tf_ints) - 1
						listOfTargets[[tfCounts[tf_id]]][beginIndex:endIndex] = tf_ints
						indVector[tfCounts[tf_id]] = endIndex + 1
					}
				}
			}
		}
		
		# Finally, compute the total number of unique interactions
		totUniqCount = 0
		for (i in 1:N)
			totUniqCount = totUniqCount + length(unique(listOfTargets[[i]]))
		
		# Return a vector of size two, containing the total number of unique
		# interactions and the total number of interactions
		return(c(totUniqCount, totCount))
		# return(listOfTargets)
	}
}




# *****************************************************************************
# Count number of genes in interactomes.
#
# ARGUMENTS:
# * selector:	What is returned depends on the value of the variable "selector":
# 	-- "single":	Returns a named vector with length(varNames) entries, one for 
# 		each networks. The i-th entry contains the number of genes in network 
#		varNames[i], counting both hub genes and targets.
# 	-- "pairwise":		Returns a N x N matrix with one row and one column for each
#		interactome (i.e., N = length(varNames)). The [i,j] entry contains the number
#		 of genes shared by the i-th and j-tj interactome. The matrix has column 
#		and row names named after the interactomes
#
# RETURN VALUE:
# Depends on the value of the argument selector, as describe above.
# *****************************************************************************
geneCounts <- function(selector = "single"){
	intCounts = sapply(varNames, getInteractomeGenes)
	if (selector == "single"){
		return(intCounts)
	}
	else if (selector == "pairwise"){
		res = matrix(0, nrow = length(varNames), ncol = length(varNames))
		rownames(res) = colnames(res) = varNames
		L = sapply(varNames, function(x){return(getInteractomeGenes(x, count=FALSE))})
		for (i in 1:(length(varNames) -1)){
			for (j in (i+1):length(varNames)){
				res[i, j] = res[j ,i] = length(intersect(L[[i]], L[[j]]))
			}
		}
		for (i in 1:length(varNames))
			res[i,i] = intCounts[i]
		return(res)
	}
	else
		stop("Function geneCounts(): No valid value for argument 'selector' provided.")
}



getMissingInteractions <- function(network){
	
	# Compile a list of unique TFs in the network
	uniqueTFs = unique(network[,1])
	logLines(paste("\tlength(uniqueTFs) = ", length(uniqueTFs)))
	
	# Each row corresponds to an EntrezId. The x-th row will be non-zero only
	# if the gene with EntrezID = x is a TF
	tfIndices = vector(mode = "integer", length = max_geneId)
	for (i in 1:length(uniqueTFs))
		tfIndices[uniqueTFs[i]] = i
	
	# Let x, y be TFs with entrez IDs gid_x, gid_y. Let i_x = tfIndices[gid_x]
	# and i_y = tfIndices[gid_y]. Then pairwiseTFs[i, j] will be set to 1 only
	# if (x, y) is an interaction in the network
	pairwiseTFs = matrix(nrow = length(uniqueTFs), ncol = length(uniqueTFs))
	pairwiseTFs[ , ] = 0
	
	# Go over all interaction (x, y) and mark those where both x and y are TFs    
	countInts = 0
	L = length(network[, 1])
	for (i in 1:L){
		i_x = tfIndices[network[i, 1]]
		i_y = tfIndices[network[i, 2]]
		if (i_x != 0 && i_y != 0){
			countInts = countInts + 1
			pairwiseTFs[i_x, i_y] = 1
		} 
	}
	
	# Again, go over all interactions and now record all interactions
	# (y, x) where x and y are TFs, (x, y) is in the network, and (y, x) is not.
	# For (y, x) use the numerical values (MI, likelihood, P-value, etc.)
	# of the (x, y) interaction
	mat = matrix(nrow = 2*countInts, ncol = 6)
	mat[ , ] = 0
	countInts = 0
	for (i in 1:L){
		i_x = tfIndices[network[i, 1]]
		i_y = tfIndices[network[i, 2]]
		if (i_x != 0 && i_y != 0 && pairwiseTFs[i_y, i_x] != 1){
			countInts = countInts + 1
			mat[countInts, 3:6] = network[i, 3:6]
			mat[countInts, 1] = network[i, 2]
			mat[countInts, 2] = network[i, 1]
		}
	}  
	
	# Finally, return only those interactions (y, x) where (x, y) is in the network
	# but (y, x) is not
	return(mat[1:countInts,])
}



panNetwork <- function(){
	
	tfCounts = vector(mode = "integer", length = max_geneId)
	
	# Count the number of interactions for each TF, excluding duplicates
	for (ind in 1:length(varNames)){
		net = get(varNames[ind])[[1]]
		netInts = get(varNames[ind])[[2]]
		L = length(net[,1])
		for (i in 1:L){
			tf_id = net[i,1]
			tf_ind = net[i,2]
			tf_ints = netInts[[tf_ind]][,1]
			count = length(tf_ints[tf_ints>0])
			tfCounts[tf_id] = tfCounts[tf_id] + count
		}
	}
	logLines("Done counting the number of targets for all TFs across all nets")
	
	# Create the results object, i.e., a list with one entry for each TF. Each entry
	# comprises a matrix of dimensions M x 6 for each TF, when M is the number of
	# targets of the TF across all networks. The first column will contain the gene id
	# of the TF target, columns 2-5 will contain the interaction metrics from the 
	# ARACNe files, and column 6 will contain a number from 1-20 providing the index
	# in vector varNames[] of the network where the interaction comes from.
	listOfTargets = list()
	
	# Summary object - matrix of size N x 2, where N is the number of TFs, eash entry 
	# corresponding to one TF. As before:
	# - mat[i,1] = gene id of i-th TF.
	# - mat[i,2] = index within the list listOfTargets corresponding to the i-th TF.
	N = length(tfCounts[tfCounts > 0])
	mat = matrix(nrow = N, ncol = 2)
	
	# Helper vector, one index for each TF, to indicate the next available position in 
	# the TFs vector
	indVector = vector(mode="integer", length = N)
	indVector[1:N] = 1
	
	j = 1
	for (i in 1:max_geneId){
		if (tfCounts[i] > 0){
			mat[j, 1] = i
			mat[j, 2] = j
			listOfTargets[[j]] = matrix(nrow = tfCounts[i], ncol = 6)
			tfCounts[i] = j
			j = j + 1
		}
	}
	logLines("Done allocating record-keeping space.")
	
	# Now go over all networks again and  populate the TF target vectors
	for (ind in 1:length(varNames)){
		net = get(varNames[ind])[[1]]
		netInts = get(varNames[ind])[[2]]
		L = length(net[,1])
		logLines(paste("L = ", L))
		for (i in 1:L){
			tf_id = net[i,1]
			tf_ind = net[i,2]
			tf_ints = netInts[[tf_ind]]
			tf_ints = tf_ints[tf_ints[,1]>0,, drop=FALSE]
			if (length(tf_ints[,1]) > 0){
				beginIndex = indVector[tfCounts[tf_id]]
				endIndex = beginIndex + length(tf_ints[,1]) - 1
				listOfTargets[[tfCounts[tf_id]]][beginIndex:endIndex, 1:5] = tf_ints[, 1:5]
				listOfTargets[[tfCounts[tf_id]]][beginIndex:endIndex, 6] = ind
				indVector[tfCounts[tf_id]] = endIndex + 1
			}
		}
	}
	
	
	# Order target lists first by p-value and then by ID. This will guarantee that target gene IDs 
	# appear in increasing order and that, for all interactions involving the same target,
	# interactions appear in order of significance, i.e., smallest p-value first.
	for (i in 1:N){
		listOfTargets[[i]] = listOfTargets[[i]][order(listOfTargets[[i]][,5]),]
		listOfTargets[[i]] = listOfTargets[[i]][order(listOfTargets[[i]][,1]),]
	}
	
	results = list()
	results[[1]] = mat
	results[[2]] = listOfTargets
	return(results)
}



# *****************************************************************************
# Compute gene-specific regulon conservation across pairs on interactomes
#
# ARGUMENTS
# * norm_method:	This arguments determines how the normalization of regulon
#		sizes is implemented for the the FET calculation. There are two choices
#		(see below for details):
#		- "interactions":	normalize based on the number of shared interactions
#				between the two interactomes being compared.
#		- "genes":			normalize based on the number of shared genes
#				between the two interactomes being compared.
#
# RETURN VALUE
# Return a 'results' object that is a list of 2 members:
# results[[1]]:	A N x N symmetric matrix (N = length(varNames)) where the entry 
#	[i,i] is zero and the entry [i, j] is the index within results[[2]] where 
#	the results of the pairwise comparison between the i-th and the j-th 
#	interactome are stored (i and j are indices within the varNames
#	variable).
# results[[2]]:	A list with choose(length(varNames), 2) elements, each 
#	corresponding to a pairwise interactome comparison. Each element is a N x 5 
# 	matrix M where N is the number of TFs that appear as hubs both in the 
#	i-th and the j-th interactome. Each row corresponds to a TF A and
#	contains the following 5 columns:
#	* The gene id of the TF A.
#	* The value round(log(P), 0) where P is the p-value of the the Fisher
#	  exact test to assess the size of the intersection of the regulons of
#	  A in the i-th and the j-th interactome.
#	* The size of the intersection of the 2 regulons.
#	* The size of the regulon of A in the i-th interactome.
#	* The size of the regulon of A in the j-th interactome.
#	
#	Within the matrix M, rows are ordered in increasing value of the 2nd
#	column, i.e., the most enriched TFs (those with the smallest log(P))
#	are listed first.
#
# When computing Fisher's exact test we use the following 2 x 2 contingency
# matrix:
#			X	Y
#			Z	W
# where:
# * X is the size of the interesection of the 2 regulons.
# * Y is the size of regulon(A) on the i-th interactome minus X.
# * Z is the size of regulon(A) on the j-th interactome minus X.
# * W is the total number of interactions in the j-th interactome -
#   the size of regulon(A) on the j-th interactome.
#
# Notice that Y, Z, W are normalized based on the size of the intersections of
# the 2 interactomes. This is a little tricky. In the typical FET setup, there
# is a single universal set of objects which is the source of elements for the 
# 2 sets whose intersection we are assessing. Here instead, the 2 sets are 
# regulons A,B (for the same TF) coming from different interactomes N1 and N2 
# which are not identical. Rather, they have both shared and non-shared 
# interactions (or genes). We cannot pretend that the universal set is the union 
# of the interactions (or genes) in N1 and N2 because that union is not accessible to both
# regulons: interactions (or genes) for regulon A come exclusively from N1 and 
# for regulon B from N2. To work around this problem we focus on an 
# artificial universal set, namely the intersection I of N1 and N2. We then 
# consider how many of A's (and B's) interactions we would expect by chance to
# find in I, based on the size of N1, N2, and I. Finally, we use FET to assess
# the actual observed intersection of A and B against the expected sizes of A 
# and B when projected on I. This is not perfect but, I think, is a reasonable
# and fair way to assess the intersection taking into account the sizes of A, B, 
# N1, N2, and I.
#
# It should also be noted that FET is run with alternative = "greater". I.e., we 
# are interested in interactomes that share more interactions than expected by chance 
# but *NOT* in interactomes that share less (the latter will get a p-value close to 1  
# and, thus, will be mostly disregarded in downstream analyses).
# *****************************************************************************
tfPairEnrichment <- function(norm_method = "interactions"){
	
	# The results objects
	mat = matrix(0, nrow=length(varNames), ncol=length(varNames))
	enrichmentScores = list()
	row.names(mat) = varNames
	colnames(mat) = varNames
	listLocation = 1
	
	for (i in 1:(length(varNames) - 1)){
		logLines(paste("net1 = ", i))
		net1 = get(varNames[i])[[1]]
		net1 = net1[order(net1[,1]),]
		net1_det = get(varNames[i])[[2]]
		
		for (j in (i+1):length(varNames)){
			logLines(paste("\tnet2 = ", j))
			count = 0
			net2 = get(varNames[j])[[1]]
			net2 = net2[order(net2[,1]),]
			net2_det = get(varNames[j])[[2]]         
			
			# Size of the intersection of the 2 interactomes, needed to normalize
			# the regulon sizes for FET calculation below.
			if (norm_method == "interactions"){
				sizeInt = pairWise[i, j]
				size1 = netSizes[i]
				size2 = netSizes[j]
			}
			else if (norm_method == "genes"){
				if (!exists("pairwiseGeneCounts"))
					pairwiseGeneCounts = geneCounts("pairwise")
				sizeInt = pairwiseGeneCounts[i, j]
				size1 = pairwiseGeneCounts[i, i]
				size2 = pairwiseGeneCounts[j, j]
			}
			
			enrichmentMat = matrix(nrow = length(intersect(net1[,1], net2[,1])), ncol = 5)
			
			n_1 = 1
			n_2 = 1
			pos = 1
			# Step down the list of TFs in the two networks that are being compared
			# and look for identical TFs. When such a pair is found we count how many 
			# interactions they have in common, excluding duplicate interactions
			while (n_1 <= length(net1[,1]) & n_2 <= length(net2[,1])){
				if (net1[n_1,1] < net2[n_2, 1])
					n_1 = n_1 + 1
				else if (net1[n_1,1] > net2[n_2, 1])
					n_2 = n_2 + 1
				else{
					ind1 = net1[n_1, 2]
					ind2 = net2[n_2, 2]
					reg1 = abs(net1_det[[ind1]][,1])
					reg2 = abs(net2_det[[ind2]][,1])
					
					#Normalized regulon sizes
					reg1_norm = round(length(reg1) * sizeInt/size1, 0)
					reg2_norm = round(length(reg2) * sizeInt/size2, 0)
					common_targets = intersect(reg1, reg2)
					common = length(common_targets)
					r1_minus_r2 = max(reg1_norm, common) - common
					r2_minus_r1 = max(reg2_norm, common) - common
					remainder = sizeInt - (common + r1_minus_r2 + r2_minus_r1)
					fet = fisher.test(rbind(c(common, r1_minus_r2), c(r2_minus_r1, remainder)), alternative = "greater")
					
					enrichmentMat[pos, 1] = net1[n_1,1]
					enrichmentMat[pos, 2] = round(log10(fet[[1]]))
					enrichmentMat[pos, 3] = common
					enrichmentMat[pos, 4] = length(reg1)
					enrichmentMat[pos, 5] = length(reg2)
					
					pos = pos + 1
					n_1 = n_1 + 1
					n_2 = n_2 + 1
					
				}
			}
			
			enrichmentMat = enrichmentMat[order(enrichmentMat[,2]),]
			rownames(enrichmentMat) = enrichmentMat[,1]
			colnames(enrichmentMat) = c("hub_gene", "log10_pval", "common", "regulon1", "regulon2")
			
			# Update the results objects
			enrichmentScores[[listLocation]] = enrichmentMat
			mat[i, j] = listLocation
			mat[j, i] = listLocation
			listLocation = listLocation + 1
		}
	}
	results = list()
	results[[1]] = mat
	results[[2]] = enrichmentScores
	return(results)
}   


### A function to calculate probability that a regulon pair share same target genes
#   if the number of total genes = l
#      the number of target genes in regulon A = m
#      the number of target genes in regulon B = n
#      the number of shared target genes between A and B = k
#
#   Prob = (C(l, k) x C (l-k, m-k) x C(l-m, n-k)) / (C(l, m) x C(l, n))
#
#   * C(a, b) = choose(a,b) = P(a, b) / P(b, b)
#   * C() = Combination, P() = Permutation
###
tfPairProbability <- function() {
	
	# The results objects
	mat = matrix(0, nrow=length(varNames), ncol=length(varNames))
	probScores = list()
	row.names(mat) = varNames
	colnames(mat) = varNames
	listLocation = 1
	
	### this library is needed to use chooseZ
	if(!require(gmp)) {
		install.packages("gmp")
		library(gmp)
	}
	
	### A function to calculate the probability if the l, m, n, k are given
	calProb <- function(l, m, n, k) {
		return(log10(as.double((chooseZ(l, k) * chooseZ(l-k, m-k) * chooseZ(l-m, n-k)) / (chooseZ(l, m) * chooseZ(l, n)))))
	}
	
	for (i in 1:(length(varNames) - 1)){
		logLines(paste("net1 = ", i))
		net1 = get(varNames[i])[[1]]
		net1 = net1[order(net1[,1]),]
		net1_det = get(varNames[i])[[2]]
		size1 = netSizes[i]
		
		for (j in (i+1):length(varNames)){
			logLines(paste("\tnet2 = ", j))
			count = 0
			net2 = get(varNames[j])[[1]]
			net2 = net2[order(net2[,1]),]
			net2_det = get(varNames[j])[[2]]         
			size2 = netSizes[j]
			
			# Size of the intersection of the 2 interactomes, needed to normalize
			# the regulon sizes for FET calculation below.
			sizeInt = pairWise[i, j]
			
			probMat = matrix(nrow = length(intersect(net1[,1], net2[,1])), ncol = 5)
			
			n_1 = 1
			n_2 = 1
			pos = 1
			# Step down the list of TFs in the two networks that are being compared
			# and look for identical TFs. When such a pair is found we count how many 
			# interactions they have in common, excluding duplicated interactions
			while (n_1 <= length(net1[,1]) & n_2 <= length(net2[,1])){
				if (net1[n_1,1] < net2[n_2, 1])
					n_1 = n_1 + 1
				else if (net1[n_1,1] > net2[n_2, 1])
					n_2 = n_2 + 1
				else{
					ind1 = net1[n_1, 2]
					ind2 = net2[n_2, 2]
					reg1 = abs(net1_det[[ind1]][,1])
					reg2 = abs(net2_det[[ind2]][,1])
					
					#Normalized regulon sizes
					reg1_norm = round(length(reg1) * sizeInt/size1, 0)
					reg2_norm = round(length(reg2) * sizeInt/size2, 0)
					common_targets = intersect(reg1, reg2)
					common = length(common_targets)
					
					### calculate the probs
					prob = calProb(total_geneNum, length(reg1), length(reg2), common)
					
					probMat[pos, 1] = net1[n_1,1]
					probMat[pos, 2] = prob
					probMat[pos, 3] = common
					probMat[pos, 4] = length(reg1)
					probMat[pos, 5] = length(reg2)
					
					pos = pos + 1
					n_1 = n_1 + 1
					n_2 = n_2 + 1
					
				}
			}
			
			probMat = probMat[order(probMat[,2]),]
			rownames(probMat) = probMat[,1]
			colnames(probMat) = c("hub_gene", "log10_pval", "common", "regulon1", "regulon2")
			
			# Update the results objects
			probScores[[listLocation]] = probMat
			mat[i, j] = listLocation
			mat[j, i] = listLocation
			listLocation = listLocation + 1
		}
	}
	
	results = list()
	results[[1]] = mat
	results[[2]] = probScores
	
	return(results)
}


# *****************************************************************************
#
# The result of the function has regulon conservation info between every possible
# hub pair in each interactome (tissue). This is different from tfPairEnrichment(),
# since this occurs in each interactome of different hubs while the tfPairEnrichment()
# takes place in different interactomes of regulons of one same hub gene.
#
# This function returns a list with length of "varNames", which means the length of
# existing interactomes. And each element in the list, there is a matrix
# that contains the regulon conservation info.
#
# Detailed descriptions of the result are:
# * tfNetEnrich is a list object such that:
#   - length(tfNetEnrich) = length(varNames).
#   - names(tfNetEnrich) = varNames
# * Let net_name be a value from varNames (e.g., net_name = "Liver") and let  X = get(net_name) . 
#	Then tfNetEnrich[[net]] is a a symmetric NxN matrix M, where N is the number of hubs
# in the interactome "net" and M[hub1, hub2] = M[hub2, hub1] = -log10(p),
# where p is the p-value of the FET for the intersection of the regulons of hub1 and hub2.
# Also, rownames(M) = colnames(M) = rownames(get(net)[[1]]) and M[i,i] = Inf.
#
# * When computing Fisher's exact test we use the following 2 x 2 contingency matrix:
#
#                 Regulon2  No-Regulon2
#                -----------------------
#   Regulon1    |     X	          Y
#   No-Regulon1 |     Z	          W
#
#   where:
#     - X is the size of the interesection of the 2 regulons : the number of shared genes
#     - Y is the size of the first regulon minus X.
#     - Z is the size of the other (second) regulon minus X.
#     - W is the total number of genes in the interactome minus (X plus Y plus Z).
#
#   The p-value of the FET will be one-sided p-value with alternative = "greater" option,
#   which means it is a test of the odds ratio being bigger than 1. I.e., we are interested
#	in interactomes that share more interactions than expected by chance but *NOT* in 
#	interactomes that share less (the latter will get a p-value close to 1 and, thus, 
#	will be mostly disregarded in downstream analyses).
#
#   * The input parameter "partialFileDir" is set to NULL
#     If it is not null, a result of every tissue will be saved in RDA partially
#     Plus, it is before filtering the hubs (p-value cutoff does not apply)
#     This can help recover the process later, since it is a long-running process
#     If a system failure happens, we can recover the results from those RDA files
#
# *****************************************************************************
tfNetEnrichment <- function(partialFileDir=NULL) {
	
	### calculate the total number of genes for each interactome
	### this will be used in FET calculation later
	totalGeneCounts <- geneCounts("single")
	
	### create the results list (a named list),  with one entry per interactome
	tfNetEnrich <- vector("list", length = length(varNames))
	names(tfNetEnrich) <- varNames
	
	### the process runs on every interactome in varNames
	for(net in varNames) {
		
		### write a log for each process run
		logLines(paste("Processing interactome ->", net))
		
		### get the interactome object
		X <- get(net)
		
		### create an empty matrix for current interactome
		tfNetEnrich[[net]] <- matrix(NA, nrow = nrow(X[[1]]), ncol = nrow(X[[1]]))
		rownames(tfNetEnrich[[net]]) <- rownames(X[[1]])
		colnames(tfNetEnrich[[net]]) <- rownames(X[[1]])
		
		### the process runs on every hub in a given interactome
		for(j in 1:(nrow(tfNetEnrich[[net]])-1)) {
			
			### write a log for each process run
			logLines(paste("\tProcessing hub # ->", j))
			
			### get the given hub
			hub1 = rownames(tfNetEnrich[[net]])[j]
			
			### make an empty vector to store p-value of FET
			### this will be forwarded to the matrix later
			### since direct save to a matrix takes ridiculously slow
			buf_vec <- NULL
			
			### Compare hub1 with all other hubs in the interactome
			for(k in (j+1):ncol(tfNetEnrich[[net]])) {
				
				### get a second hub for comparison
				hub2 <- colnames(tfNetEnrich[[net]])[k]
				
				### set the second column - p-value of FET between two regulons (conservativeness)
				reg1 <- rownames(X[[2]][[hub1]])
				reg2 <- rownames(X[[2]][[hub2]])
				common <- length(intersect(reg1, reg2))
				r1_minus_r2 <- length(reg1) - common
				r2_minus_r1 <- length(reg2) - common
				remainder <- totalGeneCounts[net] - (common + r1_minus_r2 + r2_minus_r1)
				buf_vec <- c(buf_vec, -log10(fisher.test(rbind(c(common, r1_minus_r2), c(r2_minus_r1, remainder)), alternative = "greater")$p.value))
				
			}
			
			### store the p-value of FET in the matrix
			### vectorized approach
			tfNetEnrich[[net]][hub1,(j+1):ncol(tfNetEnrich[[net]])] <- buf_vec
			
		}
		
		### tfNetEnrich[[net]][i,i] = Inf
		diag(tfNetEnrich[[net]]) <- Inf
		
		### there are some already-calculated stuffs since FET(X,Y) == FET(Y,X)
		### just copy them into the appropriate places
		tfNetEnrich[[net]][lower.tri(tfNetEnrich[[net]])] <- t(tfNetEnrich[[net]])[lower.tri(tfNetEnrich[[net]])]
		
		### sort the matrix
		tfNetEnrich[[net]] <- tfNetEnrich[[net]][order(as.integer(rownames(tfNetEnrich[[net]]))), order(as.integer(colnames(tfNetEnrich[[net]])))]
		
		### if the partialFileDir is set, save the result (for the current tissue) in RDA
		if(!is.null(partialFileDir)) {
			assign(paste0("tfNetEnrich_", net), tfNetEnrich[[net]], envir = globalenv())
			save(list = c(paste0("tfNetEnrich_", net)), file = paste0(partialFileDir, "tfNetEnrich_", net, ".rda"))
			rm(list = c(paste0("tfNetEnrich_", net)))
		}
		
		### garbage collection
		gc()
		
	}
	
	### return the result
	return(tfNetEnrich)
	
}


# *****************************************************************************
# Return the regulon of a specific TF gene, possibly filtered to include only
# interactions that clear MI and P-value thresholds.
#
# Arguments are:
# * geneHub: The gene id (either as integer or as string) or the gene symbol of 
#		the query TF gene. 
# * net:	The networks to use. This is either a subvector of varNames (i.e., a
#		character vector of interactome variable names), or a *single* ARACNe 
#		network variable (i.e., the actulal list object), or NULL. In the 
#		first two cases the function returns the regulon of the TF from the
#		specified network(s). Otherwise regulons for all networks are returned.
# * mode:	Specifies how much data to return:
#	* "all": Return the full M x 5 matrix corresponding to the query gene A
#		from net[[2]][[index(A)]].
#	* "ids": Return a vector V with the ids of the target genes in the regulon.
# * mi:		MI threshold. Returns only interactions whose MI is at least as 
#		large as this threshold.
# * pval:	P-value threshold. Returns only interactions whose P-value is at 
#		most as large as this threshold.
#
# Returns either a vector V or a matrix with 5 columns as described above. 
# If 'net' == NULL or 'net' is a character vector of length > 1, return a list
# with length(net) elements (or length(varNames), if 'net' == NULL), one for each 
# ARACNe network, each element being the vector V or matrix M described above 
# for the corresponding network.
# *****************************************************************************
getRegulon <- function(geneHub, net = NULL, mode="ids", mi = 0, pval = 1){
  geneHub = strtoi(as.entrezId(geneHub))
  if (is.na(geneHub))
    return(NULL)
  
  if (is.null(net))
    net = varNames
  if(is.character(net) && length(net) > 1)
    return(sapply(net, function(n){return(getRegulon(geneHub, n, mode, mi, pval))}))
	
	# Make sure the name of the network exists
	if (is.character(net))
		if(net %in% varNames)
			net = get(net)
		else
			return(NULL)
	
	x = net[[1]]
	ind = x[x[,1] == geneHub, 2]
	if (length(ind) > 0){
		mat = net[[2]][[ind]]
		sel = (mat[, "MI"] >= mi) & (mat[, "Pvalue"] <= pval)
		if (length(sel)== 0 | sum(sel)==0)
			return(NULL)
		mat = mat[sel, , drop=FALSE]
		if (mode == "all"){
			mat[,1] = abs(mat[,1])
			return(mat)
		}
		if (mode == "ids"){
			#res = abs(mat[,1])
			#names(res) = NULL
			#return(res)
			return(abs(mat[,1]))
		}
	}
	else
		return(NULL)
}


# *****************************************************************************
# Return genes that are the intersection or union of the regulons of multiple
# bub genes
#
# ARGUMENTS:
# * gene_hubs:	The query hub genes. This can be a vector of integers
#		representing entrez ids or a character vector containing gene 
#		symbols or entrez ids in character string form.
# * nets:		The networks to search. This variable is a character vector 
#		containing entries from varNames specifying interactome variables;
#		or it can be NULL. In the latter case, regulons from all networks in
#		varNames are used.
# * mode:	A string assuming one of the following two values:
#	* "shared":	indicates that the method should return genes that are found 
#			in the intersection of all regulons of the query hubs in "gene_hubs".
#	* "all": indicates that the method should return genes that are found in the 
#			union of all regulons of the query hubs in "gene_hubs".
#
# RETURN VALUE
# Returns a list with one entry for each network in "nets" (the list is named
# with the corresponding entries from varNames). The list entry for 
# interactome A contains a character vector with the entrez ids of the 
# genes that are either in the intersection (if mode == "shared") or the 
# union (if mode == "all") of the regulons of the hub genes in gene_hubs in 
# interactome A.
# *****************************************************************************

getManyRegulons <- function(gene_hubs, nets = NULL, mode = c("shared", "all")){
	gene_hubs = as.entrezId(gene_hubs)
	
	if (is.null(nets))
		nets = varNames
	if (!is.character(nets))
		return(NULL)
	
	res = lapply(nets, function(net){
				net = get(net)
				
				regs = sapply(gene_hubs, function(hub){
							return(getRegulon(hub, net = net))
						})
				if (mode[1] == "shared")
					fun = intersect
				else
					fun = union
				targets = Reduce(fun, regs)
				return(targets)
			})
	names(res) = nets
	return(res)
}


# *****************************************************************************
# Run pathway enrichment analysis for the targets of a hub gene that appear 
# multiple times across different interactomes.
#
# Given a query hub gene H, this method tabulates the interactions in the regulon
# of H across a user-specified set of interactomes, counting only interactions that
# clear thresholds for MI and p-value. For each target gene T belonging to such an
# interaction we compute its support, i.e., the number of interactomes contain the 
# interaction (H, T). The targets with the highest support are then analyzed using
# a GO/pathway enrichment analysis.
#
# ARGUMENTS:
# * gene:	The query hub gene, specified as a gene symbol (character string) or as 
#		an Entrez Id (either integer or character string).
# * nets:	The networks to be used for computing support for targets of the hub gene. 
#		Any subsect of the vector varNames can be passed as a value. By default, the 
#		entire vector varNames is used. As a shorthand, the value of this argument can 
#		also be set to either "TCGA" or "GTEx", in which case only the TCGA or GTEx
#		netwrorks will be used.
# * top:	The number of top-support target genes to use for the enrichment analysis.
# * wd:		Directory where to save the image files generated bt the function.
# * heatmap:	Specifies if the method makeGraphs(which = "regulon_conservation") should
#		be called also, to generate a heatmap plotting the regulon conservation of the
#		query hub gene across all pairs of TCGA and GTEx interactomes.
# * save:	if TRUE then the heatmap above is stored in an image file.
# * mi, pval:	MI and P-value thresholds to use when counting interactions. Only interactioms
#		with MI >= mi and P-value <= pval will be considered when tabulating gene target
#		support.
# * pa_thr:	the adjusted p-value used by the enrichment methods to decide if a GO term or
#		a pathway should be deemed enriched and be reported.
# * type:	either "CP" or "TB". Specifies which pathway enrichment method to use, either
#		pathwayAnalysis_CP() or pathwayAnalysis_TB().
#
# RETURN VALUE
# Writes image files in the directory 'wd' and returns the object generated by
# the call to methods pathwayAnalysis_CP() or pathwayAnalysis_TB().
# *****************************************************************************
regulonAnnotation <- function(gene, nets = varNames, top = 100, heatMap = TRUE, wd = ".", 
		pa_thr = 0.05, save=FALSE, type=c("CP", "TB"), mi = 0.2, pval = 0.0001, p_num = 40){
	cd = getwd()
	setwd(wd)
	if (heatMap)
		makeGraphs(which = "regulon_conservation", params = List(gene, "heatmap", 100), save = save)
	if (nets[1] == "TCGA"){
		nets = varNames[grep("tcga", varNames)]
	}else if (nets[1] == "GTEx") 
		nets = varNames[-grep("tcga", varNames)]
	t = tabulateRegulonInteractions(gene, mi = mi, pval = pval, nets=nets)
	t = t[1:min(length(t), top)]
	if (type[1] == "CP")
		pa = pathwayAnalysis_CP(names(t), org="human", database="GO", pv_threshold = pa_thr, displayNum = p_num)
	else
		pa = pathwayAnalysis_TB(entrezIDtoSymbol(names(t)), FDR_threshold = pa_thr, displayNum = p_num)
	setwd(cd)
	return(pa)
	
}

# *****************************************************************************
# Find all regulators R of a query target gene G, across one or multiple 
# interactomes, cosnidering only interactions (R, G) that clear MI and 
# p-value thresholds.
#
# Arguments are:
# * geneID: The query gene, represented as a gene symbol (character string) or
#		an Entrez ID (either a character string or an integer).
# * mode:	Specifies how much data to return:
#	* "ids": Return a vector V containing only the ids of the regulators.
#	* "all": Return a matrix M with one row per regulator and 5 columns. M[i,1]
#		is the gene ID of the i-th regulator. M[i, 2:5] are the interaction 
#		metrics for the interaction (regulator_i, geneID) from the ARARCNe network.
# * net:	The network(s) to use. This is either vector that is a subset of 
#		varNames or the actual ARACNe network variable. If 'net' == NULL then
#		regulators for each tumor are returned.
# * mi:		Only interactions with MI >= than this threshold are considered.
# * pval:	Only interactions with p-value <= than this threshold are considered.
# 
# Returns either a vector V or a matrix M with 5 columns as described above. 
# If 'net' == NULL return a list with one entry one for each interactome in varNames, 
# each entry being the vector V or matrix M described above for the corresponding 
# interactome. if 'net' is a subset of varNames returns a list as above, but only
# with entries corresponding to the networks in 'net'.	
# *****************************************************************************
getRegulators <- function(geneID, mode = "ids", net = NULL, mi = 0, pval = 1){
	
	geneID =  as.entrezId(geneID)
	if (is.na(geneID))
		return(NULL)
	
	if (is.null(net))
		return(sapply(varNames, function(e){return(getRegulators(geneID, mode, e, mi, pval))}))
	if (is.character(net) && (length(net) > 1))
		return(sapply(net, function(e){return(getRegulators(geneID, mode, e, mi, pval))}))
	if (is.character(net))
		net = get(net)
	
	res = matrix(nrow = length(net[[2]]), ncol = 5)
	colnames(res) = c("Hub", colnames(net[[2]][[1]])[2:5])
	k = 0
	for (i in 1:length(net[[2]])){
		if ((geneID %in% rownames(net[[2]][[i]])) && (net[[2]][[i]][geneID, "MI"] >= mi)
				&& (net[[2]][[i]][geneID, "Pvalue"] <= pval)){
			k = k+1
			res[k,1] = as.integer(names(net[[2]])[i])
			res[k, 2:5] = net[[2]][[i]][geneID, 2:5]
		}
	}
	if (k == 0)
		return(NULL)
	if (mode == "ids")
		return(res[1:k,1])
	res = res[1:k, , drop = FALSE]
	rownames(res) = res[, 1]
	return(res)
}


# *****************************************************************************
# Given a TF, assess the intersection of the regulons of that TF across  two
# interactomes. The statistics on enrichment used here are copied over from 
# the tfPairEnrich object.
#
# Arguments are:
# * gene:	A gene identifier if the form of an integer/string gene ID or a 
#		string gene symbol.
# * net1_name, net2_name:	The interactomes to use, represented as the
#		relevant strings from 'varNames'.
# 
# Returns a list L comprising 2 members:
# * L[[1]]: this is a vector containing the gene IDs of the genes that 
#		comprise the intersection of the 2 regulons.
# * L[[2]]: this is a vector containing the following values:
#	** geneID
#	** the log10 p-value of the FET test.
#	** the number of genes in the intersection of the 2 regulons (i.g., |L[[1]]|
#	** the number of genes in the first regulon
#	** the number of genes in the second regulon
# *****************************************************************************

getConservedRegulon <- function(gene, net1_name, net2_name){
	gene = as.entrezId(gene)
	if (is.na(gene))
		return(NULL)
	
	if (!is.character(net1_name) | !is.character(net2_name))
		stop("Interactome names must be strings")
	if (is.character(net1_name) & is.character(net2_name)){
		# Check that the network names exist
		if(!all(net1_name %in% varNames, net2_name %in% varNames))
			return(NULL)
		
		net1 = get(net1_name)
		net2 = get(net2_name)
	}
	
	x = net1[[1]]
	ind = x[x[,1] == gene, 2]
	if (length(ind) == 0)
		return(NULL)
	reg1 = abs(net1[[2]][[ind]][,1])
	
	x = net2[[1]]
	ind = x[x[,1] == gene, 2]
	if (length(ind) == 0)
		return(NULL)
	reg2 = abs(net2[[2]][[ind]][,1])
	common_targets = intersect(reg1, reg2)
	
	x = tfPairEnrich[[2]][[tfPairEnrich[[1]][net1_name, net2_name]]]
	results = list()
	results[[1]] = common_targets
	results[[2]] = x[x[,1] == gene	,]
	
	return(results)
}

# *****************************************************************************
# Assess how conserved the regulon of a TF is across all pairs of interactomes
# specified in the argument "nets"
#
# ARGUMENTS:
# * gene:	The gene id of the query TF, in the form of either a gene symbol
#			or Entrez ID (the latter either as string or integer)
# * nets:	A vector of length N contains the names (as character strings) 
# 			of the interactomes to use. The default is all the interactomes.
# * mode:	It assumes either of two values, "vector" or "matrix". It specifies 
#			if the results will be returned as a vector or as a	matrix.
#
# RETURN VALUE
# ---->If mode == "vector":
# Returns a vector with choose(N,2) entries, one  for each interactome pair.
# Interactomes are compared in the order in which they are listed in 'varNames',
# i.e., the first interactome is compared with the remaining N-1, followed by
# comparing the second interactome with the remaining N-2, etc. If the i-th 
# vector entry corresponds to the comparison of the interactomes A and B, then
# the value in this entry is the log(p) value of the Fisher Exact Test (FET) 
# that assesses the intersection of the regulon of "gene" in A with the regulon of
# geneID in B, i.e., the enrichment of the first regulon in genes from the 
# second regulon. For more details about the FET see the comments of the 
# function tfPairEnrichment.
#
# ----> If mode == "matrix":
# Returns the same FEt values but organized in an N x N matrix M with row names
# and column names being the values in "nets" and where M[i,j] is the FET for the
# intersection of regulons of "gene" in nets[i] and nets[j]. 
# ***************************************************************************** 
regulonConservationAcrossNets <- function(gene, nets = varNames, mode="vector"){
	gene = strtoi(as.entrezId(gene))
	
	N = length(nets)
	
	if (mode == "vector"){
		pos = vector("integer", choose(N,2))
		k = 1
		for (i in 1:(N-1))
			for (j in (i+1):N){
				pos[k] = tfPairEnrich[[1]][nets[i], nets[j]]
				k = k + 1
			}
		
		res = sapply(tfPairEnrich[[2]][pos], function(e){
					if (length(which(e[,1] == gene)) == 0)
						return (0)
					else
						return(e[which(e[,1] == gene),2])
				})
		return(res)
	}
	
	if (mode == "matrix"){
		res = matrix(-Inf, N, N)
		rownames(res) = colnames(res) = nets
		for (i in 1:(N-1))
			for (j in (i+1):N){
				fetTable = tfPairEnrich[[2]][[tfPairEnrich[[1]][nets[i], nets[j]]]]
				if (length(which(fetTable[,1] == gene)) == 0)
					res[i,j] = res[j,i] = 0
				else
					res[i,j] = res[j,i] = fetTable[which(fetTable[,1] == gene),2]
			}
		return(res)
	}
}

# *****************************************************************************
# For a query hub gene, use the tfPairEnrich object to extract the FET 
# log(pvalues) for all pairwise N = choose(M, 2) interactome regulon comparisons 
# (where M is the total number of interactomes), optionally replacing -Inf values. 
# Also, if named == TRUE names the 
# vector entries as "net1%net2"
#
# ARGUMENTS:
# * gene:	the gene id of the query hub gene, in the form of either a gene symbol
# 		or Entrez ID (the latter either as string or integer).
# * replaceInf:		If TRUE, log(pvalues) equal to -Inf are replaced with an the smallest
# 		non-infinite log(pvalue), minus 10.
#
# RETURN VALUE
# A vector of length N, containing the FET log(pvalues). Vector entries are named
# with strings of the form "A%B", where A and B	are the variable names of the two 
# interactomes being compared (e.g., "ArteryAor%Lung").
# *****************************************************************************
getRegulonFETs <- function(gene, replaceInf = FALSE){
	x = regulonConservationAcrossNets(gene, mode="matrix")
	N = nrow(x)
	res = rep(0, length(choose(N, 2)))
	k = 1
	for (i in 1:(N-1))
		for (j in (i+1):N){
			res[k] = x[i, j]
			names(res)[k] = paste(rownames(x)[i], colnames(x)[j], sep="%")
			k = k+1
		}
	if (replaceInf){
		min = min(res[res != -Inf])
		res[res == -Inf] = min-10
	}
	return(res)
}

# *****************************************************************************
# For a given query hub gene goes over the regulons of that gene across a s
# specified list of interactomes and for each encountered interaction it 
# reports how many interactomes it is seen in.
#
# ARGUMENTS:
# * gene:	The query gene identifier, in the form of an integer/string gene ID  
#			or a string gene symbol.
# * nets:	A vector of length N contains the names (as character strings) 
# 			of the interactomes to use. The default is all the interactomes.
# * mi:		mi threshold to use for choosing interactions to report.
# * pval:	p-value threshold to use for choosing interactions to report.
#
# RETURN VALUE
# A vector with one entry per interaction involving the hub gene, containing
# the number of interactomes where this is interaction is present. The vector
# is ordered in decreasing value of these numbers and its entries are named 
# with the Entrez ID of the target gene in each interaction. In compiling
# this vector we only consider interactions that clear that "mi" and "pval"
# thresholds, i.e., interactions with mutual information >= "mi" and 
# p-value <= "pval"
# ***************************************************************************** 
tabulateRegulonInteractions <- function(gene, nets = varNames, mi = 0, pval = 1){
	gene = strtoi(as.entrezId(gene))
	
	res = unlist(sapply(nets, function(e){
						if (!(gene %in% get(e)[[1]][,1]))
							return(NA)
						mat = get(e)[[2]][[as.character(gene)]]
						sel = (mat[, "MI"] >= mi) & (mat[, "Pvalue"] <= pval)
						if (length(sel)==0 | sum(sel)==0)
							return(NA)
						return(abs(mat[sel, "Target"]))
					}))
	res = res[!is.na(res)]
	return(sort(table(res), decreasing = TRUE))
}


logLines <- function(message){
	
	if (LOGGING_ON){
		writeLines(message)
		flush(stdout())
	}
}


########################################################################
#                               VIPER
########################################################################

# *****************************************************************************
# Read the VIPER data files.
#
# For each tumor stores the VIPER data in a data object A (named using the format 
# 'gbmVP'). The object comprises 2 lists:
# * A[[1]]: This is a matrix of size N x 2 where the first column contains the
#	gene IDs of VIPER regulators and the second column contains the VIPER 
#	enrichment scores (NES) for these regulators. Specifically, if A[[1]][i,1] = B
#	then the regulators B has a VIPER NES of A[[1]][i,2] in the sample whose
#	sample ID is stored in A[[2]][i]. Entries in A[[1]][,] are first ordered 
#	according to the first column (i.e., regulator gene ID) and then according to
#	the second column, highest NES first.
# * A[[2]]: This is a vector of sample IDs, as described in the bullet above.
# *****************************************************************************
readVIPER <- function(rootDir = "//isilon.c2b2.columbia.edu/ifs/scratch/c2b2/ac_lab/fmg2117/projects/pancancer/cptac/citrusFiles/"){
	for(index in 1:length(fileNamesVP)){
		logLines(paste("\nProcessing file -> ", fileNamesVP[index]))
		data = as.matrix(read.table(paste(rootDir, fileNamesVP[index], sep=""), header = TRUE))
		x = cbind(as.integer(data[,1]), as.numeric(data[,3]))
		# Sort according to gene ID and then according to Viper value
		ordering = order(x[,1], x[,2], decreasing = TRUE)
		x = x[ordering, ]
		res = list()
		res[[1]] = x
		res[[2]] = data[ordering,2]
		assign(varNamesVP[index], res, envir = globalenv())
	}
}


# *****************************************************************************
# Return the VIPER profiles of one or more regulators.
#
# * geneIDs: An integer or character vector containing the gene IDs of the query
#		regulators.
# * varVP:	The name (or the actual variable) of the query tumor, in the 
#		format 'gbmVP'
#
# Returns a matrix with one row per regulator in 'geneIDs' and one column per
# sample in the query tumor. The -ith row is the VIPER profile of the regulator
# geneIDs[i]. The rows and columns of the matrix are named with, respectively,
# the contents of 'geneIDs' and the alphabetically sorted lit of sample names
# for the query tumor.
# *****************************************************************************

getViperProfiles <- function(geneIDs, varVP){
	if (is.character(varVP))
		varDG = get(varVP)
	geneIDs = as.integer(geneIDs)
	colNames = sort(unique(varVP[[2]]))
	res = matrix(nrow = length(geneIDs), ncol = length(colNames))
	rownames(res) = geneIDs
	colnames(res) = colNames
	for (i in 1:length(geneIDs)){
		ind = which(varVP[[1]][,1] == geneIDs[i])
		val = varVP[[1]][ind,2]
		names(val) = varVP[[2]][ind]
		res[i,] = val[colNames]
	}
	return(res)
}


########################################################################
#                               DIGGIT
########################################################################


# **********************************************************************
# Reads the DIGGIT files generated for the CPTAC project and organizes
# the data into data structures convenient for exploration. The intention 
# is to run this command once, to generate the data 
# structures and then store them to an .rda file for subsequent use.
#
# One variable is created for each ARACNe network.Variables are names 
# with the standard TCGA tumor type acronyms (blca, brca, etc.) followed
# by the string "DG". Each variable A is a list with two elements:
# *  A[[1]] is a Nx3 matrix, where N is the number of TFs in the DIGGIT file, i.e., the 
#    number of unique genes B found in the second DIGGIT file column. Each row R in that 
#    matrix corresponds to a single TF B and contains the following data:
#    ** A[[1]][R,1] is the Gene ID of gene B.
#    ** A[[1]][R,2] is the index corresponding to B in the list A[[2]] -- this will be explained below.
#    ** A[[1]][R,3] is the number of genomic events reported by DIGGIT to modulate B.
# *  A[[2]] is a list with N items, again each corresponding to a TF from the DIGGIT file. The R-th element
#    of the list contains information related to the TF B for which A[[1]][,2] = R. The R-th list element 
#    A[[2]][[R]] is a matrix with dimensions M x 6 where M is the number of genomic events reported by
#    DIGGIT to modulate B. Each row S in that matrix corresponds to a genomic event XXX_C modulating B, 
#    where XXX is one of SNV, AMP, DEL, and denotes the type of genetic event (single nucleotide variation, 
#    amplification, deletion) and C is the Gene ID of the gene harboring the genomic event. The S-th row
#    contains the following data:
#    ** A[[2]][[R]][S, 1] is the gene ID of gene C.
#    ** A[[2]][[R]][S, 2] is an integer representing the type of genomic event (1 = SNV, 2 = AMP, 3 = DEL).
#    ** A[[2]][[R]][S, 3] is the "areaP" p-value for the genomic event from the DIGGIT file.
#    ** A[[2]][[R]][S, 4] is the "preppiP" p-value for the genomic event from the DIGGIT file.
#    ** A[[2]][[R]][S, 5] is the "cindyP" p-value for the genomic event from the DIGGIT file.
#    ** A[[2]][[R]][S, 6] is the "combinedP" p-value for the genomic event from the DIGGIT file.
# A few more things to note about variable A:
# * The matrix A[[1]] is ordered according to number of modulating events (i.e., the number stored in the 3rd
#   column of the table), with the largest entry first. 
# * The targets contained in the matrix A[[2]][[R]] are listed in decreasing order of significance, based on
#   "combinedP", with the lowest P-value listed first. 
#
# **********************************************************************
readDIGGIT <- function(rootDir = "/ifs/scratch/c2b2/ac_lab/fmg2117/projects/pancancer/cptac/citrusFiles/"){
	
	# Read DIGGIT networks one by one
	for(index in 1:length(fileNamesDG)){
		logLines(paste("\nProcessing file -> ", fileNamesDG[index]))
		dataDG = as.matrix(read.table(paste(rootDir, fileNamesDG[index], sep=""), header = TRUE))
		
		# Order by TF gene id.
		geneID = as.integer(dataDG[,2])
		indOrder = order(geneID)
		dataDG = dataDG[indOrder,]
		geneID = geneID[indOrder]
		uniqueTFs = sort(unique(geneID))
		N = length(uniqueTFs)
		
		mat = matrix(nrow = N, ncol = 3)
		for (i in 1:N){
			mat[i,1] = uniqueTFs[i]
			mat[i,2] = i
			mat[i,3] = 0
		}
		
		# count the number of genomic events for each TF
		j = 1
		mat[1,3] = 1
		for (i in 2:length(geneID)){
			if(geneID[i] != geneID[i-1])
				j = j+1
			mat[j, 3] =  mat[j, 3]+1
		}      
		
		listOfEvents = list()
		# populate the matrices with the actual data for the genomic events 
		# found in the DIGGIT files
		j = 1
		k = 1
		# Notice, below we use an extra as.matrix() call to defend against the case
		# of one-row matrices, which R turns automatically into vectors. The as.matrix()
		# call forces those to be saved as 1-row matrices instead, as we want.
		tmpMat = as.matrix(matrix(nrow = mat[j,3], ncol=6))
		colnames(tmpMat) = c("GeneID", "EventType", "areaP", "preppiP", "cindyP", "combinedP")
		
		# Fill in the first entry
		rowData = dataDG[1,]
		parts = strsplit(rowData[1], "_")
		eventType = diggitEventType(parts[[1]][1])
		eventGene = as.integer(parts[[1]][2])
		# if the eventType is not among those expected, then there is something wrong with
		# data file - stop.
		if (eventType == 0)
			stop(paste("The following genomic event:-->",  parts[[1]][1], "<-- found in file", 
							fileNames[index], "is not of type SNV, AMP, or DEL."))
		tmpMat[1,1] = eventGene
		tmpMat[1,2] = eventType
		for (r in 3:6)
			tmpMat[1,r] = as.numeric(rowData[r])
		
		# Fill in the rest of the entries
		M = length(dataDG[,1])
		for (i in 2:M){
			if(geneID[i] != geneID[i-1]){
				# We are just finished processing a TF, order the genomic events according
				# to the combinedP value, save, and setup the processing of the next TF.
				# In the statement below, the argument drop = FALSE protects against 
				# turning one-row matrices into vectors.
				tmpMat = tmpMat[order(tmpMat[,6]) , , drop = FALSE]
				rownames(tmpMat) = tmpMat[,1]
				listOfEvents[[j]] = tmpMat
				j = j+1
				k = 0;
				# See comment avove regarding the seemingly uncessary call to as.matrix()
				tmpMat = as.matrix(matrix(nrow = mat[j,3], ncol=6))
				colnames(tmpMat) = c("GeneID", "EventType", "areaP", "preppiP", "cindyP", "combinedP")
			}
			k = k+1
			rowData = dataDG[i,]
			parts = strsplit(rowData[1], "_")
			eventType = diggitEventType(parts[[1]][1])
			eventGene = as.integer(parts[[1]][2])
			# if the eventType is not among those expected, then there is something wrong with
			# data file - stop.
			if (eventType == 0)
				stop(paste("The following genomic event:-->",  parts[[1]][1], "<-- found in file", 
								fileNames[index], "is not of type SNV, AMP, or DEL."))
			tmpMat[k,1] = eventGene
			tmpMat[k,2] = eventType
			for (r in 3:6)
				tmpMat[k,r] = as.numeric(rowData[r])
		}
		
		# Take care of the last iterant
		tmpMat = tmpMat[order(tmpMat[,6]) , , drop = FALSE]
		rownames(tmpMat) = tmpMat[,1]
		listOfEvents[[j]] = tmpMat
		
		results = list()
		# Name the list items with the gene ID of their corresponding TF
		names(listOfEvents) = mat[,1]
		# Order mat so that that TFs with the largest regulons are listed first.
		mat = mat[order(mat[,3], decreasing=TRUE),]
		results[[1]] = mat
		results[[2]] = listOfEvents
		assign(varNamesDG[index], results, envir = globalenv())
		
	}
	
	# Finally, count the total number of genomic events in each tumor type
	tmp = sapply(varNamesDG, tmp <- function(name){x = get(name); y = x[[1]][,3]; return(sum(y))})
	assign("eventCountByTumor", tmp, envir = globalenv())
	
}


# *****************************************************************************
# Print out description of variables generated by processing the DIGGIT files
# *****************************************************************************
READMEDG <-function () {
	writeLines("This workspace contains a number of variables related to the DIGGIT files generated")
	writeLines(" in the CPTAC project - run the command ls() to see the full listing. They include,")
	writeLines("among others:\n")
	writeLines("---> One variable for each tumor, linking genomic events to the regulators they modulate,")
	writeLines("according to DIGGIT. Variables are named with the standard TCGA tumor type acronyms")
	writeLines("(blca, brca, etc.) followed by the string 'DG'. Each variable A is a list with two elements:")
	writeLines(" *  A[[1]] is a Nx3 matrix, where N is the number of TFs in the DIGGIT file, i.e., the") 
	writeLines("    number of unique genes B found in the second DIGGIT file column. Each row R in that") 
	writeLines("    matrix corresponds to a single TF B and contains the following data:")
	writeLines("    ** A[[1]][R,1] is the Gene ID of gene B.")
	writeLines("    ** A[[1]][R,2] is the index corresponding to B in the list A[[2]] -- this will be")
	writeLines("       explained below.")
	writeLines("    ** A[[1]][R,3] is the number of genomic events reported by DIGGIT to modulate B.")
	writeLines(" *  A[[2]] is a list with N items, again each corresponding to a TF from the DIGGIT file.")
	writeLines("    The R-th element of the list contains information related to the TF B for which") 
	writeLines("    A[[1]][,2] = R. The R-th list element A[[2]][[R]] is a matrix with dimensions M x 6")
	writeLines("    where M is the number of genomic events reported by DIGGIT to modulate B. Each row S") 
	writeLines("    in that matrix corresponds to a genomic event XXX_C modulating B, where XXX is one of") 
	writeLines("    SNV, AMP, DEL, and denotes the type of genetic event (single nucleotide variation, ")
	writeLines("    amplification, deletion) and C is the Gene ID of the gene harboring the genomic event.")
	writeLines("    The S-th row contains the following data:")
	writeLines("    ** A[[2]][[R]][S, 1] is the gene ID of gene C.")
	writeLines("    ** A[[2]][[R]][S, 2] is an integer representing the type of genomic event (1 = SNV, ")
	writeLines("       2 = AMP, 3 = DEL).")
	writeLines("    ** A[[2]][[R]][S, 3] is the 'areaP' p-value for the genomic event from the DIGGIT file.")
	writeLines("    ** A[[2]][[R]][S, 4] is the 'preppiP' p-value for the genomic event from the DIGGIT file.")
	writeLines("    ** A[[2]][[R]][S, 5] is the 'cindyP' p-value for the genomic event from the DIGGIT file.")
	writeLines("    ** A[[2]][[R]][S, 6] is the 'combinedP' p-value for the genomic event from the DIGGIT file.")
	writeLines(" A few more things to note about variable A:")
	writeLines(" * The matrix A[[1]] is ordered according to number of modulating events (i.e., the number")
	writeLines("   stored in the 3rd column of the table), with the largest entry first.") 
	writeLines(" * The targets contained in the matrix A[[2]][[R]] are listed in decreasing order of")
	writeLines("   significance, based on 'combinedP', with the lowest P-value listed first.\n")
	writeLines("---> A variable named 'modsPerTumor' which is a named vector listing the number of unique modulator")
	writeLines("genes harboring genomic modifications that affect the activity of VIPER regulators, per DIGGIT.\n")
	writeLines("---> A variable named 'eventsPerTumor' which is a named vector listing the number of DIGGIT")
	writeLines("associations reported for each tumor.\n")
	writeLines("---> A variable named 'altsPerTumor' which is a named vector listing the number of all")
	writeLines("unique genomic alterations reported for each tumor.\n")
	writeLines("---> A variable named 'top100FET' which is a list with one elements for each tumor")
	writeLines("(list elements are named, the the command names(top100FET) to see the names).")
	writeLines("The list element top100FET[[i]] corresponods to the i-th tumor and contains a list with one")
	writeLines("element for each of the top most modulated regulators in the tumor, up to 100 regulators. E.g.,")
	writeLines("the list R = top100FET[[\"blcaDG\"]] will contain data for the M regulators blcaDG[[1]][1:M, 1].")
	writeLines("The M list elements in R are named according to the gene IDs of the regulators, run the ")
	writeLines("command names(R) to see them. Each R[[j]] (corresponding to the j-th regulator, say G) is a")
	writeLines("matrix with two columns. Each row in the matrix corresponds to a regulator B other than G such ")
	writeLines("that (brace for this): if M(B) and M(G) are the set of modulator genes harboring DIGGIT events")
	writeLines("that affect, respectively, B and G then the intersection of M(B) and M(R) is significant, as ")
	writeLines("determined by Fisher's exact test (hence the name FET), at a level of p < 0.01 (adjusted).\n")
	writeLines("---> A variable named 'globModMatrix' listing the number of TFs that each modulator gene")
	writeLines("modulates in each tumor. It is a N x 20 matrix M with one column per tumor and one row per")
	writeLines("modulator gene. Rows are named using gene IDs and columns are named according to tumor names,")
	writeLines("using the format 'brcaDG'. M[i, j] is the number of genes modulated by gene rownames[i] in")
	writeLines("tumor colnames(j).")
	writeLines("")
}

# *****************************************************************************
# Helper function, maps the string form of a genomic event into an integer, as 
# follows:
# 
# * "SNV" --> 1
# * "AMP" --> 2
# * "DEL" --> 3 
# * Everything else --> 0
# *****************************************************************************

diggitEventType <- function(str){
	if (str == "SNV")
		return(1)  
	if (str == "AMP")
		return(2)
	if (str == "DEL")
		return(3)
	return(0)
}


# *****************************************************************************
# Summary information of DIGGIT results for a specific TF
#
# Arguments are:
# * geneID: 	The gene id of the query TF gene.
# * eventCountByTumor:	eventCountByTumor[i] is the total number of modulator
#				events for the tumor varNamesDG[i].
# * norm:		FALSE, by default. If TRUE, the result data frame will contain
#				one extra column, as described below.
# * varDG:		Character vector containing the names of the DIGGIT data objects
#				to use, in the format "gbmDG".
#
# Returns a data frame with 20 columns (one for each tumor type) and 3 (or 4)
# columns:
# * Column 1 = tumor name
# * Column 2 = number of DIGGIT genomic events modulating the TF in that tumor
# * Column 3 (optional) = if norm is TRUE, a normalized version of the number of 
#				genomic events shown in column 2. The normalization takes into 
#				account the total number of genomic events reported for a given 
#				tumor.
# * Column 4 = ranking (according to number of modulating DIGGIT genomic events) 
#	       of the TF in that tumor.
# The result data frame is ordered according to column #4.
# *****************************************************************************
tfDGsummary <- function(geneID, eventCountByTumor, norm = FALSE, varDG = varNamesDG){
	NORM_FACTOR = min(eventCountByTumor)
	tumName = vector(length = length(varDG), mode = "character")
	numEvents = vector(length = length(varDG), mode = "integer")
	normNumEvents = vector(length = length(varDG), mode = "integer")
	ranking = vector(length = length(varDG), mode = "integer")
	for (i in 1:length(varDG)){
		tumName[i] = varDG[i]
		if(length(which(get(varDG[i])[[1]][,1] == geneID)) > 0){
			numEvents[i] = get(varDG[i]	)[[1]][which(get(varDG[i])[[1]][,1] == geneID), 3]
			ranking[i] = which(get(varDG[i])[[1]][,1] == geneID)
			normNumEvents[i] = round(NORM_FACTOR*(numEvents[i] / eventCountByTumor[i]),0)
		}
		else{
			numEvents[i] = 0
			ranking[i] = 10000
			normNumEvents[i] = 0
		}
	}
	indOrder = order(ranking)
	tumName = tumName[indOrder]
	numEvents = numEvents[indOrder]
	ranking = ranking[indOrder]
	normNumEvents = normNumEvents[indOrder]
	if (norm)
		return(data.frame(tumName, numEvents, normNumEvents, ranking))
	else
		return(data.frame(tumName, numEvents, ranking))
}


# *****************************************************************************
# Return the list of DIGGIT events across all tumors for a given TF
#
# Arguments are:
# * tfGeneID: 	The gene id of the query TF gene.
# * modGeneID:	If NULL, return genomic events for all modulators. Otherwise
#				return genomic events only for the modulator with this gene ID.
# * filterNulls:	If TRUE, remove NULL entries from the results list
# * idsOnly:	If TRUE, return only the gene IDs of the modulators.
# * varDG:		Vector containing the names of DIGGIT data objects, in the format
#				"gbmDG".
#
# Returns a list with 20 entries (if filterNulls == FALSE; possibly fewer 
# otherwise), one for each tumor. The i-th entry corresponds to the tumor 
# varNamesDG[i] and, if idsOnly == FALSE, contains the value:
#		get(varNamesDG[i])[[2]][[j]] 
# where j is the entry where the DIGGIT events for gene "tfGeneID" are stored. If
# idsOnly == TRUE, the i-th entry contains the vector:
#		unique(get(varNamesDG[i])[[2]][[j]][,1]
# i.e., the list of unique genes modulating the query TF. List entries are named and  
# can be retrieved not only by index but also by providing a tumor string from 
# varNamesDG.
# *****************************************************************************
tfDGDetails <- function(tfGeneID, modGeneID = NULL, filterNulls = FALSE, idsOnly = FALSE, 
		varDG = varNamesDG){
	results = vector("list", length(varDG))
	names(results) = varDG
	for (i in 1:length(varDG)){
		sumr = get(varDG[i])[[1]]
		detl = get(varDG[i])[[2]]
		if (length(which(sumr[,1] == tfGeneID)) > 0){
			ind = sumr[which(sumr[,1] == tfGeneID), 2]
			if (is.null(modGeneID))
				if (idsOnly)
					results[[i]] = unique(detl[[ind]][,1])
				else
					results[[i]] = detl[[ind]]
			else if (length(which(detl[[ind]][,1] == modGeneID)) > 0)
				if (idsOnly)
					results[[i]] = modGeneID
				else
					results[[i]] = detl[[ind]][detl[[ind]][,1] == modGeneID,,drop=FALSE]
		}
	}
	if (filterNulls)
		return(Filter(Negate(is.null), results))
	return(results)
}

# *****************************************************************************
# Return a vector with counts of unique modulator genes harboring genomic events
# in each tumor in 'varNamesDG'. 
#
# The i-th entry lists the number of unique genes harboring genomic variants 
# (per DIGGIT) in the i-th tumor. The vector has named entries, named after the 
# strings in vadNamesDG. If sortem is FALSE then the ordering of the entries is 
# the same found in varNamesDG.	If sortem = TRUE then the entries are in decreasing 
# order of the vector values.
# *****************************************************************************
globalCountOfModulators <- function(varNamesDG, sortem = TRUE){
	uniqueIDs = vector("integer", length(varNamesDG))
	names(uniqueIDs) = varNamesDG
	for (i in 1:length(varNamesDG)){
		L = sum(get(varNamesDG[i])[[1]][,3])
		v = vector("integer", L)
		k = 0
		for (j in 1:length(get(varNamesDG[i])[[1]][,2])){
			v[(k+1):(k+length(get(varNamesDG[i])[[2]][[j]][,1]))] = get(varNamesDG[i])[[2]][[j]][,1]
			k = k + length(get(varNamesDG[i])[[2]][[j]][,1])
		}
		uniqueIDs[i] = length(unique(v))
	}
	if (sortem)
		uniqueIDs = sort(uniqueIDs, decreasing = TRUE)
	return(uniqueIDs)
}



# *****************************************************************************
# Return a vector with counts of the DIGGIT genomic events in each tumor. 
#
# The i-th entry lists the number of genomic events (per DIGGIT) in the tumor
# whose acronum (in the form "gbmDG") is found in the i-th entry of vector varDG.
# The result vector has named entries, named after the strings in varDG. If sortem 
# is FALSE then the ordering of the entries is the same found in varDG.	
# If sortem = TRUE then the entries are in decreasing order of the vector values.
# *****************************************************************************
globalCountOfEvents <- function(sortem = TRUE, varDG = varNamesDG){
	if (sortem)
		return(sort(sapply(varDG, function(e){return(sum(get(e)[[1]][,3]))}), decreasing = TRUE))
	return(sapply(varDG, function(e){return(sum(get(e)[[1]][,3]))}))
}

# *****************************************************************************
# Return the number (or list) of unique modulator genes harboring genomic events
# (per DIGGIT) for a query tumor
#
# Arguments:
# * dg_data:	the DIGGIT data object for the tumor of interest (i.e., one of
#				the variables whose name is found in the vector varDG). Or 
#				the string name of the query tumor, in the format "gbmDG". Or
#				NULL: in that case, the number (of list) of all modulator genes
#				across all tumors in varDG is returned.
# * getVector:	Determines is a count (getVector = FALSE) or a vector (getVector =
#				TRUE) of the modulator genes is returned
# * varDG:		a vector contains the names of the DIGGIT data objects to use 
#				in the case where dg_data == NULL
# *****************************************************************************
countUniqModGenes <- function(dg_data = NULL, getVector = FALSE, varDG = varNamesDG){
	tumorL = list()
	if (is.null(dg_data))
		tumorL = lapply(varDG, function(e){return(get(e))})
	else if (typeof(dg_data) == "character")
		tumorL[[1]] = get(dg_data)
	else tumorL[[1]] = dg_data
	
	tally = rep(0, max_geneId)
	for(index in 1:length(tumorL)){
		dg_data = tumorL[[index]]
		L = length(dg_data[[2]])
		for (i in 1:L){
			x = dg_data[[2]][[i]][,1]
			for (j in 1:length(x))
				tally[x[j]] = 1
		}
	}
	if (getVector)
		return(which(tally > 0))
	else
		return(length(which(tally > 0)))
}


# *****************************************************************************
# Sorted list of genes harboring modulating events for given gene ID
#
# Arguments:
# * geneID:	the ID of the query TF 
# * varDG:	character vector containing the names of DIGGIT objects, in the form
#			"gbmDG". Only tumors in this vector are examined.
#
# Returns a vector V where names(V) contains gene IDS and V[i] is the number of 
# tumors from varDG where the gene names(V)[i] is reported to harbor at least one 
# genomic event that modulates the TF. Values V[i] are in decreasing order.
# *****************************************************************************
geneTableOfModulators <- function(geneID, varDG = varNamesDG){
	y = tfDGDetails(geneID)
	union = unique(y[[1]][,1])
	for (i in 2:length(varDG))
		union = append(union, unique(y[[i]][,1]))
	return(sort(table(union), decreasing = TRUE))
}



# *****************************************************************************
# Sorted list of TFs that are modulated by alterations of a modulator gene
#
# Arguments:
# * modGeneID:	the ID of the query modulator
# * varDG:		the query tumor (either a string or the actual variable).
#
# Returns a vector V where names(V) contains gene IDS of TFs reported by DIGGIT
# to be modulated by one or more events in 'modGeneID' in the query tumor. V[i] 
# is the "combinedP" p-value for the corresponding modulation event. Values V[i] 
# are in increasing order.
# *****************************************************************************
getModulatedRegulators <- function(modGeneID, varDG){
	if (typeof(varDG) == "character")
		varDG = get(varDG)
	tmp = matrix(nrow = length(varDG[[2]]), ncol = 2)
	tmp[,] = 0
	rownames(tmp) = varDG[[1]][,1]
	for (i in 1:length(varDG[[1]][,1])){
		x = varDG[[2]][[varDG[[1]][i,2]]]
		if (length(which(x[,1] == modGeneID)) > 0){
			tmp[i,1] = 1
			tmp[i,2] = min(x[which(x[,1] == modGeneID), 6])
		}
	}
	tmp = tmp[tmp[,1] == 1,]
	res = tmp[,2]
	names(res) = rownames(tmp)
	return(sort(res))
}

# *****************************************************************************
# Count the number of TFs that each modulator gene modulates in each tumor
#
# Arguments:
# * unique:		If TRUE, a modulating effect of gene B on TF A in a tumor is 
#				counted only once, even if there are multiple alterations of B 
#				that modulate A. Otherwise, each individual alteration contributes
#				1 to the total count
# * varNamesDG: A vector containing the names of the DIGGIT data objects, one per
#				tumor to be analyzed
#
# Returns a N x length(varNamesDG) matrix M with one column per tumor and one row
# per modulator gene. Rows are named using gene IDs anc columns are named according 
# to the tumor names in 'varNamesDG'. M[i, j] is the number of genes modulated
# by gene rownames[i] in tumor colnames(j).
# *****************************************************************************
globalModulatorMatrix <- function (unique = FALSE, varDG = varNamesDG) {
	# Get a list of all the modulator genes encountered across all tumors
	modGenes = countUniqModGenes(getVector = TRUE, varDG = varDG)
	results = matrix(0, nrow = length(modGenes), ncol = length(varDG))
	rownames(results) = modGenes
	colnames(results) = varDG
	for (index in 1:length(varDG)){
		x = modulatorCounts(get(varDG[index]), unique = unique)
		for (i in 1:length(x))
			results[names(x)[i], index] = x[i] 
	}
	# Fix column names before returning by removing capitalized suffixes from the 
	# tumor acronyms and by ordering columns alphabetically
	colnames(results) = gsub("[A-Z]+", "", colnames(results))
	return (results[,sort(colnames(results))])
}


# *****************************************************************************
# Assesses the enrichment in co-modulating genes between a query TF and all
# other VIPER TFs in a given tumor
#
# Arguments:
# * geneID:		the ID of the query TF
# * varDG:		the DIGGIT data object for the tumor of interest (i.e., one of
#				the variables whose name is found in the vector varNamesDG).
# * correct:	If TRUE then the p-values from the Fisher's exact tests are 
#				Bonferroni-corrected
# * allCount:	For the computation of the test we need to know the total number
#				of genes harboring genomic events in the tumor of interest. If
#				non-zero, this argument provides this number. Otherwise the
#				number will be computed as part of the function call (which 
#				will affect performance a bit).
#
# Returns a vector V where names(V) contains gene IDS of all the TFs in 
# varNamesDG[[1]][,1] and V[i] is the p-value of Fisher's exact test that compares
# the size of the intersection between two set of genes: the set A of genes 
# harboring modulating events for 'geneID' and the set B of genes harboring
# modulating events for the TF names(V[i]).
# *****************************************************************************
coModulatedRegulators <- function(tfGeneId, varDG, correct = FALSE, allCount = 0){
	results = NULL
	tfGeneId = as.integer(tfGeneId)
	if (length(which(varDG[[1]][,1] == tfGeneId)) == 0)
		return(results)
	
	tfCount = length(varDG[[1]][,1])
	if (allCount == 0)
		allCount = countUniqModGenes(varDG)
	tfModGenes = unique(varDG[[2]][[toString(tfGeneId)]][,1])
	results = rep(1, tfCount)
	names(results) = varDG[[1]][,1]
	for (i in 1:tfCount){
		if (varDG[[1]][i,1] == tfGeneId)
			next
		otherTfGeneId = varDG[[1]][i,1]
		otherTfModGenes = unique(varDG[[2]][[toString(otherTfGeneId)]][,1])
		# Set up the counts for Fisher's exact test
		a = length(intersect(tfModGenes, otherTfModGenes))
		b = length(setdiff(tfModGenes, otherTfModGenes))
		c = length(setdiff(otherTfModGenes, tfModGenes))
		d = allCount - length(union(tfModGenes, otherTfModGenes))
		contigencyMat = matrix(c(a, c, b, d), 2, 2)
		results[i] = fisher.test(contigencyMat)$p.value
	}
	
	if (correct)
		for (i in 1:length(results))
			results[i] = min(tfCount*results[i], 1)
	return (sort(results))
}


# *****************************************************************************
# Creates (and optionally saves to disk) graphs for a bunch of variables of
# interest
# 
# Arguments:
# * which:	Indicates which graph to plot.
# * save:	Boolean, indicates if the plot should be saved to file.
# * fName:	If "save" is TRUE, this is the file name to user for storing the 
#			image. If this is NULL then use as file name the value of 'which'.  
# * varsDG: A vector containing the names of the DIGGIT data variables,
#			one per each tumor.
# * width, height, res:	Image specification arguments to control width, height,
#			and resolution.
# * params:	A list contains subfunction-specific arguments. This is up to 
#			each subfunction to interpret as needed.
# *****************************************************************************
makeGraphs <- function (which = "mods_events", save = FALSE, fName = NULL,
		varsDG = varNamesDG, width = 800, height = 800, res = 130, params=NULL){
	if (save){
		if (is.null(fName))
			fName = paste(which, ".png", sep="")
		png(fName, width = width, height = height, res = res)
	}
	
	# ******************** which = regulon_sizes  *****************************
	# For each interactome, draw a boxplot for the sizes of all regulons. 
	# 
	# ARGUMENTS
	# * params[[1]]: list of interactomes to plot, in the form of a vector of
	#		character strings, each string providing the variable name of an
	#		interactome. E.g., params[[1]] = c("Lung", "tcga_read")
	#		If params == NULL, plot bars for all interactomes.
	if (which == "regulon_sizes"){
		if (is.null(params))
			netNames = varNames
		else
			netNames = params[[1]]
		regSizes = sapply(netNames, function(x){get(x)[[1]][,3]})
		par(mfrow=c(1, 1))
		par(mar=c(8,4,4,5)+.1)
		title = "Distribution of regulon sizes per interactome"
		boxplot(regSizes, main = title, las=2)
	}
	
	
	# ******************** which = regulon_size  *****************************
	# Given a hub gene and a set of interactomes, generate a bar plot showing
	#	the regulons size for the gene in each interactome. 
	# 
	# ARGUMENTS
	# * params[[1]]: Query gene, as a gene symbol or entrez id.
	# * params[[2]]: list of interactomes to plot, in the form of a vector of
	#		character strings, each string providing the variable name of an
	#		interactome. E.g., params[[1]] = c("Lung", "tcga_read")
	#		If params == NULL, plot bars for all interactomes.
	# * params[[3]]: a numeric MI threshold. When counting the number of interactions
	#		in a regulon, only inreractions with MI >= params[[3]] will be considered.
	# * params[[4]]: a numeric p-value threshold. When counting the number of interactions
	#		in a regulon, only inreractions with p-val <= params[[4]] will be considered.
	# * params[[5]]: a boolean value; if TRUE, barplots will be plotted in 
	#		decreasing order of the regulon size. Otherwise, the will be plotted
	#		in the order specified in params[[2]]
	if (which == "regulon_size"){
		if (is.null(params) || is.na(params) || length(params) != 5){
			writeLines("Incorrect 'params[[]]' values")
			return()
		}
		hub = as.entrezId(params[[1]])
		if (is.null(params[[2]]))
			netNames = varNames
		else
			netNames = params[[2]]
		mi_thr = params[[3]]
		pval_thr = params[[4]]
		order_nets = params[[5]]
		
		reg_sizes = sapply(netNames, function(n){return(length(getRegulon(hub, net=n, mi=mi_thr, pval = pval_thr)))})
		if (order_nets)
			reg_sizes = sort(reg_sizes, decreasing = TRUE)
		par(mfrow=c(1, 1))
		par(mar=c(8,4,4,5)+.1)
		title = paste("Regulon sizes per interactome for gene ", hub, " (", entrezIDtoSymbol(hub), ")", sep="")
		ylim <- c(0, 1.1*max(reg_sizes))
		bp = barplot(reg_sizes, las=2, main = title, ylab="Regulon size", ylim=ylim)
		text(x = bp, y = reg_sizes, label = reg_sizes, pos = 3, cex = 0.8, col = "red")
	}
	
	
	# ******************** which = net_sizes  *****************************
	# Create a bar plot with one bar per interactome, each bar representing
	# the number of interactions in the relevant interactome. Optionally, 
	# graw below it a second bar plot, showing either the number of all genes
	# or the number of hub genes in each interactome.
	# 
	# ARGUMENTS
	# * params[[1]]: list of intectomes to plot, in the form of a vector of
	#		character strings, each string providing the variable name of an
	#		interactome. E.g., params[[1]] = c("Lung", "tcga_read")
	#		If params == NULL, ot params[[1]] == NULL plot bars for all interactomes.
	# * params[[2]]: one of two strings, "all" or "hubs". If provided it specifies
	#		the data to use for the second bar plot.
	if (which == "net_sizes"){
		second_plot = FALSE
		if (is.null(params) | is.null(params[[1]]))
			netNames = varNames
		else
			netNames = params[[1]]
		if (!is.null(params) & length(params) == 2){
			par(mfrow=c(2, 1))
			second_plot = TRUE
		}
		else{
			par(mfrow=c(1, 1))
			
		}
		
		nets = sort(netSizes[netNames])
		par(mar=c(8,4,4,5)+.1)
		title = "Number of interactions per interactome (sorted, smallest to largest)"
		colors = rep("blue", length(nets))
		colors[grep("tcga", names(nets))] = "red"
		barplot(nets, xlab="", main=title, las = 2, col=colors)
		legend('topleft', legend = c("TCGA", "GTEx"), col=c("red", "blue"), lty=1)
		
		if (second_plot){
			if (params[[2]] == "hubs"){
				t = sapply(netNames, function(x){getInteractomeGenes(x, hubs_only=TRUE)})
				title = "Number of hub genes per interactome"
			}
			else if (params[[2]] == "all"){
				t = sapply(netNames, getInteractomeGenes)
				title = "Number of genes per interactome"
			}
			else
				return()
			t=t[names(sort(netSizes[netNames]))]
			barplot(t, main=title, las=2)
		}
	}
	
	# ******************** which = mods_events  *****************************
	# Plot 2 number series, each comprising 20 counts, one for each tumor: 
	# (1) the number of genes harboring DIGGIT events.
	# (2) the total number of repored DIGGIT associations 
	if (which == "mods_events"){
		par(mar=c(5,4,4,5)+.1)
		plot(modsPerTumor, xlab = "Tumors", ylab = "# of Genes harboring DIGGIT mutations", xaxt="n", yaxt="n")
		axis(1, at=seq(1,length(varsDG)), labels=gsub("DG", "", names(modsPerTumor)), las = 2)
		axis(3, at=seq(1,length(varsDG)), labels=gsub("DG", "", names(modsPerTumor)), las = 2)
		axis(2, at=seq(1, max(modsPerTumor), length.out=5), 
				labels = as.character(round(seq(1, max(modsPerTumor), length.out=5))), 
				col = "black", las=2)
		par(new=TRUE)
		plot(eventsPerTumor[names(modsPerTumor)],type="p",col="red",xaxt="n",yaxt="n",xlab="",ylab="")
		axis(4, at=seq(1, max(eventsPerTumor), length.out=5), 
				labels = as.character(round(seq(1, max(eventsPerTumor), length.out=5))), 
				col = "red", las=2, cex=0.75)
		abline(v=seq(1:length(varsDG)), lty=3)
		legend('topright', legend = c("# Genes", "# Assoc"), col=c("black", "red"), lty=1)
		mtext("# of DIGGIT associations", side=4, line = 4)
	}
	
	# ********************* which = events_regs ********************************-
	# Plot 2 number series, each comprising 20 counts, one for each tumor: 
	# (1) the number of VIPER regulators
	# (2) the total number of reported DIGGIT associations 
	if (which == "events_regs"){
		par(mar=c(5,4,4,5)+.1)
		#series = sort(sapply(varsDG, function (e) {return(length(get(e)[[1]][,1]))}), decreasing = TRUE)
		plot(regsPerTumor, xlab = "Tumors", ylab = "# of VIPER regulators", xaxt="n", yaxt="n")
		axis(1, at=seq(1,length(regsPerTumor)), labels=gsub("DG", "", names(regsPerTumor)), las = 2)
		axis(3, at=seq(1,length(regsPerTumor)), labels=gsub("DG", "", names(regsPerTumor)), las = 2)
		axis(2, at=seq(1, max(regsPerTumor), length.out=5), 
				labels = as.character(round(seq(1, max(regsPerTumor), length.out=5))), 
				col = "black", las=2)
		par(new=TRUE)
		plot(eventsPerTumor[names(regsPerTumor)],type="p",col="red",xaxt="n",yaxt="n",xlab="",ylab="")
		axis(4, at=seq(1, max(eventsPerTumor), length.out=5), 
				labels = as.character(round(seq(1, max(eventsPerTumor), length.out=5))), 
				col = "red", las=2, cex=0.75)
		abline(v=seq(1:length(varsDG)), lty=3)
		legend('topright', legend = c("# Regs", "# Assoc"), col=c("black", "red"), lty=1)
		mtext("# of DIGGIT associations", side=4, line = 4)
	}
	
	# ********************* which = median_events ********************************
	if (which == "median_events"){
		par(mar=c(5,5,4,4)+.1)
		plot(eventsPerTumor, xlab = "Tumors", ylab = "", xaxt="n", yaxt="n")
		mtext("# of DIGGIT associations", side=2, line = 4)
		axis(1, at=seq(1,length(eventsPerTumor)), labels=gsub("DG", "", names(eventsPerTumor)), las = 2)
		axis(3, at=seq(1,length(eventsPerTumor)), labels=gsub("DG", "", names(eventsPerTumor)), las = 2)
		axis(2, at=seq(1, max(eventsPerTumor), length.out = 5), 
				labels = as.character(round(seq(1, max(eventsPerTumor), length.out = 5))), 
				las=2, cex=0.75)
		par(new=TRUE)
		medians = sapply(names(eventsPerTumor), function(e) {return(median(get(e)[[1]][,3]))})
		plot(medians,type="p",col="red",xaxt="n",yaxt="n",xlab="",ylab="")
		axis(4, at=seq(0, max(medians), ceiling(max(medians)/5)), 
				labels = as.character(seq(0, max(medians), ceiling(max(medians)/5))), 
				col = "red", las=2, cex=0.75)
		abline(v=seq(1:length(medians)), lty=3)
		legend('topright', legend = c("# Assoc", "Median/gene"), col=c("black", "red"), lty=1)
		mtext("Median # of associations per regulator", side=4, line = 3)
		
	}
	
	# ********************* which = viper_regs ********************************
	# Plot the number of VIPER regulators for each tumor
	if (which == "viper_regs"){
		series = sort(regsPerTumor, decreasing = TRUE)
		plot(series, xlab = "Tumors", ylab = "# of VIPER regulators", xaxt="n")
		axis(1, at=seq(1,length(series)), labels=gsub("DG", "", names(series)), las = 2)
		axis(3, at=seq(1,length(series)), labels=gsub("DG", "", names(series)), las = 2)
		abline(v=seq(1:length(series)), lty=3)
	}
	
	# ********************* which = events_per_gene ********************************
	# Generate one plot for each tumor, each showing the density of the number of DIGGIT events
	# affecting every gene. Graphs are plotted in in decreasing order of the total number 
	# of DIGGIT associations
	if (which == "events_per_gene"){
		par(mfrow=c(ceiling(sqrt(length(varsDG))), ceiling(sqrt(length(varsDG)))))
		orderToDraw = names(eventsPerTumor)
		for (i in 1:length(orderToDraw)){
			xlab = paste("N = ", length(get(orderToDraw[i])[[1]][,3]), 
					"    Median = ", median(get(orderToDraw[i])[[1]][,3]))
			plot(density(get(orderToDraw[i])[[1]][,3]), main=gsub("DG", "", orderToDraw[i]), xlab=xlab)
		}
	}
	
	# ********************* which = enriched_activity ********************************
	# Generate one plot for each tumor, showing the enrichment of top-modulated 
	# regulators in highly active regulators (red line), highly suppressed modulators 
	# (blue line), and several random choices of modulator (black lines) 
	if (which == "enriched_activity"){
		par(mfrow=c(ceiling(sqrt(length(varsDG))), ceiling(sqrt(length(varsDG)))))
		how_many_random = 5
		for (i in 1:length(varsDG)){
			viper = get(gsub("DG", "VP", varsDG[i]))[[1]]
			plot(incrementalFETs(varsDG[i], viper, mode="pos")$pvals, t="l", 
					main = gsub("DG", "", varsDG[i]), ylab = "-Log(pval)", xlab="", col="red")
			lines(incrementalFETs(varsDG[i], viper, mode="neg")$pvals, t="l", col="blue")
			for (j in 1:how_many_random)
				lines(incrementalFETs(varsDG[i], viper, mode="pos", random=TRUE)$pvals, t="l")
			
		}
	}
	
	# ********************* which = events_boxplots ********************************
	# Generate one box plots for each tumor, for the counts of DIGGIT events
	# per VIPER regulator.
	if (which == "events_boxplots"){
		par(mar=c(5,5,4,7)+.1)
		medians = sort(sapply(names(eventsPerTumor), function(e) {return(median(get(e)[[1]][,3]))}))
		L = sapply(varsDG, function(e){return(get(e)[[1]][,3])})
		L = L[names(medians)]
		names(L) = gsub("DG", "", names(L))
		title = paste("Counts of DIGGIT events associated with VIPER regulators, ordered by median.\n",
				"(red line: total # of DIGGIT events per tumor)")
		boxplot(L, las=2, main=title, xlab="", ylab="")
		par(new=TRUE)
		plot(eventsPerTumor[names(medians)], type="l", col = "red", axes=FALSE, xlab="", 
				ylab = "Median counts of regulator events")
		axis(4, at=seq(1, max(eventsPerTumor), length.out = 5), 
				labels = as.character(round(seq(1, max(eventsPerTumor), length.out = 5))), 
				col="red", las=2, cex=0.75)
		abline(v=seq(1:20), lty=3)
		mtext("Total # of DIGGIT associations", side=4, line = 5)
	}
	
	
	# ********************* which = events_barplots ********************************
	# Generate one bar plots for each tumor, showing the counts of the 3 
	# types of DIGGIT events: SNVs, AMPs, DELs
	if (which == "events_barplots"){
		par(mar=c(5,5,4,7)+.1)
		tmp = sapply(varsDG, function(x){return(rowSums(sapply(get(x)[[2]], 
											function(e){return(c(sum(e[,2] == 1), sum(e[,2] == 2), sum(e[,2] == 3)))})))})
		title = "Counts of DIGGIT events per Tumor"
		ylab = "# of DIGGIT associations"
		mp = barplot(tmp[,names(sort(eventsPerTumor, decreasing=TRUE))], col = c("black", "red", "blue"), 
				xaxt="n", yaxt="n",	legend.text = c("SNV", "AMP", "DEL"), main = title)
		axis(1, at=mp, labels=gsub("DG", "", names(sort(eventsPerTumor, decreasing=TRUE))), las = 2)
		axis(2, at=seq(1, max(eventsPerTumor), length.out = 5), 
				labels = as.character(round(seq(1, max(eventsPerTumor), length.out = 5))), 
				col="black", las=2)
	}
	
	# ********************* which = density_per_event ********************************
	# Generate one plot for each tumor, each showing the 3 density graphs, one each for 
	# the -log(combinedP) value of DIGGIT associations involving either SNV, AMP, or DEL 
	# type genomic events
	if (which == "density_per_event"){
		par(mfrow=c(ceiling(sqrt(length(varsDG))), ceiling(sqrt(length(varsDG)))))
		orderToDraw = names(eventsPerTumor)
		for (i in 1:length(orderToDraw)){
			snv = unlist(sapply(get(orderToDraw[i])[[2]], function(x){return(x[x[,2] == 1,6])}))
			amp = unlist(sapply(get(orderToDraw[i])[[2]], function(x){return(x[x[,2] == 2,6])}))
			del = unlist(sapply(get(orderToDraw[i])[[2]], function(x){return(x[x[,2] == 3,6])}))
			xlab = paste("SNV# = ", length(snv), " AMP# = ",  length(amp), "  DEL# = ",  length(del))
			plot(density(-log(snv)), ylim=c(0, 0.5), main=gsub("DG", "", orderToDraw[i]), xlab=xlab)
			lines(density(-log(amp)), col="red")
			lines(density(-log(del)), col="blue")
		}
	}
	
	# ********************* which = chrom_diff *************************************
	# Generate one histograms for each tumor, showing the frequency of the the 
	# values:
	#			chromosome(mod) - chromosome(reg)
	# for each DIGGIT event involving a regulator gene 'reg' and a modulator gene 
	# 'mod'. Our goal is to observe if there is any systematic (and thus, possibly
	# alarming) preference for modulators targeting regulators in the same chromosome.
	# This, e.g., could indicate that what DIGGIT picks up are associations between
	# modulator-regulator pairs where the modulator and the regulator are part of the 
	# same genomic aberration.	
	if (which == "chrom_diff"){
		par(mfrow=c(ceiling(sqrt(length(varsDG))), ceiling(sqrt(length(varsDG)))))
		# Print plots in increasing order of number of DIGGIT events
		names = names(sort(eventsPerTumor))
		chr_map = rep(0, max_geneId)
		for (index in 1:length(names)){
			dg_data = get(names[index])
			t = vector(mode="numeric", length=sum(dg_data[[1]][,3]))
			k = 1
			for (i in 1:length(dg_data[[1]][,1])){
				rid = dg_data[[1]][i,1]
				r_ind = dg_data[[1]][i,2]
				mods = dg_data[[2]][[r_ind]]
				rid_chr = chr_map[rid]
				if (!is.na(rid_chr) && rid_chr == 0){
					rid_chr = getChromosome(rid)
					chr_map[rid] = rid_chr
				}
				for (j in 1:length(mods[,1])){
					mid_chr = chr_map[mods[j,1]]
					if (!is.na(mid_chr) && mid_chr == 0){
						mid_chr = getChromosome(mods[j,1])
						chr_map[mods[j,1]] = mid_chr
					}
					t[k] = mid_chr - rid_chr
					k = k+1
				}
			}
			title = paste(gsub("DG", "", names[index]),"=", length(unlist(t)))
			x_lab = "Chr(mod)-Chr(reg)"
			y_lab = "# DIGGITs"
			hist(unlist(t), main = title, ylab = y_lab, xlab=x_lab)
		}
	}
	# ********************* which = distr_high_positive_vipers ************************
	# We create a plot demonstrating that regulators that are highly active across
	# multiple tumors have regulons which are well conserved across tumors. Specifically,
	# within each tumor, we sort all VIPER regulators in order of activity, from highest
	# to lowest and select the 150 most active regulators in each tumor (as the activity of
	# a regulator in a tumor we use its median VIPER activity across all samples). We then 
	# combine this data to identify the vector "top" of regulators that appear in the top-150 
	# lists of at least five tumors. For each gene G in "top" we then retrieve the -Log(p)
	# FET enrichment values for the regulon of G from fPairEnrich[[2]], which indicate how
	# well conserved the regulon of G is across all 190 pairs of interactomes. The -Log(p) 
	# values for all genes G in "top" are the aggregated and their distribution is plotted
	# (using a red-colored line). We then repeat this process for (1) the union "tfids" of 
	# all VIPER regulators from all tumors (plotted using a black line) and (2) a number "rep" of 
	# random samples from "tfids" where the size of each sample is same as the number of
	# genes in "top".
	if (which == "distr_high_positive_vipers"){
		tfids = unique(unlist(sapply(varNamesVP, function(x){return(unique(get(x)[[1]][,1]))})))
		xAxis = "Log(p) enrichment score"
		plotTitle = "Density of Log(p) values for the enrichment in common targets among the regulons 
				of VIPER regulators, across all possible pairs of tumors. Black line is density for
				all regulators. Red line is only for the list L of most active regulators in the top-150 list.
				Green lines are for multiple random choices of |L| regulators from the full set."
		plot(density(unlist(sapply(tfPairEnrich[[2]], 
										function(x){
											return(x[, 2])
										}))), xlab = xAxis, main = plotTitle, cex.main=0.8)
		L = lapply(varNamesVP, function(name){
					net = get(name)
					inc = length(unique(net[[2]]))
					x= nrow(net[[1]])/inc
					y = sapply(seq(1,x), function(k){return(median(net[[1]][((k-1)*inc+1):(k*inc), 2]))})
					names(y) = unique(net[[1]][,1])
					return(sort(y, decreasing=T))
				})
		names(L) = varNamesVP
		t_150 = sort(table(unlist(lapply(L, function(x, top=150){return(strtoi(names(x)[1:top]))}))), decreasing=T)
		min = 5
		max = length(varNamesVP)
		res = rep("", length(max:min))
		for (i in max:min){
			t = names(t_150[t_150 == i])
			if (length(t) > 0)
				res[max-i+1] = toString(entrezIDtoSymbol(t))
		}		
		top = geneSymbolToEntrezId(unlist(strsplit(paste(res[(res != "")], collapse=", "), ", ", fixed=T)))
		lines(density(unlist(sapply(tfPairEnrich[[2]], 
										function(x){
											return(x[which(x[,1] %in% top), 2])
										}))), col="red")
		repNum = 200
		for (i in 1:repNum){
			topSample = sample(tfids, length(top))
			lines(density(unlist(sapply(tfPairEnrich[[2]], 
											function(x){
												return(x[which(x[,1] %in% topSample), 2])
											}))), col="green")
		}
		
	}
	
	# ***************************** which = cluster_networks ********************************
	# Cluster interactomes based on "similarity". A variety of distacne metrics is used, aiming
	# to capture how well interactions are conserved among pairs of interactomes. Similarity 
	# is assessed using either all hubs, TFs only, TFs and co-TFs only, or siganling hubs only,
	# thus resulting in four different cluster plots.
	#
	# @FIXME: the setting in params[[3]] is currently relevant only when params[[1]] == "fet".
	# Should implement more generally.
	#
	# ARGUMENTS:
	# params[[1]]: a string. If equal to "fet" then the pre-computed FET p-values stored in
	#		the object "tfPairEnrich" are used for compuitng inter-interactome distances. If
	#		equal to "prob" then the p-values stored in the object "tfPairProb" are used instead.
	#		In all other cases, distances are simply based on the number of shared interactions
	#		between two interactomes.
	# params[[2]]: a string, specifying the type of plot. If "mds", MDS plots are created.
	#		Otherwise, hierarchical clustering dendrograms are drawn.
	# params[[3]]: optional parameter. If provided, it should be a vector of gene symbols
	#		(character strings) or a vector of entrez ids (integers or character strings)
	#		correponding to hub genes (if non-hub genes are included, they are just ignored).
	#		In that case, when computing distances, interactions involving these hub genes are
	#		not take into consideration. in the current implementation, the params[[3]] option
	#		is taken into account only when params[[1]] == "fet".
	if (which == "cluster_networks"){
		if((is.null(params)) || is.na(params) || length(params) < 2 || length(params) > 3)
			stop("Wrong number of parameters.")
		if(params[[1]] == "fet") {
			if (length(params) != 3)
				buf <- calFet(NULL)
			else
				buf <- calFet(params[[3]])
			
			par(mfrow=c(2,2))
			graphTitles = c("Distances using all FET",
					"Distances using FET from TF hubs",
					"Distances using FET from TF and co-TF hubs",
					"Distances using FET from signaling hubs")
			
			for(j in 1:length(buf)) {
				M <- min(buf[[j]])
				buf[[j]] = buf[[j]]-M+1
				
				d = as.dist(buf[[j]])
				
				if(params[[2]] == "mds") {
					fit <- cmdscale(d,eig=TRUE, k=2)
					x_c <- fit$points[,1]
					y_c <- fit$points[,2]
					
					plot(x_c, y_c, main=graphTitles[j],	type="n")
					text(x_c, y_c, labels = labels(d), cex=.7, pos=3)
				} else {
					plot(hclust(d), main=graphTitles[j])
				}
			}
			
		} else if(params[[1]] == "prob") {
			buf <- calNetProb()
			
			par(mfrow=c(2,2))
			graphTitles = c("Distances using the probability from all hubs",
					"Distances using the probability from TF hubs",
					"Distances using the probability from TF and co-TF hubs",
					"Distances using the probability from signaling hubs")
			
			for(j in 1:length(buf)) {
				M <- min(buf[[j]])
				buf[[j]] = buf[[j]]-M+1
				
				d = as.dist(buf[[j]])
				
				if(params[[2]] == "mds") {
					fit <- cmdscale(d,eig=TRUE, k=2)
					x_c <- fit$points[,1]
					y_c <- fit$points[,2]
					
					plot(x_c, y_c, main=graphTitles[j],	type="n")
					text(x_c, y_c, labels = labels(d), cex=.7, pos=3)
				} else {
					plot(hclust(d), main=graphTitles[j])
				}
			}
		} else {
			# The following 4 variables are "expensive" to generate, so do that only once.
			if (!exists("pairWise"))
				assign("pairWise", interactionCounts("pairwise"), envir = globalenv())	
			if (!exists("pairWiseTFonly"))
				assign("pairWiseTFonly", interactionCounts("pairwise", is.tf), envir = globalenv())
			if (!exists("pairWiseTFandCoTF"))
				assign("pairWiseTFandCoTF", interactionCounts("pairwise", function(x){return (is.tf(x) || is.cotf(x))}), 
						envir = globalenv())
			if (!exists("pairWiseSignalOnly"))
				assign("pairWiseSignalOnly", interactionCounts("pairwise", is.sign), envir = globalenv())
			pwiseCounts = list(pairWise, pairWiseTFonly, pairWiseTFandCoTF, pairWiseSignalOnly)
			# par(mfrow=c(1,length(pwiseCounts)))
			par(mfrow=c(2,2))
			graphTitles = c("Distances using all interactions",
					"Distances using only interactions from TF hubs",
					"Distances using only interactions from TF and co-TF hubs",
					"Distances using only interactions from signaling hubs")
			for (j in 1:length(pwiseCounts)){
				x = pwiseCounts[[j]]
				N = nrow(x)
				for (i in 1:N)
					x[i,i]=0
				max = max(x)
				x = max - x
				for (i in 1:N)
					x[i,i]=0
				d = as.dist(x)
				if ((!is.null(params)) && (params[[1]] == "mds")){
					fit <- cmdscale(d,eig=TRUE, k=2)
					x_c <- fit$points[,1]
					y_c <- fit$points[,2]
					# This is a hack; when computing distances using TFs or TFs+co-TFs, the MDS diagram
					# is plotted in the opposite orientation. To make all images look the same, we just
					# flit the x axis.
					# if (j %in% c(2,3))
					#	x_c = -x_c
					plot(x_c, y_c, main=graphTitles[j],	type="n")
					text(x_c, y_c, labels = labels(d), cex=.7, pos=3)
				} else{
					plot(hclust(d), main=graphTitles[j])
				}
			}
		}
	}
	
	
	# ********************* which = regulon_sizes_per_gene_type ************************
	# Create one plot per ARACNe network, each graphing regulon sizes in decreasing order, 
	# with vertical lines corresponding to hub type (TF = red, coTF = blue, signaling = green)
	# type genomic events. ++++NOTE: The ARACNe networks to plot must be passed as the value
	# of the function argument 'varsGS', in the form of a string vector such as:
	# 		c("brca", "blca", "gbm", "lgg")
	if (which == "regulon_sizes_per_gene_type"){
		if (!is.character(varsDG))
			stop()
		nets =  lapply(varsDG, function(x){return(get(x))})
		#par(mfrow=c(ceiling(sqrt(length(varsDG))), ceiling(sqrt(length(varsDG)))))
		par(mfrow=c(ceiling(length(varsDG)/3), min(3, length(varsDG))))
		for (i in 1:length(varsDG)){
			plotReg(nets[[i]])
			title(main = varsDG[i])
		}
	}
	
	
	# ********************* which = exclusive_interactions ************************
	# Process the file generated by a call to the method:
	#			oneOffs("exclusive_interactions", params)
	# which produces a sorted table T for each hub gene R seen it both GTEx and TCGA
	# networks; each table entry corresponds to gene G such that the interaction 
	# (R, G) appears only in GTEx (or TCGA) networks and T[G] is the support of 
	# (R, G), i.e., the number of networks where it is present. Additionally, the call
	# to oneOffs() will generate similar tables for randomly reassigning networks to 
	# the GTEx or TCGA category, to be used as nulls.
	# The code that follows will generate two plots, one for GTEx and the other for
	# TCGA. For GTEx, it will look at interactions involving target G and count the 
	# number of hub genes R such that (1) (R, G) in an interactions that is 
	# exclusive to GTEx, and (2) has support larger than a user-specified value "min".
	# It will the plot the density plot for these counts and overlay the same count 
	# as computed from the randomly reassigned networks. The plot for TCGA is computed
	# in the exact same manner, but counting interactions that are exclusive to GTEx.
	# Our goal here is to identify true interactions (hence the requirement to have
	# large support and to clear MI and p-value thresholds) that are exclusive to 
	# either GTEx or TCGA (hence, putatively plan a role in suppressing or promoting
	# tumorigenesis) and whose regulation is coordinated by an unexpectedly high number
	# of regulators. The characteristics, we postulate, make for ineresting interactions.
	# 
	# The code below expects the following values to be passed in the params[[]] argument:
	# * params[[1]]:	Character string containing the full pathname of the file 
	#			produced by the call to oneOffs.
	# * params[[2]]:	Integer providing the value "min" which will be used as the minimum
	#			support for the interactions to be considered.
	if (which == "exclusive_interactions"){
		
		# ****************** tabulateCommonRegulators **************************
		# Takes as input a list like L_gtex of L_tcga generated by method computeDiffinteractions().
		# Each list member correspond to a hub gene and contains its top exclusive interactions, tabulated
		# according to their support. The methods first identifies interactions that have support 
		# in at least "min" networks. It then examines the targets involved in these interactions
		# and tabulates them according to how many hub genes interact with each target at that level of 
		# suppor or higher. Essentially, we are looking for regulatory modules which:
		# 1.	regulate the same target
		# 2.	invovle "real" interactions (given the high support)
		# 3.	are exclusive to GTEx or TCGA, indicating functional significance.
		#
		# ARGUMENTS:
		# * L:		A list formatted as the variables  L_gtex of L_tcga generated by method 
		#			computeDiffinteractions()
		# * min:	Minimum support threshold. Only consider interactions with support in at least
		#			"min" networks.
		#
		# RETURN VALUE
		# A sorted table T with one entry per target gene G involved in an exclusive interaction (R,G). 
		# The value T[G] is the number of hub genes R for which (R,G) has support at least "min".
		# Only genes G for which T[G] > 0 are reported. The entries in the table are sorted from highest
		# to lowest value.
		tabulateCommonRegulators <- function(L, min = 15){
			LL = unlist(sapply(L, function(x){
								if (x[1] < min)
									return(NA)
								else
									return(names(x[x >= min]))}))
			LL = LL[!is.na(LL)]
			LL = sort(table(LL), decreasing = TRUE)
		}
		
		load(params[[1]])
		min = params[[2]]
		iterations <- length(ls(pattern="L_random_"))
		par(mfrow=c(1,2))
		# Generate the GTEx plot
		title = "Number of regulators for high-support, GTEx-exclusive interactions"
		plot(density(tabulateCommonRegulators(L_actual[[1]], min)), col="red", 
				xlab = paste("Number of regulators, min support =", min), ylab = "density", 
				ylim = c(0, 0.8), main=title)
		for (i in 1:iterations){
			L = get(paste("L_random_", i, sep = ""))
			t = tabulateCommonRegulators(L[[1]], min)
			if (length(t) > 1)
				lines(density(t))
		}
		legend('topright', legend = c("GTEx", "Random"), col=c("red", "black"), lty=1)
		
		# Generate the TCGA plot
		title = "Number of regulators for high-support, TCGA-exclusive interactions"
		plot(density(tabulateCommonRegulators(L_actual[[2]], min)), col="red", 
				xlab = paste("Number of regulators, min support =", min), ylab = "density", 
				ylim = c(0, 0.8), main=title)
		for (i in 1:iterations){
			L = get(paste("L_random_", i, sep = ""))
			t = tabulateCommonRegulators(L[[2]], min)
			if (length(t) > 1)
				lines(density(t))
		}
		legend('topright', legend = c("TCGA", "Random"), col=c("red", "black"), lty=1)
	}
	
	
	# ***************************** which = cluster_vipers ********************************
	# Cluster tissues based on viper matrices
	# params[[1]]: a character vector of variable names of viper matrices
	# params[[2]]: a method when making one vector for one tissue (mean or median)
	# params[[3]]: a method for measuring correlation (euclidean, pearson, spearman, kendall)
	# e.g., params=list("varNamesVP", "median", "pearson")
	
	if (which == "cluster_vipers"){
		if(!is.null(params) && length(params) > 2) {
			
			### get variable names
			viperMatNames <- get(as.character(params[[1]]))
			
			### since all the viperMats have different set of hubs, get the commons
			common_hubs <- rownames(get(viperMatNames[1]))
			for(i in 2:length(viperMatNames)) {
				common_hubs <- intersect(common_hubs, rownames(get(viperMatNames[i])))
			}
			
			
			### make one vector for one tissue then combine all of them
			viperMat <- matrix(NA, length(common_hubs), length(viperMatNames))
			rownames(viperMat) <- common_hubs
			colnames(viperMat) <- substring(varNamesVP, 6)
			if(as.character(params[[2]]) == "mean") {
				for(i in 1:length(viperMatNames)) {
					viperMat[,i] <- apply(get(viperMatNames[i])[common_hubs,], 1, mean)
				}
			} else if(as.character(params[[2]]) == "median") {
				for(i in 1:length(viperMatNames)) {
					viperMat[,i] <- apply(get(viperMatNames[i])[common_hubs,], 1, median)
				}
			}
			else {
				stop("Error: incorrect params[[2]]")
			}
			
			### make a distance matrix
			if(as.character(params[[3]]) == "euclidean") {
				d <- dist(t(viperMat))
			} else if(as.character(params[[3]]) == "pearson") {
				d <- as.dist(1-cor(viperMat, method = "pearson"))
			} else if(as.character(params[[3]]) == "spearman") {
				d <- as.dist(1-cor(viperMat, method = "spearman"))
			} else if(as.character(params[[3]]) == "kendall") {
				d <- as.dist(1-cor(viperMat, method = "kendall"))
			} else {
				stop("Error: incorrect params[[3]]")
			}
			
			### make a plot
			par(mfrow=c(1,2), oma = c(0, 0, 2, 0))
			# hierarchical clustering
			plot(hclust(d), main="Hierarchical Clustering")
			# MDS
			fit <- cmdscale(d,eig=TRUE, k=2)
			x_c <- fit$points[,1]
			y_c <- fit$points[,2]
			plot(x_c, y_c, main="MDS",	type="n")
			text(x_c, y_c, labels = labels(d), cex=.7, pos=3)
			### common title
			fTitle <- paste(toupper(strsplit(viperMatNames[1], split = "_", fixed = TRUE)[[1]][2]),
					"clustering", params[[2]], params[[3]], sep = "_")
			mtext(fTitle, outer = TRUE, cex = 1.5)
			
		} else {
			stop("required params do not exist")
		}
	}
	
	
	# ***************************** which = viper_mds_per_tissue ********************************
	# Generates MDS plots per tissue
	# may find sub-types in each tissue
	# params[[1]]: a character vector of variable names of viper matrices
	# params[[2]]: a method for measuring correlation (euclidean, pearson, spearman, kendall)
	# e.g., params=list("varNamesVP", "euclidean")
	
	if (which == "viper_mds_per_tissue"){
		if(!is.null(params) && length(params) > 1) {
			
			### get variable names
			viperMatNames <- get(as.character(params[[1]]))
			
			### iteratively perform MDS plot for each tissue
			for(i in 1:length(viperMatNames)) {
				### get viper matrix of each tissue
				viperMat <- get(viperMatNames[i])
				
				### make a distance matrix
				if(as.character(params[[2]]) == "euclidean") {
					d <- dist(t(viperMat))
				} else if(as.character(params[[2]]) == "pearson") {
					d <- as.dist(1-cor(viperMat, method = "pearson"))
				} else if(as.character(params[[2]]) == "spearman") {
					d <- as.dist(1-cor(viperMat, method = "spearman"))
				} else if(as.character(params[[2]]) == "kendall") {
					d <- as.dist(1-cor(viperMat, method = "kendall"))
				} else {
					stop("Error: incorrect params[[2]]")
				}  
				
				### plot the MDS
				fit <- cmdscale(d,eig=TRUE, k=2)
				x_c <- fit$points[,1]
				y_c <- fit$points[,2]
				plot(x_c, y_c, main="MDS",	type="n")
				text(x_c, y_c, labels = labels(d), cex=.7, pos=3)
				### common title
				fTitle <- paste("Clustering", viperMatNames[i])
				mtext(fTitle, outer = TRUE, cex = 1.5)
			}
			
		} else {
			stop("required params do not exist")
		}
	}
	
	# ***************************** which = regulon_conservation ********************************
	# Generates images that assess how well the regulon of a query hub gene is conserved across
	# interactomes, utilizing the FET p-values in the object tfPairEnrich. A heatmap, a bar plot,
	# or a hierarchical clustering plot is generated, depending on the value of params[[2]].
	#
	# ARGUMENTS
	# * params[[1]]: the query hub gene, specified in any of the following ways:
	#		- gene symbol, entered as a character string.
	#		- entrez id, entered either as an integer or a character string.
	# * params[[2]]: one of four character strings: "heatmap" or "barplot" or "hierclust" or "mds"
	# params[[3]]: a positive integer, needed only when params[[2]] == "barplot". If not provided,
	#		a default value of 100 is used.
	#
	# RESULTS
	# If params[[2]] == "heatmap", it plots a square heatmap with on row and column per interactome,
	# showing how strong the regulon conservation is for every pair of interactomes. If params[[2]] ==
	# "barplot", it orders the FET pvalues for the query gene across all pairs of interactomes and 
	# selects the top params[[3]] most significant. For each of those it retrieves the relevant pair
	# A, B of interactomes and counts how many times each interactome is seen among the top params[[3]]
	# such pairs. It then plots a bar plot with these counts. If params[[2]] == "hierclust",
  # the distances between interactome pairs are represented as a dendrogram. Lastly, if params[[2]] ==
  # "mds", a MDS plot will be generated based on the distances between interactome pairs.
	if (which == "regulon_conservation") {
		
	  if(is.null(params) | (length(params) != 2 & length(params) !=3))
	    stop("Wrong number of parameters provided")
	  
		plotHeatMap <- function(gene){
			disp = 350
			mat = regulonConservationAcrossNets(gene, mode="matrix")
			for (i in 1:nrow(mat))
				for (j in 1:ncol(mat))
					if (mat[i,j] == -Inf){
						mat[i,j] = 0
					}else{
						#mat[i,j] = mat[i,j] + disp
						mat[i,j] = (mat[i,j] + disp)*log(mat[i,j] + disp+1)
					}
			rd<-as.dist(mat)
			rc<-hclust(rd)
			title = paste("Clustering of nets according to regulon conservation for gene:", 
					entrezIDtoSymbol(as.entrezId(gene)))
			par(cex.main=.8)
			heatmap(mat, Rowv=as.dendrogram(rc), symm=TRUE, main=title)
		}
		
		plotBarPlot <- function(gene, top = 100){
			par(mar=c(7,5,4,2)+.1)
			t=names(sort(getFets(gene, TRUE, TRUE))[1:top])
			t1 = sapply(t, function(x){return(strsplit(x, "%")[[1]][1])})
			t2 = sapply(t, function(x){return(strsplit(x, "%")[[1]][2])})
			t = sort(table(c(t1, t2)), decreasing=TRUE)
			title = paste("Network counts in top", top, "network pairs according to regulon conservation for gene:", gene)
			ylab = "Number of pairs appearing"
			mp = barplot(t, xaxt="n", main = title, ylab = ylab, cex.main = 0.8)
			axis(1, at=mp, labels=names(t), las = 2, cex.axis=0.8)
		}
		
		plotHclust <- function(gene) {
		  disp = 350
		  mat = regulonConservationAcrossNets(gene, mode="matrix")
		  for (i in 1:nrow(mat))
		    for (j in 1:ncol(mat))
		      if (mat[i,j] == -Inf){
		        mat[i,j] = 0
		      }else{
		        mat[i,j] = (mat[i,j] + disp)*log(mat[i,j] + disp+1)
		      }
		  rd<-as.dist(mat)
		  plot(hclust(rd), main = paste("Hierarchical Clustering of nets according to regulon conservation for gene:", 
		                                entrezIDtoSymbol(as.entrezId(gene))),
		       xlab = "", sub = "")
		}
		
		plotMDS <- function(gene) {
		  disp = 350
		  mat = regulonConservationAcrossNets(gene, mode="matrix")
		  for (i in 1:nrow(mat))
		    for (j in 1:ncol(mat))
		      if (mat[i,j] == -Inf){
		        mat[i,j] = 0
		      }else{
		        mat[i,j] = (mat[i,j] + disp)*log(mat[i,j] + disp+1)
		      }
		  rd<-as.dist(mat)
		  
		  ### get MDS points
		  fit <- cmdscale(rd,eig=TRUE, k=2)
		  Dimension1 <- fit$points[,1]
		  Dimension2 <- fit$points[,2]
		  plot(Dimension1, Dimension2, col = "white",
		       main = paste("MDS plot of nets according to regulon conservation for gene:", 
		                    entrezIDtoSymbol(as.entrezId(gene))))
		  text(Dimension1, Dimension2, labels = labels(rd), cex=.7, pos=3)
		}
		
		gene = params[[1]]
		if (params[[2]] == "heatmap")
			plotHeatMap(gene)
		if (params[[2]] == "barplot")
			if (length(params)==2)
				plotBarPlot(gene)
			else
				plotBarPlot(gene, top = params[[3]])
		if (params[[2]] == "hierclust")
		  plotHclust(gene)
		if (params[[2]] == "mds")
		  plotMDS(gene)
	}  
	
	
	# ***************************** which = viper_heatmap ********************************
	# Make a hierarchically clustered heatmap with the viper matrices
	# may find sub-types in each tissue
	# Needs all_64_vipermats.rda
	# params[[1]]: a character vector of variable names of viper matrices
	# params[[2]]: a method for measuring correlation (euclidean, pearson, spearman, kendall)
	# e.g., params=list("varNamesVP", "euclidean")
	
	if (which == "viper_heatmap"){
		if(!is.null(params) && length(params) > 1) {
			
			### get variable names
			viperMatNames <- get(as.character(params[[1]]))
			
			### iteratively perform MDS plot for each tissue
			for(i in 1:length(viperMatNames)) {
				### get viper matrix of each tissue
				viperMat <- get(viperMatNames[i])
				
				### make a distance matrix
				if(as.character(params[[2]]) == "euclidean") {
					d <- dist(t(viperMat))
				} else if(as.character(params[[2]]) == "pearson") {
					d <- as.dist(1-cor(viperMat, method = "pearson"))
				} else if(as.character(params[[2]]) == "spearman") {
					d <- as.dist(1-cor(viperMat, method = "spearman"))
				} else if(as.character(params[[2]]) == "kendall") {
					d <- as.dist(1-cor(viperMat, method = "kendall"))
				} else {
					stop("Error: incorrect params[[2]]")
				}
				
				### heatmap
				heatmap(as.matrix(d), main = paste0("Heatmap_", viperMatNames[i]),
						Rowv = as.dendrogram(hclust(d, method = "average")),
						symm = T, labRow = F, labCol = F)
			}
			
		} else {
			stop("required params do not exist")
		}
	}
	
	
	# ***************************** which = global_chromosome_target_enrichment ***************
	# Draw boxplots for each of the TCGA interactomes. For each interactome, the boxplot plots
	# the values -log10(Pval(H)), for all hub genes H in the interactome, where Pval(H) is the 
	# is the p-value of the chi-square test assessing if the targets in the regulon of H are 
	# localized preferentially on a specific interactome.  
	#
	# The call requires that the object "chr_enrich_all62") has been loaded to memory
	
	if (which == "global_chromosome_target_enrichment"){
		if (!exists("chr_enrich_all62"))
			stop("Object \"chr_enrich_all62\" not loaded.")
		par(mar=c(8,4,4,5)+.1)
		boxplot(-log10(chr_enrich_all62), las=2, ylab="-log10 chi-square p-value", 
				main="Boxplot of hub chi-square p-values per tissue")
	}
	
	# ***************************** which = hub_FET_enrichment_boxplot ***************
	# Given a hub gene H, generate a boxplot for the choose(length(varNames), 2) 
	# log10(FET pvalues) in tfPairEnrich[[2]] corresponding to H. Prior to plotting
	# the log p-values are ordered in increasing order (i.e., smallest p-value first).
	# Then, for pairwise interactions involving two TCGA or two GTEx networks, we use
	# -log10(p-value), i.e., a positive value. This way, when displaying the box plot
	# we can easily separate the significance of conservation within collection (positive
	# bars) vs. conservation across collections (negative bars) and thus easily assess
	# visually the dominant conservation mode (within or across collections).
	#
	# ARGUMENTS:
	# * params[[1]]:	The query hub gene, provided as gene symbol (string) or a entrez id
	#		(string or integer).
	# * params[[2]]:	Optional string argument. If present, it must assume the value of
	#		"GTEX" or "TCGA". In that case, a second barplot is plotted below the one 
	#		described above. This second barplot is identical except that log(pvalues) for 
	#		pairwise comparisons involving two GTEx inteactomes (if params[[2]] == "GTEX")
	#		or two TCGA interactomes (if params[[2]] == "TCGA") are zeroed out. This way
	#		for hubs where the conservation mode is dominantly within collections, we can 
	#		easily see which collection is the one involcing the stronger conservation. E.g.,
	#		in the case of TCGA dominant, the leading porton of the bottom boxplot will be 
	#		missing a lot of the bars that are present in the top boxplot.
	#
	# The call requires that the object "chr_enrich_all62") has been loaded to memory
	if (which == "hub_FET_enrichment_boxplot"){
		if (is.null(params) || is.na(params) || length(params) < 1 || length(params) > 2)
			stop("Wrong number of parameters")
		
		# Generate variables needed by the code that follows; do only once.
		if (!exists("reg_exclusivity_scores"))
			oneOffs(which = "generate_exclusivity_scores", params=list("FET_SUM"))
		
		hub = params[[1]]
		exclude = ifelse(length(params) == 2, params[[2]], NA)
		
		hub = as.entrezId(hub)
		if (!(hub %in% names(reg_exclusivity_scores)))
			return("Query gene is not a hub gene.")
		if (!is.na(exclude))
			par(mfrow=c(2, 1))

		minp = min(reg_fet_cons_mat[reg_fet_cons_mat != -Inf])
		fets = reg_fet_cons_mat[hub,]
		fets = -replace(fets, fets==-Inf, minp)
		fets = sort(fets, decreasing=TRUE)
		for (i in 1:length(fets))
			if (interactome_pair_map[names(fets)[i]] == "BOTH")
				fets[i] = -fets[i]
		max_y = max(fets)
		min_y = min(fets)
		barplot(fets, main = paste("FET p-value conservation plot for gene:",  entrezIDtoSymbol(hub)), 
				xaxt = "n", ylim=c(min_y, max_y), xlab = "Red lines every 250 pairwise interactome comparisons", 
				ylab = "-/+log10(pval) - same/different collection", cex.lab=0.8)
		for(l in seq(250, 1800, 250))
			abline(v = l, col = "red")
		if (!is.na(exclude)){
			for (i in 1:length(fets))
				if (interactome_pair_map[names(fets)[i]] == exclude)
					fets[i] = 0
			barplot(fets, main = paste("Excluding FET p-values where both interactomes are from ", exclude),
					xaxt = "n", ylim=c(min_y, max_y), ylab = "-/+log10(pval) - same/different collection", cex.lab=0.8)
			for(l in seq(250, 1800, 250))
				abline(v = l, col = "red")
		}
	}
  
	if(save)
		dev.off()
}

### Calculates FET score based on tfPairEnrich
calFet <- function(exclude = NULL) {
	if(exists("tfPairEnrich")) {
		
		n <- nrow(tfPairEnrich[[1]])
		fet_all <- matrix(0, n, n)
		fet_tf <- matrix(0, n, n)
		fet_tfcotf <- matrix(0, n, n)
		fet_signal <- matrix(0, n, n)
		
		rownames(fet_all) <- rownames(tfPairEnrich[[1]])
		colnames(fet_all) <- colnames(tfPairEnrich[[1]])
		rownames(fet_tf) <- rownames(tfPairEnrich[[1]])
		colnames(fet_tf) <- colnames(tfPairEnrich[[1]])
		rownames(fet_tfcotf) <- rownames(tfPairEnrich[[1]])
		colnames(fet_tfcotf) <- colnames(tfPairEnrich[[1]])
		rownames(fet_signal) <- rownames(tfPairEnrich[[1]])
		colnames(fet_signal) <- colnames(tfPairEnrich[[1]])
		
		### FET sum (1 = tf only, 2 = tf & cotf, 3 = signal only)
		sumFet <- function(enrich, digit, exclude = NULL) {
			keep = list(c("TF"), c("TF", "co-TF"), c("SIGN"))
			hubs = rownames(enrich)[geneTypeMap[rownames(enrich)] %in% keep[[digit]]]
			if (!is.null(exclude))
				hubs = setdiff(hubs, exclude)
			logps = enrich[hubs, "log10_pval"]
			minp <- min(logps[logps != -Inf])
			logps = replace(logps, logps==-Inf, minp)
			return(sum(logps))
		}
		
		if (!is.null(exclude))
			exclude = as.entrezId(exclude)
		for(i in 1:(n-1)) {
			for(j in (i+1):n) {
				fet_tf[j,i] <- sumFet(tfPairEnrich[[2]][tfPairEnrich[[1]][j,i]][[1]], 1, exclude)
				fet_tfcotf[j,i] <- sumFet(tfPairEnrich[[2]][tfPairEnrich[[1]][j,i]][[1]], 2, exclude)
				fet_signal[j,i] <- sumFet(tfPairEnrich[[2]][tfPairEnrich[[1]][j,i]][[1]], 3, exclude)
				fet_all[j,i] <- fet_tfcotf[j,i] + fet_signal[j,i]
			}
		}
		
		fet <- list(fet_all, fet_tf, fet_tfcotf, fet_signal)
		
		return(fet)
	} else {
		writeLines("tfPairEnrich does not exist.")
	}
}


### Calculate Probability that a regulon pair share same target genes
#   if the number of total genes = l
#      the number of target genes in regulon A = m
#      the number of target genes in regulon B = n
#      the number of shared target genes between A and B = k
#
#   Prob = (C(l, k) x C (l-k, m-k) x C(l-m, n-k)) / (C(l, m) x C(l, n))
#
#   * C(a, b) = choose(a,b) = P(a, b) / P(b, b)
#   * C() = Combination, P() = Permutation
###
calNetProb <- function(exclude = NULL) {
  if(exists("tfPairProb")) {
    
    n <- nrow(tfPairProb[[1]])
    prob_all <- matrix(0, n, n)
    prob_tf <- matrix(0, n, n)
    prob_tfcotf <- matrix(0, n, n)
    prob_signal <- matrix(0, n, n)
    
    rownames(prob_all) <- rownames(tfPairProb[[1]])
    colnames(prob_all) <- colnames(tfPairProb[[1]])
    rownames(prob_tf) <- rownames(tfPairProb[[1]])
    colnames(prob_tf) <- colnames(tfPairProb[[1]])
    rownames(prob_tfcotf) <- rownames(tfPairProb[[1]])
    colnames(prob_tfcotf) <- colnames(tfPairProb[[1]])
    rownames(prob_signal) <- rownames(tfPairProb[[1]])
    colnames(prob_signal) <- colnames(tfPairProb[[1]])
    
    ### prob sum (1 = tf only, 2 = tf & cotf, 3 = signal only)
    sumProb <- function(enrich, digit, exclude = NULL) {
      keep = list(c("TF"), c("TF", "co-TF"), c("SIGN"))
      hubs = rownames(enrich)[geneTypeMap[rownames(enrich)] %in% keep[[digit]]]
      if (!is.null(exclude))
        hubs = setdiff(hubs, exclude)
      logps = enrich[hubs, "log10_prob"]
      minp <- min(logps[logps != -Inf])
      logps = replace(logps, logps==-Inf, minp)
      return(sum(logps))
    }
    
    if (!is.null(exclude))
      exclude = as.entrezId(exclude)
    for(i in 1:(n-1)) {
      for(j in (i+1):n) {
        prob_tf[j,i] <- sumProb(tfPairProb[[2]][tfPairProb[[1]][j,i]][[1]], 1, exclude)
        prob_tfcotf[j,i] <- sumProb(tfPairProb[[2]][tfPairProb[[1]][j,i]][[1]], 2, exclude)
        prob_signal[j,i] <- sumProb(tfPairProb[[2]][tfPairProb[[1]][j,i]][[1]], 3, exclude)
        prob_all[j,i] <- prob_tfcotf[j,i] + prob_signal[j,i]
      }
    }
    
    prob <- list(prob_all, prob_tf, prob_tfcotf, prob_signal)
    
    return(prob)
  } else {
    writeLines("tfPairProb does not exist.")
  }
}


# *****************************************************************************
# Assesses the enrichment in co-modulating genes between each of the the 
# top-modulated TFs in a tumor and all other TFs in that tumor. It does so
# for all tumors.
#
# Arguments:
# * numTFs:		the number of top-modulated TFs to use in each tumor. If the tunor
#				does not have that many TFs being modulated, then all the available
#				TFs are used.
# * cutoff:		significance cutoff, to decide the minimum enrichment required 
#				for reporting co-modulated reguators. Specifically, for a given 
#				top-modulated TF A, we only return TFs B such that the modulators
#				of B are enriched in modulators of A at a (corrected) p-value 
#				(per Fisher's exact test) < cutoff. Correction here is at the level
#				of individual TFs: the FET p-value for the comparison between A and 
#				all possible Bs is multiplied by the number of all possible TFs B 
#				in a tumor. Howver, we do not correct for the fact that we are testing
#				for multiple TFs A and mulitple tumors. In that case we should further
#				mutiply p-values by numTFs * 20 (20 is the total number of tumors).
# * removeNulls: If TRUE then we only report TFs A for which there is are least one TF
#				B such that the modulators of A are signficantly enriched in the
#				modulators of B.  Otherwise a TF A will be reported even if there is no
#				enrichment.
# *varNamesDG:	vector containing the names of the DIGGIT data objects.
#
# Returns a list with one entry for  each tumor in 'varNamesDG'. Each entry is itself 
# a list L with 'numTFs' named members, each member corresponding to a TF A among the 
# numTFs most modulated TFs in the tumor. L[[A]] contains a vector V with as many elements 
# as the number of TFs B whose modulators are found to be enriched in modulators  
# of A. V contains named entries, with names(V) being the gene ID of such TFs B 
# and V[i] being the corrected p-value of Fisher's exact test that compares
# the size of the intersection of the gene modulator of A and names(V[i]). Ae mentioned
# above, if removeNulls == TRUE and length(v) == 0, then L[[A]] is purged from the 
# final results.
# *****************************************************************************
globalCoModulatedRegulators <- function(numTFs, cutOff = 0.01, removeNulls = TRUE, varDG = varNamesDG){
	
	results = vector("list", length(varDG))
	names(results) = varDG
	gCounts = globalCountOfModulators(varDG, FALSE)
	
	for (index in 1:length(varDG)){
		# logLines(paste("Processing tunor -->", varDG[index]))
		varDG = get(varDG[index])
		num = min(numTFs, length(varDG[[1]][,1]))
		tumorRes = vector("list", num)
		names(tumorRes) = varDG[[1]][1:num, 1]
		for (i in 1:num){
			tfGeneID = varDG[[1]][i,1]
			v = coModulatedRegulators(tfGeneID, varDG, TRUE, gCounts[index])
			if (!is.null(v)){
				v = v[v < cutOff]
				if (length(v) > 0)
					tumorRes[[i]] = v
			}
		}
		results[[index]] = tumorRes
	}
	if(removeNulls)
		results = sapply(results, function(e){return(e[!sapply(e, is.null)])})
	return(results)
}


# *****************************************************************************
# For every tumor in varDG, get all TFs genes that are both highly modulated and
# are themselves highly active modulators of other TFs 
# 
# Arguments:
# * top:	for each tunor the function constructs 2 data sets: the set A of the 
#			'top' TFs modulated by the highest number of DIGGIT events; and the set B
#			of 'top' highest modulator genes, i.e., genes modulating the highest 
#			number of TFs. The result object contians the intersection of A and B.
# * varDG:	character vector containing the variable names of DIGGIT data objects
#			corresponding to the tumors to be analyzed, in the form "gbmDG".
#
# Returns a list with one element for every tumor in varDG (list members are
# named, in the format 'gbmDG'). Each list member is a characted vector 
# containing the symbol names of the genes in the intersection of A and B.
# *****************************************************************************
getTopRegsAndMods <- function (top = 100, varDG = varNamesDG){
	# get a table listing for each modulator gene (rows) the number of TFs 
	# modulated by the gene in each tumor (columns)
	z = globalModulatorMatrix(TRUE, varNamesDG = varDG)
	L = list()
	for (i in 1:length(varDG)) 
		L[[i]] = z[,i]
	#y = sapply(L, function(e){return(strtoi(names(sort(e, decreasing = TRUE)[1:top])))})
	#y = sort(table(y), decreasing = TRUE)
	#entrezIDtoSymbol(names(y[1:20]))
	
	# For reach tumor report the genes that are both among the 'top' most modulated 
	# regulators and the TOP most active modulators
	res = sapply(seq(length(varDG)), function(i){
				x = intersect(strtoi(names(sort(L[[i]], decreasing = TRUE)[1:top])),
						get(varDG[i])[[1]][1:top, 1]); 
				if(length(x) == 0) 
					return(NULL) 
				return(entrezIDtoSymbol(x))})
	names(res) = varDG
	return(res)
}


# *****************************************************************************
# Counts the number of genes modulated by each modulator gene (with or without
# duplicates).
#
# Arguments:
# * varDG:	A DIGGIT data object corresponding to a tumor; or the character 
#			string name of the tumor, in the format 'brcaDG'
# * unique:	This controls what happens if multiple alterations in a modulator 
#			are identified to be associated with TF A. If unique = TRUE, all of
#			these alterations count as only 1 event. Otherwise, each alteration
#			separetely adds 1 to the total count for B
#
# Returns a named vector V where V[i] is the number of genes modulated by the 
# modulator with gene ID names(V)[i]. Vector entries are ordered from highest
# to smallest.
# *****************************************************************************
modulatorCounts <- function(varDG, unique = FALSE){
	if (typeof(varDG) == "character")
		varDG = get(varDG)
	resMat = matrix(0, max_geneId, 2)
	for (i in 1:length(varDG[[2]])){
		if (unique)
			tmp = unique(varDG[[2]][[i]][,1])
		else
			tmp = varDG[[2]][[i]][,1]
		for (j in 1:length(tmp)){
			resMat[tmp[j], 2] = resMat[tmp[j], 2] + 1
			resMat[tmp[j], 1] = tmp[j]
		}
	}
	resMat = resMat[resMat[,2] > 0,]
	resMat = resMat[order(resMat[,2], decreasing = TRUE),]
	results = resMat[,2]
	names(results) = resMat[,1]
	return(results)
}


# *****************************************************************************
# Test the enrichment in the intersection of the set of top-modulated regulators 
# and the set of highly active regulators, within a query tumor. Generally, one 
# would expect that highly modulated regulators are also also active in a large 
# proportion of the samples (otherwise, why bother modulating them?).
#
# Arguments:
# * varDG:	A DIGGIT data object corresponding to the query tumor; or the character 
#			string name of the tumor, in the format 'brcaDG'
# * vipers:	A table containing the VIPER data for the query tumor, compiled bu
#			reading one of Federicos VIPER files, e.g., viperValues-brca.dat.
#			The table has dimensions M x 2 where M is the number of lines in the
#			VIPER file. Eacf table row corresponds to a line from that file and
#			contains the first and third column from the line, i.e., the gene ID
#			of the VIPER regulator and the corresponding activity value. We don't
#			keep track of the 2nd column in the data file, which contains the ID
#			of the sample whose regulator activity is recorded in that line.
# * min:	Cutoff for filering our low-activity regulators, see below.
# * mode:	A string specifying how the 'min' cutoff is applied. If it has a value
#			of "both" then we will filter out all regulators with activity less than
#			abs(min). If it has a value of "pos" then we will filter out all regulators
#			with activity less than 'min'. Otherwise we will filter out all regulators
#			with activity > -min. 
# * top:	The number of most modulated regulators to use in the calculations.
# * random:	If TRUE then instead of intersecting the top modulated regulators with the 
#			most active regulators we instead intersect the top modulated regulators
#			with a random set of regulators. This argument was introduced in order
#			to generate empirical nulls: after computing the function with a given
#			set of parameters and random = FALSE, we can rerun it a number of times
#			with the same parameters and random = TRUE to get a bunch of nulls.
#
# The function works as follows:
# 1. Remove from 'vipers' all rows where the regulator gene ID is not among
#    the gene IDs found in varDG[[1]][,1], i.e., remove all regulator that are not
#	 found to be signficantly modulated by DIGGIT.
# 2. Filter out from 'vipers' all rows corresponding to regulators with
#    low activity, using the arguments 'min' and 'mode'.
# 3. Tabulate the remaining rows into a table tmp on a per-regulator basis, by:
#    * Counting the number of rows for each regulator (i.e. the number of samples
#      where a regulator is highly active).
#    * Sorting regultors in decreasing order, according to their count.
# 4. If random=TRUE, replace tmp with a same-size list of randomly chosen regulators
#    from varDG[[1]][,1].
# 5. Let A be the set of top modulated regulators, i.e., the set varDG[[1]][1:top,1]
# 6. For 1 < i < length(tmp), let B_i be the set of the first i entries in tmp. I.e.,
#    the set of the i most active regulators.
# 7. Let count_i be the intersection of A and B_i. And let pval_i be the p-value of 
#	 the Fisher's exact test forthe intersection of A and B_i.
#
#
# The function returns a data frame with 'length(tmp)' rows where the i-th row is 
# the pair (-log(pval_i), count_i).
# *****************************************************************************
incrementalFETs <- function (varDG, vipers, min = 3, mode = "both", top = 100, random = FALSE) {
	if (typeof(varDG) == "character")
		varDG = get(varDG)
	# Retain only the vipers involving regulators modulated by DIGGIT events
	tmp = vipers[vipers[,1] %in% varDG[[1]][,1],]
	# Filter out vipers where regulator activity is "low" and tabulate the results,
	# counting for each regulator the number of samples where its activity deviates 
	# signficantly from its mean (at least 'min' SDs).
	if (mode == "both")
		tmp = sort(table((tmp[abs(tmp[,2]) > min, 1])), decreasing = TRUE)
	else if (mode == "pos")
		tmp = sort(table((tmp[tmp[,2] > min, 1])), decreasing = TRUE)
	else 
		tmp = sort(table((tmp[tmp[,2] < -min, 1])), decreasing = TRUE)
	if (random){
		L = length(tmp)
		tmp = unique(vipers[vipers[,1] %in% varDG[[1]][,1],1])
		tmp = tmp[sample(seq(1:length(tmp)), L)]
		names(tmp) = tmp
	}
	pvals = vector(length = length(tmp), mode = "numeric")
	counts = vector(length = length(tmp), mode = "integer")
	
	N = length(varDG[[1]][,1])
	for (i in 1:length(tmp)){
		n_reg = min(top, length(varDG[[1]][,1]))
		counts[i] = length(intersect(strtoi(names(tmp))[1:i], varDG[[1]][1:n_reg,1]))
		dgOnly = length(setdiff(varDG[[1]][1:n_reg,1], strtoi(names(tmp))[1:i]))
		vpOnly = length(setdiff(strtoi(names(tmp))[1:i], varDG[[1]][1:n_reg,1]))
		none = N - dgOnly - vpOnly + counts[i]
		pvals[i] = -log(fisher.test(matrix(c(counts[i], vpOnly, dgOnly, none), 2, 2), alternative = "greater")$p.value)
	}
	return(data.frame(pvals, counts))
}


# *****************************************************************************
# Save to an RDA file the variables storing the primary and derived data from
# the ARACNe and DIGGIT files.
#
# Arguments:
# * mode:	Specify if saving should be done for ARACNe (mode == "ARACNe") or
#			DIGGIT (mode == "DIGGIT") files.
# * fNames:	File name used for the RDA file where variables will be stored.  
#
# *****************************************************************************
saveToRDA <- function (mode = "ARACNe", fName = NULL){
	
	if (mode == "ARACNe"){
		vars = c(varNames, "varNames", "pairWise", "netSizes", "tfPairEnrich", "tfPairProb", "tfNetEnrich", "README")
		if (is.null(fName))
			fName = "aracne.rda"
	}
	else if (mode == "DIGGIT"){
		vars = c(varNamesDG, "varNamesDG", "eventsPerTumor", "modsPerTumor", "altsPerTumor", "regsPerTumor", 
				"top100FET", "globModMatrix", "READMEDG")
		if (is.null(fName))
			fName = "tcga_diggit.rda"
	}
	else if (mode == "VIPER"){
		vars = varNamesVP
		if (is.null(fName))
			fName = "tcga_viper.rda"
	}
	
	save(list = vars, file = fName)
}


# *****************************************************************************
# For each tumor count (1) the number of hub genes in the corresponsing ARACNE 
# network, (2) the number of interactions in the same network, (3) the number of 
# regulators reported in the tumor's VIPER files, (4) the number of DIGGIT 
# associations, (5) the number of regulators modulated by DIGGIT events
#
# Arguments:
# * save:	Indicates if the results table should be stored to file (save = TRUE)
#			or not (save = FALSE).
# * fNames:	File name used for saving the results table, if save == TRUE.  
#
# Returns a table with 20 rows (one for each tumor) and 5 columns, each column
# containing the counts described above.
# *****************************************************************************
getSummaryCounts <- function(save = FALSE, fName = NULL){
	# Get the number of hubs and the number of interactions in each ARACNe network
	numRegsAR = sapply(varNames, function(e){return(length(get(e)[[1]][,1]))})
	numIntsAR = sapply(varNames, function(e){return(sum(get(e)[[1]][,3]))})
	
	# Get the number of VIPER regulators in each tumor and fix the columns names
	numRegsVP = sapply(varNamesVP, function(e){return(length(unique(get(e)[[1]][,1])))})
	names(numRegsVP) = gsub("VP", "", names(numRegsVP)) 
	
	# In each tumor get the total number of DIGGIT events and the number of distinct
	# regulators modulated by them (and fix columns names)
	numEventsDG = sapply(varNamesDG, function(e){return(sum(get(e)[[1]][,3]))})
	numRegsDG = sapply(varNamesDG, function(e){return(length(get(e)[[1]][,1]))})
	names(numEventsDG) = gsub("DG", "", names(numEventsDG))
	names(numRegsDG) = gsub("DG", "", names(numRegsDG))
	
	# Harmonize order of column names before merging.
	names = sort(names(numRegsAR))
	numRegsAR = numRegsAR[names]
	numIntsAR = numIntsAR[names]
	numRegsVP = numRegsVP[names]
	numEventsDG = numEventsDG[names]
	numRegsDG = numRegsDG[names]
	res = cbind(numRegsAR, numIntsAR, numRegsVP, numEventsDG, numRegsDG)
	if(save){
		if (is.null(fName))
			fName = "summary_counts_cptac.xlsx"
		new_names = c("ARACNe Hubs", "ARACNe Interactions", "VIPER Regulators", "DIGGIT Events", "DIGGIT Regulators")
		old_names = colnames(res)
		colnames(res) = new_names
		write.xlsx(res, file = fName, col.names = TRUE)
		colnames(res) = old_names
	}
	return(res)	
}


# *****************************************************************************
# Find the most frequently modulated regulators across all tumors. Specifically,
# for each tumor, identify the 'top' regulators with the most DIGGIT events. 
# Then merge these lists and count, for each regulator, in how many tumors it
# appears in.
#
# Arguments:
# * top:	Number of top most modulated regulators to consider, in each tumor.
#			If the number N of modulated regulators in a tumor is less that 
#			'top', then use	N 
# * varDG:	Character vector containing the names of the DIGGIT objects 
#			corresponding to the tumors to use, in the form "gbmDG".
# * filterSNV: If TRUE then SNVs are exluded from the analysis.
#
# Returns a named vector with one entry for each regulator in the merged list.
# Entries are named with the gene IDs of the corresponding regulators. The entry
# corresponding to a regulator A contains the number of tumors in varDG where A 
# is found to be among the most modulated regulators. Vector entries are listed 
# in decreasing order of that number.
# *****************************************************************************
getTopGlobalRegulators <- function(top = 50, varDG = varNamesDG, filterSNV = FALSE){
	varDG_tmp = paste(varDG, "tmp", sep="_") 
	for(i in 1:length(varDG)){
		assign(varDG_tmp[i], get(varDG[i])[[1]])
		if (filterSNV){
			x = get(varDG[i])[[1]]
			y = get(varDG[i])[[2]]
			for (j in 1:nrow(x))
				x[j,3] = sum(y[[x[j,2]]][,2] != 1)			
			x = x[order(x[,3], decreasing = TRUE), ]
			assign(varDG_tmp[i], x)
		}
	}
	topTFs = matrix(NA, nrow=top, ncol=length(varDG))
	colnames(topTFs) = varDG
	for(i in 1:length(varDG)){
		N = min(top, length(get(varDG_tmp[i])[, 1]))
		topTFs[1:N, i] = get(varDG_tmp[i])[1:N, 1]
	}
	return(sort(table(as.vector(topTFs)), decreasing = TRUE))
}


# *****************************************************************************
# Given a list of regulators R this method goes over the DIGGIT results in each
# tumor found in the global variable varNames DG and identifies the set M 
# of all the modudulators thar modulate one or more genes in R. It returns a 
# matrix |M| x |R| where rows are named after the genes in M and columns after 
# the genes in R and the [m, r]-th cell contains # the number of tumors where 
# the m-th modulator is reported to regulate the r-th # regulator.
#
# Arguments:
# * regulators:	Vector containing the gene IDs of the query regulators.
# * sort:		If TRUE, the rows in the result matrix are ordered, first by
#				maximum number of tumors in a row and then, in case of ties, 
#				by maximum row sums. Otherwise there is no ordering.
# * geneSymbols:	If FALSE, the row names and column names of the results
#		table are the gene IDs of the relevant genes. If TRUE, then they are 
#		the gene symbols. *Attention*: using geneSymbols = TRUE can impose a
#		significant performance premium. If gene symbols are desired for the
#		top modulators/regulators, it is best to do the mapping outside the 
#		function.
#
# Returns an |M| x |R| table, as described above.
# *****************************************************************************
modsVsRegsMatrix <- function(regulators, sort = TRUE, geneSymbols = FALSE){
	L = list()
	for (i in 1:length(regulators))
		L[[i]] = table(unlist(tfDGDetails(regulators[i], filterNulls = TRUE, idsOnly = TRUE)))
	M = List()
	for (i in 1:length(L))
		M[[i]] = strtoi(names(L[[i]]))
	res = matrix(0, nrow = length(unique(unlist(M))), ncol = length(regulators))
	rownames(res) = unique(unlist(M))
	colnames(res) = regulators
	
	for (i in 1:length(L)){
		vec = L[[i]]
		for (j in 1:length(vec))
			res[names(vec)[j], i] = vec[j]
	}
	
	if (geneSymbols){
		rownames(res) = entrezIDtoSymbol(unique(unlist(M)))
		colnames(res) = entrezIDtoSymbol(regulators)
	}
	
	if (sort){
		rs = rowSums(res)
		rc = apply(res, 1, max)
		ordering = order(rc, rs, decreasing = TRUE)
		res = res[ordering, ]
	}
	return(res)
}


# *****************************************************************************
# Identify all tumors where a gene is found among the tumors' top regulators.
#
# Arguments:
# * geneID:	The gene ID of the query gene.
# * top:	The number of top regulators to inspect.
# * varNamesDG: Character vector containing the names of DIGGIT data objects 
#		corresponding to the tumors to be analyzed, in the form "gbmDG".
#
# Returns a named vector V where names(V) == varNamesDG and V[i] is TRUE if 
# the query gene is among the 'top' regulators in the corresponding tumor, i.e.,
# if the gene ID of the query gene is found in the following set:
# 			get(varNamesDG[i])[[1]][1:top,1]
# *****************************************************************************
getTumorsForRegulator <- function(geneID, top = 50, varNamesDG = varNamesDG){
	return(sapply(varNamesDG, 
					function(e){
						x = get(e)[[1]][1:top,1]
						if (length(which(x == geneID)) > 0)
							return(TRUE)
						else
							return(FALSE)
					}))
}


# *****************************************************************************
# A "catch-all" method, to include miscellaneous code that was used for various 
# computations but which was not deemed worthy of separate "function" status.
#
# Arguments:
# * which:		Specifies which sub-function to run.
# * params:		A list contains subfunction-specific arguments. This is up to 
#				each subfunction to interpret as needed.
# *****************************************************************************		
oneOffs<- function (which = "freq_mods", params=NULL){
	if (which == "freq_mods"){
		# Get all regulators that are found in the top 150 most modulated genes in one
		# or more tumors.
		x150 = getTopGlobalRegulators(150)
		# Get the modulator genes the modulate regulators found in the top-150 lists of 
		# at least 7 tumors
		x = modsVsRegsMatrix(strtoi(names(x150[x150>6])), sort=TRUE)
		# For each regulator, identify the 30 modulators reported to modulated that regulator in the 
		# most number of tumors. The tabulate the results and for each modulator report the number of
		# regulators in whose top-30 list they appear.
		t_30 = sort(table(sapply(seq(1, length(colnames(x))), 
								function(i){return(rownames(x[order(x[,i], decreasing = TRUE),][1:30,]))})), 
				decreasing = TRUE)
		# As above, but identify the 50 (instead of of 30) top modulators per regulator.
		t_50 = sort(table(sapply(seq(1, length(colnames(x))), 
								function(i){return(rownames(x[order(x[,i], decreasing = TRUE),][1:50,]))})), 
				decreasing = TRUE)
		# As above, but identify the 100 (instead of of 50) top modulators per regulator.
		t_100 = sort(table(sapply(seq(1, length(colnames(x))), 
								function(i){return(rownames(x[order(x[,i], decreasing = TRUE),][1:100,]))})), 
				decreasing = TRUE)
		
		# For t_30, t_50, t_100: create an excel spreadsheet with one worksheet for each of these
		# 3 variables. Each worksheet lists the top 30 modulators based on the number of regulators
		# they are found to frequently modulate. Each worksheet line contains 2 columns: the number of
		# regulators frequently modulated, followed by the modulators that do the modulaating in these
		# regulators.
		fileName = "Reccurent_mods.xlsx"
		r = t_30[1:30]
		L = max(r) - min(r) + 1
		tmp = matrix(nrow = L, ncol = 2)
		tmp[,1] = seq(min(r), max(r))
		for (i in min(r):max(r)){
			rx = r[r==i]
			if (length(rx) > 0)
				tmp[i-min(r)+1,2] = toString(entrezIDtoSymbol(names(rx)))
			else
				tmp[i-min(r)+1,2] = ""
		}
		write.xlsx(apply(tmp,2,rev), file=fileName, sheetName="t_30")
		r = t_50[1:30]
		L = max(r) - min(r) + 1
		tmp = matrix(nrow = L, ncol = 2)
		tmp[,1] = seq(min(r), max(r))
		for (i in min(r):max(r)){
			rx = r[r==i]
			if (length(rx) > 0)
				tmp[i-min(r)+1,2] = toString(entrezIDtoSymbol(names(rx)))
			else
				tmp[i-min(r)+1,2] = ""
		}
		write.xlsx(apply(tmp,2,rev), file=fileName, append=TRUE, sheetName = "t_50")
		r = t_100[1:30]
		L = max(r) - min(r) + 1
		tmp = matrix(nrow = L, ncol = 2)
		tmp[,1] = seq(min(r), max(r))
		for (i in min(r):max(r)){
			rx = r[r==i]
			if (length(rx) > 0)
				tmp[i-min(r)+1,2] = toString(entrezIDtoSymbol(names(rx)))
			else
				tmp[i-min(r)+1,2] = ""
		}
		write.xlsx(apply(tmp,2,rev), file=fileName, append=TRUE, sheetName = "t_100")
		
		# Finally, take the union of the 30 top modulators from each of the lists
		# t_30, t_50, t_100 and compile a matrix where rows correspond to these 
		# modulators and rows correspond to to the top regulators in names(x150[x150>6]).
		# The [m ,r]-th entry in that matrix will be the number of tunors where the m-th
		# regulators is found to regulate the r-th regulator. Write this matrix as a 
		# new worksheet in the spreadsheet used above.
		n = union(names(t_30[1:30]), union(names(t_50[1:30]), names(t_100[1:30])))
		y = x[ n,]
		colnames(y) = entrezIDtoSymbol(colnames(y))
		rownames(y) = entrezIDtoSymbol(rownames(y))
		ordering = order(apply(y,1,max), rowSums(y), decreasing = TRUE)
		y = y[ordering,]
		write.xlsx(y, file=fileName, append=TRUE, sheetName = "Details")
	}
	if (which == "enriched_regulons"){
		fileName = "enriched_regulons.xlsx"
		how_many = 20
		steps = c(10, 20, 30)
		res = matrix("", nrow = how_many, ncol = length(steps))
		colnames(res) = paste("N = ", steps)
		for (i in 1:length(steps)){
			t = sort(table(unlist(sapply(tfPairEnrich[[2]], function(e){
												return(e[1:steps[i], 1])
											}))), decreasing = TRUE)
			t = t[1:how_many]
			names(t) =  entrezIDtoSymbol(names(t))
			for (j in 1:how_many)
				res[j,i] = paste(names(t)[j], " (", t[j], ")", sep = "")
		}
		write.xlsx(res, file=fileName)
	}
	# ----- Identify regulators that are highly active across multiple tumors. -----
	# This code expects that the params argument to have the following entries:
	# * params[[1]]:	The minimum number of top-N tunor lists that a regulator needs to
	#					appear, in order to be reported.
	# * params[[2]]:	A character value (either "P", "N", or "B") prescribing how to order
	#					regulators accrording to their average activity.
	# * params[[3]]:	A function name, prescribing which statistic to use for computing the average
	#					regulator actvity. Possible values can be "mean" or "median"
	# The code works as follows: First we order regulators in order of their mean VIPER activity across 
	# all samples, based on the value of the parameter params[[2]]. Specifically:
	# * If params[[2]] == "P", then we order from most positive to most negative.
	# * If params[[2]] == "N", then we order from most negative to most positive.
	# * If params[[2]] == "B", then we consider absoulte values and order from higher to lower.
	# After ordering the reguators as described above we then take the N top regulators in each tumor 
	# (N = 50, 100, 150)  and look for repeated occurrences, i.e., regulators that appear in multiple top-N 
	# lists. Finally, we report regulators that are found in at least "min" top-N lists, where min is the
	# values of the argument params[[1]]
	if (which == "top_active_across_tumors"){
		# check parameter values
		if (!(params[[2]] %in% c("P", "N", "B")))
			stop("params[[2]] not set correclty")
		
		L = lapply(varNamesVP, function(name){
					net = get(name)
					inc = length(unique(net[[2]]))
					x= nrow(net[[1]])/inc
					FUN = params[[3]]
					y = sapply(seq(1,x), function(k){return(FUN(net[[1]][((k-1)*inc+1):(k*inc), 2]))})
					names(y) = unique(net[[1]][,1])
					if (params[[2]] == "N")
						y = -y
					if (params[[2]] == "B")
						y = abs(y)
					return(sort(y, decreasing=T))
				})
		names(L) = varNamesVP
		t_50 = sort(table(unlist(lapply(L, function(x, top=50){return(strtoi(names(x)[1:top]))}))), decreasing=T)
		t_100 = sort(table(unlist(lapply(L, function(x, top=100){return(strtoi(names(x)[1:top]))}))), decreasing=T)
		t_150 = sort(table(unlist(lapply(L, function(x, top=150){return(strtoi(names(x)[1:top]))}))), decreasing=T)
		min = params[[1]]
		max = length(varNamesVP)
		res = matrix("", nrow = length(max:min), ncol=4)
		colnames(res) = c("Num_of_Tumors", "N_equals_50", "N_equals_100", "N_equals_150")
		res[,1] = max:min
		for (i in max:min){
			t = names(t_50[t_50 == i])
			if (length(t) > 0)
				res[max-i+1, 2] = toString(entrezIDtoSymbol(t))
			t = names(t_100[t_100 == i])
			if (length(t) > 0)
				res[max-i+1, 3] = toString(entrezIDtoSymbol(t))
			t = names(t_150[t_150 == i])
			if (length(t) > 0)
				res[max-i+1, 4] = toString(entrezIDtoSymbol(t))
		}
		fName = "temp_vip.xlsx"
		write.xlsx(res, file = fName, row.names = FALSE)
	}
	if (which == "top_modulated_regulators"){
		noSnv = FALSE
		if (!is.null(params) && (length(params) > 1) && (params[[2]] == TRUE))
			noSnv = TRUE
		t_50 = getTopGlobalRegulators(top = 50, filterSNV = noSnv)
		t_50 = t_50[t_50 > params[[1]]]
		names(t_50) = entrezIDtoSymbol(names(t_50))
		t_100 = getTopGlobalRegulators(top =100, filterSNV = noSnv)
		t_100 = t_100[t_100 > params[[1]]]
		names(t_100) = entrezIDtoSymbol(names(t_100))
		t_150 = getTopGlobalRegulators(top =150, filterSNV = noSnv)
		t_150 = t_150[t_150 > params[[1]]]
		names(t_150) = entrezIDtoSymbol(names(t_150))
		minNum = min(c(t_50, t_100, t_150))
		maxNum = max(c(t_50, t_100, t_150))
		res = matrix(nrow = maxNum - minNum + 1, ncol=3)
		colnames(res) = c("N=50", "N=100", "N=150")
		rownames(res) = maxNum:minNum
		for (ind in c("50", "100", "150")){
			x = get(paste("t", ind, sep="_"))
			for (i in minNum:maxNum)
				res[toString(i), paste("N", ind, sep="=")] = paste(names(x[x == i]), collapse = ",")
		}
		write.xlsx(res, file="top_modulated_regulators.xlsx")
		return(res)
	}
	
	# ******************** which = unique_regulons  *****************************
	# Reports all gene hubs that appear in only one regulon. 
	# 
	# Returns a matrix with one row per such hub gene, listing the gene Entrez ID,
	# the gene symbol, the network its reqgulon appears in, the size of that regulon
	# and the gene description. If params is NOT null then params[[1]] is expected 
	# to be the name of an Excel file where to write out the results.
	if(which == "unique_regulons"){
		hubs_seen_at_least_twice = unique(unlist(sapply(tfPairEnrich[[2]], function(x){return(x[,1])})))
		all_hubs = unique(unlist(sapply(varNames, function(x){return(get(x)[[1]][,1])})))
		unique_hubs = setdiff(all_hubs, hubs_seen_at_least_twice)
		
		# Find which network each unique hub appears in
		t = sapply(unique_hubs, function(gid){
					for (net in varNames){
						hubs = get(net)[[1]][,1]
						if (gid %in% hubs)
							return(net)
					}
				})
		res = matrix(nrow = length(unique_hubs), ncol = 5)
		colnames(res) = c("Entrez_ID", "Gene_Symbol", "Network", "Regulon_Size", "Gene_Description")
		rownames(res) = unique_hubs
		res[,1] = unique_hubs
		res[,2] = entrezIDtoSymbol(unique_hubs)
		res[,3] = t
		res[,4] = sapply(1:length(unique_hubs), function(i){
					return(get(t[i])[[1]][toString(unique_hubs[i]), 3])
				})
		temp <- entrezIDtoDescription(unique_hubs)
		temp[sapply(temp, is.null)] <- "NA"
		res[,5] = unlist(temp)
		
		# order results in decreasing order of hub regulon size
		res = res[order(strtoi(res[, 4]), decreasing = TRUE), ]
		# Write results to file, if requested.
		if (!is.null(params)){
			f_name = params[[1]]
			write.xlsx(res, file=f_name, row.names=FALSE)
		}
		
		return(res)
	}
	
	# ******************** which = 2_6_only  *****************************
	# Reports all gene hubs that appear in only 2-6 networks
	#
	# Returns a matrix with one row per such hub gene, listing the gene Entrez ID,
	# the gene symbol, the network its reqgulon appears in, the size of that regulon
	# and the gene description.
	# params[[1]]: the number of interactions that should be used for filtering out 
	#				 false-positive regulons. If NULL, no filtering will be performed.
	# params[[2]]: the file name of the output Excel file
	
	if(which == "2_6_only"){
		all_hubs = unique(unlist(sapply(varNames, function(x){return(get(x)[[1]][,1])})))
		m <- data.frame(matrix(0, length(all_hubs), 8))
		colnames(m) <- c("Entrez_ID", "Gene_Symbol", "Network_Names", "Regulon_Size", "Median_FET", "Median_Prob", "Gene_Description", "Count")
		rownames(m) <- all_hubs
		m$Entrez_ID <- all_hubs
		m$Gene_Symbol <- entrezIDtoSymbol(all_hubs)
		temp <- entrezIDtoDescription(all_hubs)
		temp[sapply(temp, is.null)] <- "NA"
		m$Gene_Description <- unlist(temp)	
		
		if(is.null(params) || is.null(params[[1]])) {
			N <- 0
		} else {
			N <- as.numeric(params[[1]])
		}
		
		### Now it's time to check how many regulons in all interactomes
		for(net in varNames) {
			temp <- get(net)[[1]][,1]
			temp2 <- get(net)[[1]][,3]
			
			for(i in 1:length(temp)) {
				if(temp2[i] > N) {
					m[as.character(temp[i]),"Count"] <- m[as.character(temp[i]),"Count"]+1
					m[as.character(temp[i]), "Network_Names"] <- paste0(m[as.character(temp[i]), "Network_Names"],net,",")
					m[as.character(temp[i]), "Regulon_Size"] <- paste0(m[as.character(temp[i]), "Regulon_Size"],temp2[i],",")
				}
			}
		}
		
		for(i in 1:nrow(m)) {
			m[i,3] <- substr(m[i,3], 2, nchar(m[i,3])-1)
			m[i,4] <- substr(m[i,4], 2, nchar(m[i,4])-1)
		}
		
		### Extract hubs which have count >=2 and count <= 6
		two_Six <- m[intersect(which(m$Count >= 2), which(m$Count <= 6)),]
		
		if(nrow(two_Six) > 0) {
			### A function to find statistic value from tfPairEnrich[[2]] or tfPairProb[[2]]
			findValue <- function(pairProb, geneID) {
				rownames(pairProb) <- pairProb[,1]
				pairProb[which(pairProb[,2] == -Inf),2] <- min(pairProb[which(pairProb[,2] != -Inf),2])
				
				return(pairProb[geneID,2])
			}
			
			for(i in 1:nrow(two_Six)) {
				v <- strsplit(two_Six[i,3], split=",", fixed=TRUE)[[1]]
				
				f <- 0
				p <- 0
				cnt <- 0
				for(j in 1:(length(v)-1)) {
					for(k in (j+1):length(v)) {
						f[cnt+1] <- findValue(tfPairEnrich[[2]][tfPairEnrich[[1]][v[j],v[k]]][[1]], as.character(two_Six[i,1]))
						p[cnt+1] <- findValue(tfPairProb[[2]][tfPairProb[[1]][v[j],v[k]]][[1]], as.character(two_Six[i,1]))
						cnt <- cnt+1
					}
				}
				
				two_Six[i,5] <- median(f)
				two_Six[i,6] <- median(p)
			}
			
			### order results in increasing order of column 'Count'
			two_Six <- two_Six[order(strtoi(two_Six[,8]), decreasing=FALSE),]
			
			### Write results to file, if requested.
			if(!is.null(params) && length(params) > 1) {
				f_name = params[[2]]
				write.xlsx(two_Six, file=f_name, row.names=FALSE)
			}
			
			return(two_Six)
		} else {
			writeLines("There are no 2-6 regulons")
		}
	}
	
	
	# ******************** which = shared_hubs  *****************************
	# Reports all gene hubs that are shared in all networks
	#
	# Returns a list that contains information of shared regulons
	# Also creates pathway analysis result for shared regulons
	#
	# params[[1]]: the number of interactions that should be used for filtering out 
	#				       false-positive regulons. If NULL, no filtering will be performed.
	# params[[2]]: the file name of the output Excel file
	# params[[3]]: the number of regulons to be used for pathway analysis
	# params[[4]]: the file name of the pathway analysis result file
	
	if(which == "shared_hubs"){
		
		all_hubs = unique(unlist(sapply(varNames, function(x){return(get(x)[[1]][,1])})))
		m <- data.frame(matrix(0, length(all_hubs), 8))
		colnames(m) <- c("Entrez_ID", "Gene_Symbol", "Network_Names", "Regulon_Size", "Median_FET", "Median_Prob", "Gene_Description", "Count")
		rownames(m) <- all_hubs
		m$Entrez_ID <- all_hubs
		m$Gene_Symbol <- entrezIDtoSymbol(all_hubs)
		temp <- entrezIDtoDescription(all_hubs)
		temp[sapply(temp, is.null)] <- "NA"
		m$Gene_Description <- unlist(temp)	
		
		if(is.null(params) || is.null(params[[1]])) {
			N <- 0
		} else {
			N <- as.numeric(params[[1]])
		}
		
		### Now it's time to check how many hubs in all interactomes
		for(net in varNames) {
			temp <- get(net)[[1]][,1]
			temp2 <- get(net)[[1]][,3]
			
			for(i in 1:length(temp)) {
				if(temp2[i] > N) {
					m[as.character(temp[i]),"Count"] <- m[as.character(temp[i]),"Count"]+1
					m[as.character(temp[i]), "Network_Names"] <- paste0(m[as.character(temp[i]), "Network_Names"],net,",")
					m[as.character(temp[i]), "Regulon_Size"] <- paste0(m[as.character(temp[i]), "Regulon_Size"],temp2[i],",")
				}
			}
		}
		
		for(i in 1:nrow(m)) {
			m[i,3] <- substr(m[i,3], 2, nchar(m[i,3])-1)
			m[i,4] <- substr(m[i,4], 2, nchar(m[i,4])-1)
		}
		
		### Extract hubs which have count = total number of networks
		shared_hubs <- m[m$Count == length(varNames),]
		
		if(nrow(shared_hubs) > 0) {
			### A function to find statistic value from tfPairEnrich[[2]] or tfPairProb[[2]]
			findValue <- function(pairProb, geneID) {
				rownames(pairProb) <- pairProb[,1]
				pairProb[which(pairProb[,2] == -Inf),2] <- min(pairProb[which(pairProb[,2] != -Inf),2])
				
				return(pairProb[geneID,2])
			}
			
			for(i in 1:nrow(shared_hubs)) {
				v <- strsplit(shared_hubs[i,3], split=",", fixed=TRUE)[[1]]
				
				f <- 0
				p <- 0
				cnt <- 0
				for(j in 1:(length(v)-1)) {
					for(k in (j+1):length(v)) {
						f[cnt+1] <- findValue(tfPairEnrich[[2]][tfPairEnrich[[1]][v[j],v[k]]][[1]], as.character(shared_hubs[i,1]))
						p[cnt+1] <- findValue(tfPairProb[[2]][tfPairProb[[1]][v[j],v[k]]][[1]], as.character(shared_hubs[i,1]))
						cnt <- cnt+1
					}
				}
				
				shared_hubs[i,5] <- median(f)
				shared_hubs[i,6] <- median(p)
			}
			
			### Write results to file, if requested.
			if(!is.null(params) && length(params) > 1) {
				f_name = params[[2]]
				write.xlsx(shared_hubs, file=f_name, row.names=FALSE)
			}
		} else {
			writeLines("There are no such regulons")
		}
		
		### order results in increasing order of column 'Med_FET' and 'Med_Prob'
		shared_hubs <- shared_hubs[order(shared_hubs[,5], shared_hubs[,6]),]
		
		if(is.null(params) || is.null(params[[3]])) {
			N2 <- 0
		} else {
			N2 <- as.numeric(params[[3]])
		}
		
		### select the predefined number of regulons only
		shared_hubs <- shared_hubs[1:N2,]
		
		### get gene symbols without NA
		geneList <- shared_hubs$Gene_Symbol
		geneList <- geneList[which(!is.na(geneList))]
		
		### run pathway analysis
		if(!is.null(params) && length(params) > 3) {
			pathwayAnalysis_TB(geneList = geneList, title = as.character(params[[4]]), fName = as.character(params[[4]]))
		} else {
			writeLines("params[[4]] does not exist")
		}
		
		return(shared_hubs)
	}
	
	
	# ******************** which = exclusive_hubs  *****************************
	# Reports all gene hubs that are exclusive between GTEx and TCGA
	#
	# Returns a list that contains information of exclusive regulons
	# Also creates pathway analysis result for the regulons
	#
	# params[[1]]: the excel file path of shared regulons of GTEx
	# params[[2]]: the excel file path of shared regulons of TCGA
	# params[[3]]: the directory path of the pathway analysis result files
	
	if(which == "exclusive_hubs"){
		
		if(!is.null(params) && length(params) > 2) {
			
			shared_hubs_gtex <- read.xlsx2(as.character(params[[1]]), sheetIndex = 1, stringsAsFactors = FALSE)
			shared_hubs_tcga <- read.xlsx2(as.character(params[[2]]), sheetIndex = 1, stringsAsFactors = FALSE)
			
			rownames(shared_hubs_gtex) <- shared_hubs_gtex$Entrez_ID
			rownames(shared_hubs_tcga) <- shared_hubs_tcga$Entrez_ID
			
			gtex_only <- shared_hubs_gtex[setdiff(shared_hubs_gtex$Entrez_ID, shared_hubs_tcga$Entrez_ID),]
			tcga_only <- shared_hubs_tcga[setdiff(shared_hubs_tcga$Entrez_ID, shared_hubs_gtex$Entrez_ID),]
			
			geneList <- gtex_only$Gene_Symbol
			geneList <- geneList[which(!is.na(geneList))]
			
			pathwayAnalysis_TB(geneList = geneList, title = "GTEx_Exclusive", fName = paste0(params[[3]], "GTEx_Exclusive.pdf"))
			
			geneList <- tcga_only$Gene_Symbol
			geneList <- geneList[which(!is.na(geneList))]
			
			pathwayAnalysis_TB(geneList = geneList, title = "TCGA_Exclusive", fName = paste0(params[[3]], "TCGA_Exclusive.pdf"))
			
		} else {
			writeLines("required params do not exist")
		}
	}
	
	
	# ******************** which = totally_exclusive_hubs  *****************************
	# Reports all gene hubs that are totally exclusive between GTEx and TCGA
	# Loading All_Aracne.rda (64 networks = 36 GTExs + 28 TCGAs) is needed before this analysis 
	#
	# Returns a list that contains information of totally exclusive regulons
	# Also creates pathway analysis result for the regulons
	#
	# params[[1]]: the index in varNames that separates GTEx and TCGA (e.g., if GTEx is 1:N, enter N)
	# params[[2]]: the excel file path of shared regulons of GTEx
	# params[[3]]: the excel file path of shared regulons of TCGA
	# params[[4]]: the directory path of the pathway analysis result files
	
	if(which == "totally_exclusive_hubs"){
		
		if(!is.null(params) && length(params) > 3) {
			N <- as.numeric(params[[1]])
			
			shared_hubs_gtex <- read.xlsx2(as.character(params[[2]]), sheetIndex = 1, stringsAsFactors = FALSE)
			shared_hubs_tcga <- read.xlsx2(as.character(params[[3]]), sheetIndex = 1, stringsAsFactors = FALSE)
			
			rownames(shared_hubs_gtex) <- shared_hubs_gtex$Entrez_ID
			rownames(shared_hubs_tcga) <- shared_hubs_tcga$Entrez_ID
			
			gtex_hubs <- as.character(unique(unlist(sapply(varNames[1:N], function(x){return(get(x)[[1]][,1])}))))
			tcga_hubs <- as.character(unique(unlist(sapply(varNames[(N+1):length(varNames)], function(x){return(get(x)[[1]][,1])}))))
			
			gtex_only <- setdiff(shared_hubs_gtex$Entrez_ID, tcga_hubs)
			tcga_only <- setdiff(shared_hubs_tcga$Entrez_ID, gtex_hubs)
			
			
			geneList <- entrezIDtoSymbol(gtex_only)
			write.xlsx2(geneList, paste0(params[[4]], "GTEx_Totally_Exclusive.xlsx"), row.names=TRUE, col.names=FALSE)
			geneList <- geneList[which(!is.na(geneList))]
			pathwayAnalysis_TB(geneList = geneList, title = "GTEx_Totally_Exclusive", fName = paste0(params[[4]], "GTEx_Totally_Exclusive.pdf"))
			
			geneList <- entrezIDtoSymbol(tcga_only)
			write.xlsx2(geneList, paste0(params[[4]], "TCGA_Totally_Exclusive.xlsx"), row.names=TRUE, col.names=FALSE)
			geneList <- geneList[which(!is.na(geneList))]
			pathwayAnalysis_TB(geneList = geneList, title = "TCGA_Totally_Exclusive", fName = paste0(params[[4]], "TCGA_Totally_Exclusive.pdf"))
			
		} else {
			writeLines("required params do not exist")
		}
	}
	
	# ******************** which = TCGA_spreadsheet  *****************************
	# Load TCGA sample matrices, combine them, and save them as an Excel file
	# cancer-expmat.rda is needed. It contains all the TCGA expression information
	#
	# params[[1]]: the input file path - cancer-expmat.rda
	# params[[2]]: the output Excel file path
	
	if(which == "TCGA_spreadsheet"){
		if(!is.null(params) && length(params) > 1) {
			
			load(as.character(params[[1]]))
			
			tissues <- unique(samples[,1])
			
			wb <- createWorkbook()
			
			for(i in 1: length(tissues)) {
				sheet <- createSheet(wb, sheetName = paste("TCGA", toupper(tissues[i]), sep="_"))
				
				d <- samples[which(samples[,1] == tissues[i]),]
				d <- as.data.frame(d[,-1])
				
				addDataFrame(d, sheet, row.names = FALSE)
			}
			
			saveWorkbook(wb, as.character(params[[2]]))
		} else {
			writeLines("required params do not exist")
		}  
	}
	
	# ******************** which = compare_two_networks  *****************************
	# get common interactions and their p-values from the two networks
	# then draw cor & line graph of their p-values for comparison
	#
	# params[[1]]: network1 name. i.g., "Spleen"
	# params[[2]]: network2 name. i.g., "spleen2"
	# params[[3]]: the cor & line graph path to be stored
	
	if(which == "compare_two_networks"){
		if(!is.null(params) && length(params) == 3) {
			
			### load library
			if(!require(ggplot2)) {
				install.packages("ggplot2")
				library(ggplot2)
			}
			if(!require(VennDiagram)) {
				install.packages("VennDiagram")
				library(VennDiagram)
			}
			if(!require(gridExtra)) {
				install.packages("gridExtra")
				library(gridExtra)
			}
			
			### get two networks
			n1 <- get(params[[1]])[[2]]
			n2 <- get(params[[2]])[[2]]
			
			### shared hubs between the two networks
			shared_hubs <- intersect(names(n1), names(n2))
			sharedHubIdx1 <- which(names(n1) %in% shared_hubs)
			sharedHubIdx2 <- which(names(n2) %in% shared_hubs)
			
			### shared targets between the two networks for each hubs
			sharedTargetIdx1 <- list()
			sharedTargetIdx2 <- list()
			for(i in 1:length(shared_hubs)) {
				temp <- intersect(abs(n1[[sharedHubIdx1[i]]][,1]), abs(n2[[sharedHubIdx2[i]]][,1]))
				
				sharedTargetIdx1[[i]] <- which(abs(n1[[sharedHubIdx1[i]]][,1]) %in% temp)
				sharedTargetIdx2[[i]] <- which(abs(n2[[sharedHubIdx2[i]]][,1]) %in% temp)
			}
			
			
			columnBind <- function(hub, a, b) {
				x <- abs(matrix(a, length(a)/3, 3))
				y <- abs(matrix(b, length(b)/3, 3))
				
				x <- matrix(x[order(x[,1]),], length(a)/3, 3)
				y <- matrix(y[order(y[,1]),], length(a)/3, 3)
				
				if(nrow(x) == 1) {
					z <- c(hub, x, y[-1])
				} else {
					z <- cbind(hub, x, y[,-1])
				}
				
				return(z)
			}
			
			tempCnt <- 1
			while(length(sharedTargetIdx1[[tempCnt]]) < 1) {
				tempCnt <- tempCnt + 1
			}
			intrcn <- columnBind(shared_hubs[tempCnt], n1[[sharedHubIdx1[tempCnt]]][sharedTargetIdx1[[tempCnt]],-c(3,4)], n2[[sharedHubIdx2[tempCnt]]][sharedTargetIdx2[[tempCnt]],-c(3,4)])
			for(i in 2:length(shared_hubs)) {
				if(length(sharedTargetIdx1[[i]]) > 0) {
					intrcn <- rbind(intrcn, columnBind(shared_hubs[i], n1[[sharedHubIdx1[i]]][sharedTargetIdx1[[i]],-c(3,4)], n2[[sharedHubIdx2[i]]][sharedTargetIdx2[[i]],-c(3,4)]))
				}
			}
			intrcn <- data.frame(intrcn)
			colnames(intrcn) <- c("Gene1", "Gene2", "MI1", "PV1", "MI2", "PV2")
			idx <- sapply(intrcn, is.factor)
			intrcn[idx] <- lapply(intrcn[idx], function(x) as.numeric(as.character(x)))
			
			### correlation plot
			ggplot(intrcn, aes(x=MI1, y=MI2)) +
					ggtitle("Cor of MIs between the two networks", subtitle=sprintf("P.Cor = %s", round(cor(intrcn$MI1, intrcn$MI2), 5))) +
					labs(x = params[[1]], y=params[[2]]) +
					theme_classic(base_size = 16) +
					geom_point() +
					geom_smooth(method = lm, color="red", se=FALSE)
			ggsave(paste0(params[[3]], params[[1]], "_", params[[2]], "_Cor_MI.png"), width = 10, height = 10)
			
			ggplot(intrcn, aes(x=PV1, y=PV2)) +
					ggtitle("Cor of PVs between the two networks", subtitle=sprintf("P.Cor = %s", round(cor(intrcn$PV1, intrcn$PV2), 5))) +
					labs(x = params[[1]], y=params[[2]]) +
					theme_classic(base_size = 16) +
					geom_point() +
					geom_smooth(method = lm, color="red", se=FALSE)
			ggsave(paste0(params[[3]], params[[1]], "_", params[[2]], "_Cor_PV.png"), width = 10, height = 10)
			
			### make data for point plot
			mi_data <- data.frame(rbind(cbind(Group="1", MI=intrcn[,3]), cbind(Group="2", MI=intrcn[,5])))
			mi_data <- cbind(Interaction=c(rownames(mi_data)[1:nrow(intrcn)], rownames(mi_data)[1:nrow(intrcn)]), mi_data)
			mi_data <- mi_data[order(as.numeric(as.character(mi_data$Interaction))),]
			
			pv_data <- data.frame(rbind(cbind(Group="1", PV=-log10(intrcn[,4])), cbind(Group="2", PV=-log10(intrcn[,6]))))
			pv_data <- cbind(Interaction=c(rownames(pv_data)[1:nrow(intrcn)], rownames(pv_data)[1:nrow(intrcn)]), pv_data)
			pv_data <- pv_data[order(as.numeric(as.character(pv_data$Interaction))),]
			
			### point graph
			ggplot(mi_data, aes(x=Interaction, y=MI, group=Group)) +
					ggtitle(paste0(params[[1]], "_", params[[2]], "_Point_Plot_MI")) +
					theme_classic(base_size = 16) +
					theme(axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank()) +
					labs(x = paste(nrow(intrcn), "Interactions"), y="Mutual-Informations") +
					geom_point(aes(color=Group)) +
					scale_color_manual(labels = c(params[[1]], params[[2]]), values = c("blue", "orange")) +
					scale_size_manual(values = c(2, 2))
			ggsave(paste0(params[[3]], params[[1]], "_", params[[2]], "_Point_MI.png"), width = 10, height = 10)
			
			ggplot(pv_data, aes(x=Interaction, y=PV, group=Group)) +
					ggtitle(paste0(params[[1]], "_", params[[2]], "_Point_Plot_PV")) +
					theme_classic(base_size = 16) +
					theme(axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank()) +
					labs(x = paste(nrow(intrcn), "Interactions"), y="P-values") +
					geom_point(aes(color=Group)) +
					scale_color_manual(labels = c(params[[1]], params[[2]]), values = c("blue", "orange")) +
					scale_size_manual(values = c(2, 2))
			ggsave(paste0(params[[3]], params[[1]], "_", params[[2]], "_Point_PV.png"), width = 10, height = 10)
			
			
			### hubs venn diagram
			### get hub and degree info
			net1 <- data.frame(get(params[[1]])[[1]])
			net2 <- data.frame(get(params[[2]])[[1]])
			
			### draw a venn diagram of hubs
			v <- venn.diagram(list(net1$X1, net2$X1), category.names = c(params[[1]], params[[2]]), cat.cex = 1.5, cex = 1.5, filename = NULL)
			
			### save the diaram as png
			png(paste0(params[[3]], params[[1]], "_", params[[2]], "_Venn_Diagram_Hubs.png"))
			grid.arrange(gTree(children=v), top="Shared Hubs", bottom="")
			dev.off()
			
			
			### hubs degree correlation plot
			hubs_degree <- data.frame(cbind(net1[shared_hubs,3], net2[shared_hubs,3]))
			rownames(hubs_degree) <- shared_hubs
			
			### correlation plot
			ggplot(hubs_degree, aes(x=X1, y=X2)) +
					ggtitle("Cor of degrees between the two networks", subtitle=sprintf("P.Cor = %s", round(cor(hubs_degree$X1, hubs_degree$X2), 5))) +
					labs(x = params[[1]], y=params[[2]]) +
					theme_classic(base_size = 16) +
					geom_point() +
					geom_smooth(method = lm, color="red", se=FALSE)
			ggsave(paste0(params[[3]], params[[1]], "_", params[[2]], "_Cor_Degree.png"), width = 10, height = 10)
			
		} else {
			writeLines("required params do not exist")
		}  
	}
	
	# ******************** which = relatively_exclusive_hubs  *****************************
	# Reports all gene hubs that are relatively exclusive between GTEx and TCGA
	# Loading All_62_Aracne.rda (62 networks = 36 GTExs + 26 TCGAs) is needed before this analysis
	#
	# params[[1]]: the index in varNames that separates GTEx and TCGA (e.g., if GTEx is 1:N, enter N)
	# params[[2]]: the number of interactions that should be used for filtering out 
	#				       false-positive regulons. If NULL, no filtering will be performed.
	# params[[3]]: the ratio of hub existence in GTEx (e.g., 20) [0-100]
	# params[[4]]: the ratio of hub existence in TCGA (e.g., 80) [0-100]
	#
	# For the ratio values, 80 means hubs appeared in 80% or more of the networks
	# And 20 means hubs appeared in 20% or less of the networks
	# ratio R >= 50: [R-100], ratio R < 50: [0-R]  
	#
	# params[[5]]: the output text file directory - contains result regulon list
	# params[[6]]: the output pathway file directory
	#
	# function usage example: oneOffs(which="relatively_exclusive_hubs", params=list(36, 10, 20, 80, "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/results/", "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/pathway/"))
	
	if(which == "relatively_exclusive_hubs"){
		if(!is.null(params) && length(params) == 6) {
			
			all_hubs <- unique(unlist(sapply(varNames, function(x){return(get(x)[[1]][,1])})))
			
			### GTEx
			gtex_m <- data.frame(matrix(0, length(all_hubs), 6))
			colnames(gtex_m) <- c("Entrez_ID", "Gene_Symbol", "Network_Names", "Regulon_Size", "Gene_Description", "Count")
			rownames(gtex_m) <- all_hubs
			gtex_m$Entrez_ID <- all_hubs
			gtex_m$Gene_Symbol <- entrezIDtoSymbol(all_hubs)
			temp <- entrezIDtoDescription(all_hubs)
			temp[sapply(temp, is.null)] <- "NA"
			gtex_m$Gene_Description <- unlist(temp)	
			
			### TCGA
			tcga_m <- gtex_m
			
			if(is.null(params[[2]])) {
				N <- 0
			} else {
				N <- as.numeric(params[[2]])
			}
			
			### GTEx
			### Now it's time to check how many hubs in all interactomes
			for(net in varNames[1:params[[1]]]) {
				temp <- get(net)[[1]][,1]
				temp2 <- get(net)[[1]][,3]
				
				for(i in 1:length(temp)) {
					if(temp2[i] > N) {
						gtex_m[as.character(temp[i]),"Count"] <- gtex_m[as.character(temp[i]),"Count"]+1
						gtex_m[as.character(temp[i]), "Network_Names"] <- paste0(gtex_m[as.character(temp[i]), "Network_Names"],net,",")
						gtex_m[as.character(temp[i]), "Regulon_Size"] <- paste0(gtex_m[as.character(temp[i]), "Regulon_Size"],temp2[i],",")
					}
				}
			}
			
			for(i in 1:nrow(gtex_m)) {
				if(gtex_m[i,3] == 0) {
					gtex_m[i,3] <- NA
					gtex_m[i,4] <- NA
				} else {
					gtex_m[i,3] <- substr(gtex_m[i,3], 2, nchar(gtex_m[i,3])-1)
					gtex_m[i,4] <- substr(gtex_m[i,4], 2, nchar(gtex_m[i,4])-1)
				}
			}
			
			gtex_m <- gtex_m[order(gtex_m$Count),]
			
			
			### TCGA
			### Now it's time to check how many regulons in all interactomes
			for(net in varNames[(params[[1]]+1):length(varNames)]) {
				temp <- get(net)[[1]][,1]
				temp2 <- get(net)[[1]][,3]
				
				for(i in 1:length(temp)) {
					if(temp2[i] > N) {
						tcga_m[as.character(temp[i]),"Count"] <- tcga_m[as.character(temp[i]),"Count"]+1
						tcga_m[as.character(temp[i]), "Network_Names"] <- paste0(tcga_m[as.character(temp[i]), "Network_Names"],net,",")
						tcga_m[as.character(temp[i]), "Regulon_Size"] <- paste0(tcga_m[as.character(temp[i]), "Regulon_Size"],temp2[i],",")
					}
				}
			}
			
			for(i in 1:nrow(tcga_m)) {
				if(tcga_m[i,3]== 0) {
					tcga_m[i,3] <- NA
					tcga_m[i,4] <- NA
				} else {
					tcga_m[i,3] <- substr(tcga_m[i,3], 2, nchar(tcga_m[i,3])-1)
					tcga_m[i,4] <- substr(tcga_m[i,4], 2, nchar(tcga_m[i,4])-1)
				}
			}
			
			tcga_m <- tcga_m[order(tcga_m$Count),]
			
			
			### Extract hubs which have count = pre-defined relative number of networks
			if(params[[3]] >= 50) {
				gtexCnt <- floor(params[[1]]*params[[3]]/100)
				gtex_shared_hubs <- gtex_m[gtex_m$Count >= gtexCnt,]
			} else {
				gtexCnt <- ceiling(params[[1]]*params[[3]]/100)
				gtex_shared_hubs <- gtex_m[gtex_m$Count <= gtexCnt,]
			}
			
			if(params[[4]] >= 50) {
				tcgaCnt <- floor((length(varNames)-params[[1]])*params[[4]]/100)
				tcga_shared_hubs <- tcga_m[tcga_m$Count >= tcgaCnt,]
			} else {
				tcgaCnt <- ceiling((length(varNames)-params[[1]])*params[[4]]/100)
				tcga_shared_hubs <- tcga_m[tcga_m$Count <= tcgaCnt,]
			}
			
			shared_hubs <- data.frame(Entrez_ID=intersect(gtex_shared_hubs$Entrez_ID, tcga_shared_hubs$Entrez_ID))
			shared_hubs$Gene_Symbol <- entrezIDtoSymbol(shared_hubs$Entrez_ID)
			shared_hubs$GTEX_Count <- gtex_shared_hubs[as.character(shared_hubs$Entrez_ID),6]
			shared_hubs$TCGA_Count <- tcga_shared_hubs[as.character(shared_hubs$Entrez_ID),6]
			
			### save the hub info
			write.table(shared_hubs, file = paste0(params[[5]], "Shared_hubs_GTEx_", params[[3]], "_TCGA_", params[[4]], ".txt"), sep = "\t", row.names = FALSE)
			
			### pathway analysis
			geneList <- shared_hubs$Gene_Symbol[which(!is.na(shared_hubs$Gene_Symbol))]
			pathwayAnalysis_TB(geneList = geneList,
					title = paste0("Shared_hubs_GTEx_", params[[3]], "(", gtexCnt, ")_TCGA_", params[[4]], "(", tcgaCnt, ")"),
					fName = paste0(params[[6]], "Pathway_shared_hubs_GTEx_", params[[3]], "_TCGA_", params[[4]], ".pdf"))
			
		} else {
			writeLines("required params do not exist")
		}
	}
	
	
	# ******************** which = viper_activity_one  *****************************
	# Generate Viper activity profiles using gene expression and Aracne network
	# Gene expression, Aracne network, and outputPath are required as "params"
	# Generates Viper activity matrix for one tissue
	#
	# The code below expects the following varaibles to be in memory: 
	# * expression matrices: emat_*
	# * regulon objects: tcga_<tissue>, <GTEX_tissue>
	#
	# There is one more option: number of samples for null model
	# Default value for the number is 15, if you don't specify it
	#
	# params[[1]]: Gene expression variable name: a character string, this needs to be the name
	#		of one of the emat_* variables.
	# params[[2]]: A character string, this must be the name of the regulon variable corresponding
	#		to params[[1]].
	# params[[3]]: Name of string vector variable containing the names of the emat_* variables. E.g.,
	#		"matNames" or another variable with smilar content (matNames is found in RDA_Files/all_64_expmat.rda.
	# params[[4]]: An integer N. The reference dataset for the VIPER computations will be constructed by 
	#		randomly choosing N samples from each expression matrix in params[[3]].
	# params[[5]]: A character string (either "zscore", "ttest", or "mean"). This is the method used for making 
	#		the viper signature and will be passed to the viperSignature call as the value of the argument "method".
	# params[[6]]: an integer. This is the number of permutations for the viper signature, it will be passed to 
	#		the viperSignature call as the value of the argument "per". 
	# params[[7]]: A character string specifying the directory where to store the results file.
	#
	# e.g., params <- list("emat_gtex_Spleen", "regulon_gtex_Spleen", "emat_gtex_names", 15, "zscore", 1000, "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/viper/")
	
	if(which == "viper_activity_one"){
		if(!is.null(params) && length(params) > 6) {
			
			### load library
			if(!require(viper)) {
				source("https://bioconductor.org/biocLite.R")
				biocLite("viper")
				library(viper)
			}
			
			### get expression
			exp <- get(as.character(params[[1]]))
			
			### get Aracne regulon
			regulon <- get(as.character(params[[2]]))
			
			### a function to combine two matrices with different rows
			### put 0 for undefined rows
			combineMats <- function(A, B) {
				if(is.null(dim(A))) {
					return(B)
				} else if(is.null(dim(B))) {
					return(A)
				} else {
					genes <- union(rownames(A), rownames(B))
					C <- data.frame(matrix(0, length(genes), (ncol(A)+ncol(B))))
					C <- cbind(data.frame(matrix(min(A), length(genes), ncol(A))), data.frame(matrix(min(B), length(genes), ncol(B))))
					rownames(C) <- genes
					colnames(C)[1:ncol(A)] <- colnames(A)
					colnames(C)[(ncol(A)+1):ncol(C)] <- colnames(B)
					C[rownames(A),1:ncol(A)] <- A
					C[rownames(B),(ncol(A)+1):ncol(C)] <- B
					
					return(C)
				}
			}
			
			### make a reference matrix
			matNames <- get(as.character(params[[3]]))
			N <- as.integer(params[[4]])
			set.seed(1234)
			refMat <- 0
			for(i in 1:length(matNames)) {
				if(matNames[i] != as.character(params[[1]])) {
					tempMat <- get(matNames[i])
					refMat <- combineMats(refMat, tempMat[,sample(ncol(tempMat),N)]) 
				}
			}
			
			### make the same rows for the test matrix and the reference
			genes <- union(rownames(refMat), rownames(exp))
			# refMat
			refMat2 <- refMat
			refMat <- data.frame(matrix(min(refMat), length(genes), ncol(refMat2)))
			rownames(refMat) <- genes
			colnames(refMat) <- colnames(refMat2)
			refMat[rownames(refMat2),] <- refMat2
			# testMat
			testMat <- data.frame(matrix(min(exp), length(genes), ncol(exp)))
			rownames(testMat) <- genes
			colnames(testMat) <- colnames(exp)
			testMat[rownames(exp),] <- exp
			
			### viper
			sigMethod <- as.character(params[[5]])
			permutNum <- as.integer(params[[6]])
			vpsig <- viperSignature(eset = as.matrix(testMat), ref = as.matrix(refMat), method = sigMethod, per = permutNum, seed = 1, verbose = FALSE)
			vpres <- viper(vpsig, regulon, verbose = FALSE)
			#vpres <- cbind(Hubs=rownames(vpres), vpres)
			
			fileName <- paste(params[[2]], "ViperMat", params[[4]], params[[5]], params[[6]], sep = "_")
			write.table(vpres, file = paste0(params[[7]], fileName, ".txt"), sep = "\t", row.names = TRUE, quote = FALSE)
		} else {
			writeLines("required params do not exist")
		}  
	}
	
	
	# ******************** which = viper_activity_all  *****************************
	# Generate Viper activity profiles using gene expression and Aracne network
	# Gene expression, Aracne network, and outputPath are required as "params"
	# Generates Viper activity matrix for all the tissues in params[[1]] and params[[2]]
	#
	# Gene expression: GTEx_36_EMat_Viper.rda
	# Aracne network: GTEx_36_Regulons.rda
	#
	# There is one more option: number of samples for null model
	# Default value for the number is 15, if you don't specify it
	#
	# params[[1]]: Expression matrices variable name
	# params[[2]]: Aracne regulons variable name
	# params[[3]]: number of samples for random selection for reference (null) model
	# params[[4]]: method when making viper signature (should be either "zscore", "ttest", or "mean")
	# params[[5]]: permutation number for viper signature
	# params[[6]]: output path
	#
	# e.g., params <- list("emat_gtex_names", "gtexRegNames", 50, "zscore", 1000, "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/viper/")
	
	if(which == "viper_activity_all") {
		if(!is.null(params) && length(params) > 5) {
			
			### load library
			if(!require(viper)) {
				source("https://bioconductor.org/biocLite.R")
				biocLite("viper")
				library(viper)
			}
			
			### a function to combine two matrices with different rows
			combineMats <- function(A, B) {
				if(is.null(dim(A))) {
					return(B)
				} else if(is.null(dim(B))) {
					return(A)
				} else {
					genes <- union(rownames(A), rownames(B))
					C <- data.frame(matrix(0, length(genes), (ncol(A)+ncol(B))))
					C <- cbind(data.frame(matrix(min(A), length(genes), ncol(A))), data.frame(matrix(min(B), length(genes), ncol(B))))
					rownames(C) <- genes
					colnames(C)[1:ncol(A)] <- colnames(A)
					colnames(C)[(ncol(A)+1):ncol(C)] <- colnames(B)
					C[rownames(A),1:ncol(A)] <- A
					C[rownames(B),(ncol(A)+1):ncol(C)] <- B
					
					return(C)
				}
			}
			
			emat_variableNames <- get(params[[1]])
			reg_variableNames <- get(params[[2]])
			
			for(i in 1:length(emat_variableNames)) {
				### get expression
				exp <- get(emat_variableNames[i])
				
				### get Aracne regulon
				regulon <- get(reg_variableNames[i])
				
				### make a reference matrix
				N <- as.integer(params[[3]])
				set.seed(1234)
				refMat <- 0
				for(j in 1:length(emat_variableNames)) {
					if(emat_variableNames[j] != get(params[[1]])[i]) {
						tempMat <- get(emat_variableNames[j])
						refMat <- combineMats(refMat, tempMat[,sample(ncol(tempMat),N)]) 
					}
				}
				
				### make the same rows for the test matrix and the reference
				genes <- union(rownames(refMat), rownames(exp))
				# refMat
				refMat2 <- refMat
				refMat <- data.frame(matrix(min(refMat), length(genes), ncol(refMat2)))
				rownames(refMat) <- genes
				colnames(refMat) <- colnames(refMat2)
				refMat[rownames(refMat2),] <- refMat2
				# testMat
				testMat <- data.frame(matrix(min(exp), length(genes), ncol(exp)))
				rownames(testMat) <- genes
				colnames(testMat) <- colnames(exp)
				testMat[rownames(exp),] <- exp
				
				### viper
				sigMethod <- as.character(params[[4]])
				permutNum <- as.integer(params[[5]])
				vpsig <- viperSignature(eset = as.matrix(testMat), ref = as.matrix(refMat), method = sigMethod, per = permutNum, seed = 1, verbose = FALSE)
				vpres <- viper(vpsig, regulon, verbose = FALSE)
				
				fileName <- paste(reg_variableNames[i], "ViperMat", params[[3]], params[[4]], params[[5]], sep = "_")
				write.table(vpres, file = paste0(params[[6]], fileName, ".txt"), sep = "\t", row.names = TRUE, quote = FALSE)
			}
		} else {
			writeLines("required params do not exist")
		}  
	}
	
	
	# ******************** which = meaningful_hubs  *****************************
	# The upgraded version of relatively_exclusive_hubs
	# This function can also consider significance of hubs
	# e.g., degree, FET & Prob scores
	#
	# Loading All_64_Aracne.rda (64 networks = 36 GTExs + 28 TCGAs) is needed before this analysis
	#
	# params[[1]]: the index in varNames that separates GTEx and TCGA (e.g., if GTEx is 1:N, enter N)
	# params[[2]]: the number of interactions that should be used for filtering out 
	#				       false-positive regulons. If NULL, no filtering will be performed.
	# params[[3]]: approach type (either "degree" or "fet" or "prob")
	# params[[4]]: approach type2 ("median" or "mean")
	# params[[5]]: the output file directory
	# params[[6]]: the number of top hubs for pathway analysis
	#
	# e.g., params <- list(36, 10, "fet", "median", "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/results/hubs/", 100)
	
	if(which == "meaningful_hubs"){
		if(!is.null(params) && length(params) > 5) {
			
			### GTEx
			gtex_hubs <- unique(unlist(sapply(varNames[1:params[[1]]], function(x){return(get(x)[[1]][,1])})))
			gtex_m <- data.frame(matrix(0, length(gtex_hubs), 12))
			colnames(gtex_m) <- c("Entrez_ID", "Gene_Symbol", "Network_Names", "Regulon_Size", "Median_Degree", "Mean_Degree", "Median_FET", "Mean_FET", "Median_Prob", "Mean_Prob", "Gene_Description", "Count")
			rownames(gtex_m) <- gtex_hubs
			gtex_m$Entrez_ID <- gtex_hubs
			gtex_m$Gene_Symbol <- entrezIDtoSymbol(gtex_hubs)
			temp <- entrezIDtoDescription(gtex_hubs)
			temp[sapply(temp, is.null)] <- "NA"
			gtex_m$Gene_Description <- unlist(temp)
			
			if(is.null(params[[2]])) {
				N <- 0
			} else {
				N <- as.numeric(params[[2]])
			}
			
			### Now it's time to check how many hubs in all interactomes
			for(net in varNames[1:params[[1]]]) {
				temp <- get(net)[[1]][,1]
				temp2 <- get(net)[[1]][,3]
				
				for(i in 1:length(temp)) {
					if(temp2[i] > N) {
						gtex_m[as.character(temp[i]),"Count"] <- gtex_m[as.character(temp[i]),"Count"]+1
						gtex_m[as.character(temp[i]), "Network_Names"] <- paste0(gtex_m[as.character(temp[i]), "Network_Names"],net,",")
						gtex_m[as.character(temp[i]), "Regulon_Size"] <- paste0(gtex_m[as.character(temp[i]), "Regulon_Size"],temp2[i],",")
					}
				}
			}
			
			for(i in 1:nrow(gtex_m)) {
				if(gtex_m[i,3] == 0) {
					gtex_m[i,3] <- NA
					gtex_m[i,4] <- NA
				} else {
					gtex_m[i,3] <- substr(gtex_m[i,3], 2, nchar(gtex_m[i,3])-1)
					gtex_m[i,4] <- substr(gtex_m[i,4], 2, nchar(gtex_m[i,4])-1)
				}
			}
			
			removeIdx <- which(is.na(gtex_m[,3]))
			if(length(removeIdx) > 0) {
				gtex_m <- gtex_m[-removeIdx,]
			}
			
			for(i in 1:nrow(gtex_m)) {
				str <- as.numeric(strsplit(gtex_m[i,4], split = ",", fixed = TRUE)[[1]])
				gtex_m[i,"Median_Degree"] <- median(str)
				gtex_m[i,"Mean_Degree"] <- mean(str)
			}
			
			
			### A function to find statistic value from tfPairEnrich[[2]] or tfPairProb[[2]]
			findValue <- function(pairProb, geneID) {
				rownames(pairProb) <- pairProb[,1]
				pairProb[which(pairProb[,2] == -Inf),2] <- min(pairProb[which(pairProb[,2] != -Inf),2])
				
				return(pairProb[geneID,2])
			}
			
			### set Median_FET & Median_Prob
			for(i in 1:nrow(gtex_m)) {
				v <- strsplit(gtex_m[i,3], split=",", fixed=TRUE)[[1]]
				
				f <- 0
				p <- 0
				cnt <- 0
				if(length(v) > 1) {
					for(j in 1:(length(v)-1)) {
						for(k in (j+1):length(v)) {
							f[cnt+1] <- findValue(tfPairEnrich[[2]][tfPairEnrich[[1]][v[j],v[k]]][[1]], as.character(gtex_m[i,1]))
							p[cnt+1] <- findValue(tfPairProb[[2]][tfPairProb[[1]][v[j],v[k]]][[1]], as.character(gtex_m[i,1]))
							cnt <- cnt+1
						}
					}
				}
				
				gtex_m[i,"Median_FET"] <- median(f)
				gtex_m[i,"Median_Prob"] <- median(p)
				gtex_m[i,"Mean_FET"] <- mean(f)
				gtex_m[i,"Mean_Prob"] <- mean(p)
			}
			
			gtex_m <- gtex_m[order(gtex_m$Count),]
			
			### save the hub info
			write.table(gtex_m, file = paste0(params[[5]], "meaningful_hubs_GTEx_", params[[2]], ".txt"), sep = "\t", row.names = FALSE)
			
			
			### TCGA
			tcga_hubs <- unique(unlist(sapply(varNames[(params[[1]]+1):length(varNames)], function(x){return(get(x)[[1]][,1])})))
			tcga_m <- data.frame(matrix(0, length(tcga_hubs), 12))
			colnames(tcga_m) <- c("Entrez_ID", "Gene_Symbol", "Network_Names", "Regulon_Size", "Median_Degree", "Mean_Degree", "Median_FET", "Mean_FET", "Median_Prob", "Mean_Prob", "Gene_Description", "Count")
			rownames(tcga_m) <- tcga_hubs
			tcga_m$Entrez_ID <- tcga_hubs
			tcga_m$Gene_Symbol <- entrezIDtoSymbol(tcga_hubs)
			temp <- entrezIDtoDescription(tcga_hubs)
			temp[sapply(temp, is.null)] <- "NA"
			tcga_m$Gene_Description <- unlist(temp)
			
			### Now it's time to check how many regulons in all interactomes
			for(net in varNames[(params[[1]]+1):length(varNames)]) {
				temp <- get(net)[[1]][,1]
				temp2 <- get(net)[[1]][,3]
				
				for(i in 1:length(temp)) {
					if(temp2[i] > N) {
						tcga_m[as.character(temp[i]),"Count"] <- tcga_m[as.character(temp[i]),"Count"]+1
						tcga_m[as.character(temp[i]), "Network_Names"] <- paste0(tcga_m[as.character(temp[i]), "Network_Names"],net,",")
						tcga_m[as.character(temp[i]), "Regulon_Size"] <- paste0(tcga_m[as.character(temp[i]), "Regulon_Size"],temp2[i],",")
					}
				}
			}
			
			for(i in 1:nrow(tcga_m)) {
				if(tcga_m[i,3]== 0) {
					tcga_m[i,3] <- NA
					tcga_m[i,4] <- NA
				} else {
					tcga_m[i,3] <- substr(tcga_m[i,3], 2, nchar(tcga_m[i,3])-1)
					tcga_m[i,4] <- substr(tcga_m[i,4], 2, nchar(tcga_m[i,4])-1)
				}
			}
			
			removeIdx <- which(is.na(tcga_m[,3]))
			if(length(removeIdx) > 0) {
				tcga_m <- tcga_m[-removeIdx,]
			}
			
			for(i in 1:nrow(tcga_m)) {
				str <- as.numeric(strsplit(tcga_m[i,4], split = ",", fixed = TRUE)[[1]])
				tcga_m[i,"Median_Degree"] <- median(str)
				tcga_m[i,"Mean_Degree"] <- mean(str)
			}
			
			### set Median_FET & Median_Prob
			for(i in 1:nrow(tcga_m)) {
				v <- strsplit(tcga_m[i,3], split=",", fixed=TRUE)[[1]]
				
				f <- 0
				p <- 0
				cnt <- 0
				if(length(v) > 1) {
					for(j in 1:(length(v)-1)) {
						for(k in (j+1):length(v)) {
							f[cnt+1] <- findValue(tfPairEnrich[[2]][tfPairEnrich[[1]][v[j],v[k]]][[1]], as.character(tcga_m[i,1]))
							p[cnt+1] <- findValue(tfPairProb[[2]][tfPairProb[[1]][v[j],v[k]]][[1]], as.character(tcga_m[i,1]))
							cnt <- cnt+1
						}
					}
				}
				
				tcga_m[i,"Median_FET"] <- median(f)
				tcga_m[i,"Median_Prob"] <- median(p)
				tcga_m[i,"Mean_FET"] <- mean(f)
				tcga_m[i,"Mean_Prob"] <- mean(p)
			}
			
			tcga_m <- tcga_m[order(tcga_m$Count),]
			
			### save the hub info
			write.table(tcga_m, file = paste0(params[[5]], "meaningful_hubs_TCGA_", params[[2]], ".txt"), sep = "\t", row.names = FALSE)
			
			
			### get shared hubs between GTEx and TCGA
			shared_hubs <- intersect(rownames(gtex_m), rownames(tcga_m))
			shared_gtex <- gtex_m[shared_hubs,]
			shared_tcga <- tcga_m[shared_hubs,]
			
			
			### make a specific file
			if(as.character(params[[3]]) == "degree") {
				
				if(as.character(params[[4]]) == "median") {
					v <- abs(shared_gtex$Median_Degree - shared_tcga$Median_Degree)
					
					shared_gtex <- shared_gtex[order(-v),]
					shared_tcga <- shared_tcga[order(-v),]
					
					shared_m <- data.frame(cbind(Entrez_ID=shared_gtex$Entrez_ID, Gene_Symbol=as.character(shared_gtex$Gene_Symbol),
									GTEX_Degree=shared_gtex$Median_Degree, TCGA_Degree=shared_tcga$Median_Degree,
									Diff=abs(shared_gtex$Median_Degree - shared_tcga$Median_Degree)))
				} else if(as.character(params[[4]]) == "mean") {
					v <- abs(shared_gtex$Mean_Degree - shared_tcga$Mean_Degree)
					
					shared_gtex <- shared_gtex[order(-v),]
					shared_tcga <- shared_tcga[order(-v),]
					
					shared_m <- data.frame(cbind(Entrez_ID=shared_gtex$Entrez_ID, Gene_Symbol=as.character(shared_gtex$Gene_Symbol),
									GTEX_Degree=shared_gtex$Mean_Degree, TCGA_Degree=shared_tcga$Mean_Degree,
									Diff=abs(shared_gtex$Mean_Degree - shared_tcga$Mean_Degree)))
				}
				
			} else if(as.character(params[[3]]) == "fet") {
				
				if(as.character(params[[4]]) == "median") {
					v <- abs(shared_gtex$Median_FET - shared_tcga$Median_FET)
					
					shared_gtex <- shared_gtex[order(-v),]
					shared_tcga <- shared_tcga[order(-v),]
					
					shared_m <- data.frame(cbind(Entrez_ID=shared_gtex$Entrez_ID, Gene_Symbol=as.character(shared_gtex$Gene_Symbol),
									GTEX_Degree=shared_gtex$Median_FET, TCGA_Degree=shared_tcga$Median_FET,
									Diff=abs(shared_gtex$Median_FET - shared_tcga$Median_FET)))
				} else if(as.character(params[[4]]) == "mean") {
					v <- abs(shared_gtex$Mean_FET - shared_tcga$Mean_FET)
					
					shared_gtex <- shared_gtex[order(-v),]
					shared_tcga <- shared_tcga[order(-v),]
					
					shared_m <- data.frame(cbind(Entrez_ID=shared_gtex$Entrez_ID, Gene_Symbol=as.character(shared_gtex$Gene_Symbol),
									GTEX_Degree=shared_gtex$Mean_FET, TCGA_Degree=shared_tcga$Mean_FET,
									Diff=abs(shared_gtex$Mean_FET - shared_tcga$Mean_FET)))
				}
				
			} else if(as.character(params[[3]]) == "prob") {
				
				if(as.character(params[[4]]) == "median") {
					v <- abs(shared_gtex$Median_Prob - shared_tcga$Median_Prob)
					
					shared_gtex <- shared_gtex[order(-v),]
					shared_tcga <- shared_tcga[order(-v),]
					
					shared_m <- data.frame(cbind(Entrez_ID=shared_gtex$Entrez_ID, Gene_Symbol=as.character(shared_gtex$Gene_Symbol),
									GTEX_Degree=shared_gtex$Median_Prob, TCGA_Degree=shared_tcga$Median_Prob,
									Diff=abs(shared_gtex$Median_Prob - shared_tcga$Median_Prob)))
				} else if(as.character(params[[4]]) == "mean") {
					v <- abs(shared_gtex$Mean_Prob - shared_tcga$Mean_Prob)
					
					shared_gtex <- shared_gtex[order(-v),]
					shared_tcga <- shared_tcga[order(-v),]
					
					shared_m <- data.frame(cbind(Entrez_ID=shared_gtex$Entrez_ID, Gene_Symbol=as.character(shared_gtex$Gene_Symbol),
									GTEX_Degree=shared_gtex$Mean_Prob, TCGA_Degree=shared_tcga$Mean_Prob,
									Diff=abs(shared_gtex$Mean_Prob - shared_tcga$Mean_Prob)))
				}
				
			}
			
			shared_m <- shared_m[0:as.integer(params[[6]]),]
			
			write.table(shared_m, file = paste0(params[[5]], "meaningful_hubs_", params[[4]], "_", params[[3]], "_", as.integer(params[[6]]), ".txt"), sep = "\t", row.names = FALSE)
			
			### pathway analysis
			geneList <- shared_m$Gene_Symbol[which(!is.na(shared_m$Gene_Symbol))]
			pathwayAnalysis_TB(geneList = geneList, displayNum = 50,
					title = paste0("meaningful_hubs_", params[[4]], "_", params[[3]], "_", as.integer(params[[6]])),
					fName = paste0(params[[5]], "meaningful_hubs_", params[[4]], "_", params[[3]], "_", as.integer(params[[6]]), ".pdf"))
			
		} else {
			writeLines("required params do not exist")
		}  
	}
	
	
	# ******************** which = exclusive_interactions  *****************************
	# Here we generate a file with data to be used in the methods:
	#			make_graphs("exclusive_interactions")
	#			oneOffs("exclusive_interactions_GO")
	#
	# Our goal is to identify interactions with high support exclusively in GTEx or TCGA networks.
	# The support of an interaction is defined as the number of networks that contain it
	# (considering only inteactions that meet the user-provided MI and p-value thresholds
	# described below). Our assumption is that interactions with high support are "real".
	# Further, it seems reasonable to assume that if such interactions appear only in the
	# GTEx netrworks (i.e., the have zero support in the TCGA networks) they could
	# participate in tumor suppresing processes. Analogously, high-support interactions
	# that are exclusive to the TCGA networks could participate in tumor promoting 
	# processes.
	#
	# The code below assumes that params[[]] contains the following mandatory elements:
	# * params[[1]]:	an integer indicating the maximum number "top" of exlusive 
	#			interactions per hub gene to report.
	# * params[[2]]:	a real number providing the MI threshold. Interactions with MI value 
	#			less then this threshold are not considered in the processing (i.e., they 
	#			are treated as if they do not exist).
	# * params[[3]]:	a real number providing the p-value threshold. Interactions with P-value 
	#			higher then this threshold are not considered in the processing (i.e., they 
	#			are treated as if they do not exist).
	# * params[[4]]:	an integer indicating the number of network permuations to use for 
	# 			generating null distributions. See documentations of method computeDiffinteractions()
	#			below.
	# * params[[5]]:	a string providing the name of the file where to store the results.
	#
	# The following params[[]] values are optional: they can either be all ommitted (in which case
	# defaults will be used) or the must all be present.
	# * params[[6]]:	If this parameter is not provided, the analysis below will operate on the 
	#			network variables whose names are provided in the variable "varNames". Further, it
	#			will assume that varNames[1:36] are the GTEx networks, which varNames[37:64] are the
	#			TCGA networks. If this parameter is provided, it should be a vector like varNames
	#			containing the names of interactome variables to use for the analysis:
	# * params[[7]]:	A vector of integers, providing the indices within params[[6]] that contain
	#			the names of the GTEx networks.
	# * params[[8]]:	A vector of integers, providing the indices within params[[6]] that contain
	#			the names of the TCGA networks.
	#
	#
	# The code below tabulates the exclusive interactions in the GTEx and TCGA networks as well as in
	# params[[4]] random shufflings and saves the results in the file params[[5]]. See documentation
	# of the methods below for exact details. 
	# *****ATTENTION*****: This is a coputatationally intensive processing, computing the null for 
	# each permutation can take 5-10 minutes. So, if params[[4]] is 100 or so, this will take several
	# hours. This is why we save the results in a file for further processing. At the time of writing
	# these comments, the file generated by running this methods is stored here:
	#			/ifs/archive/shares/af_lab/GTEx/IntermediateFiles/Aris/ExclusiveInteractions/randomized_exclusive.rda
	
	if(which == "exclusive_interactions"){
		
		README <- function(){
			writeLines("The file contains data generated to test if interactions that appear only in the GTEx")
			writeLines("(or only in the TCGA) networks are a fluke of they are statistically more significant")
			writeLines("than what would be expected by chance. The files contains the following variables\n")
			writeLines("* L_actual:\t\tThis is a list with two elements:")
			writeLines("\t* The first is a list with one entry per hub gene seen in both GTEx and TCGA")
			writeLines("\t\tnetworks. For each such hub genes H it contains a table with targets of H")
			writeLines("\t\tgene seen in at least one GTEx network but never seen in any TCGA network. For")
			writeLines("\t\teach such target T the table entry is the number N of GTEx interactomes where")
			writeLines("\t\tthe interaction (H,T) is seen. The table is ordered in decreasing size of N.")
			writeLines("\t\tOnly the top 100 genes are returned. Only interactions with MI above a given ")
			writeLines("\t\tthreshold and p-value below a given threshold are considered in the counting,")
			writeLines("\t\tall other interactions are treated as if the do not exist. In the run that Aris")
			writeLines("\t\tdid he used the following thresholds, respectively: MI >= 0.1 and p-value <= 0.0001.")
			writeLines("\t\tThese are conservative, they remove from consideration almost half interactions in")
			writeLines("\t\teach interactome.")
			writeLines("\t* The second is similar to the first but reporting interactions seen exclusively")
			writeLines("\t\tin TCGA networks.")
			writeLines("* L_random_<i>:\t\tSame structure as L_actual bur computed on random GTEx/TCGA")
			writeLines("\t\tnetworks, acquired by shuffling networks in either the GTEx or TCGA group.")
			writeLines("\t\t100 random assigments were used, i.e., <i> runs from 1 - 100.")
			writeLines("")
		}
		
		# *****************************************************************************
		# Identify targets of hub gene "gene" that are present only in its GTEx 
		# regulons but in none of its TCGA networks (and vise versa).
		# 
		# ARGUMENTS
		# * gene:	the query hub gene.
		# * top:	specifies how many of the targers that meet the criteria to return.
		# * comp:	specifies the direction of the comparison. A value of "GTEx" will
		#			return interactions that appear in GTEx but not in TCGA. A value 
		#			of "TCGA" will return interactions that appear in TCGA but not in GTEx.
		# * mi:		mutual information threshold.
		# * pval:	p-valune treshold.
		#
		# RETURN VALUE
		# Assuming that comp="GTEx", the method will visit the regulons of "gene" in
		# all GTEx networks and will collect the targets of all interactions with
		# mututal information >= "mi" and p-value <= "pval". It will then do the same
		# for all regulons of "gene" in TCGA. It will then take a diff of the two sets
		# of targets, retaining targets that appear only in GTEx regulons but in none 
		# of the TCGA networks. It will then tabulate the results, returning a table
		# with one entry for each gene in the diff set, where the the value of the 
		# table is the number of GTEx networks found to contain the ("gene", target)
		# interaction. Table entries are sorted in decreasing order and named with the
		# Entrez Id of the respective target. If the diff set is empty, then NULL
		# will be returned. Otherwise, the "top" interactions with the highest support
		# will be returned. The behavior is symetrical when comp="TCGA", i.e., the
		# method returns a table for targets of "gene" that appear only in TCGA
		# networks.
		# *****************************************************************************
		getTopDiffInteractions <- function(gene, nets, top , comp = c("GTEx", "TCGA"), mi = 0, pval = 1){
			x_g = tabulateRegulonInteractions(gene, nets = nets[gtex_indices], mi = mi, pval = pval)
			x_t = tabulateRegulonInteractions(gene, nets = nets[tcga_indices], mi = mi, pval = pval)
			if (comp == "GTEx"){
				if (length(x_g) == 0)
					return(NULL)
				if (length(x_t) == 0)
					return(x_g[1:min(length(x_g), top)])
				res = setdiff(names(x_g), names(x_t))
				if (length(res) == 0)
					return(NULL)
				return(x_g[res][1:min(length(res), top)])
			}
			else if (comp == "TCGA"){
				if (length(x_t) == 0)
					return(NULL)
				if (length(x_g) == 0)
					return(x_t[1:min(length(x_t), top)])
				res = setdiff(names(x_t), names(x_g))
				if (length(res) == 0)
					return(NULL)
				return(x_t[res][1:min(length(res), top)])
			}
		}
		
		
		# *****************************************************************************
		# For every hub gene seen in at least one GTEx network and at least one TCGA
		# network compile two sorted tables: one listing interactions with high
		# support in GTEx networks but zeso support in TCGA networks; and one listing
		# interactions with high support in TCGA networks but zeso support in GTEx 
		# networks. It does so for either the real netwrorks, or permuted versions (to
		# be used in generating appropriate nulls).
		#
		# ARGUMENTS:
		# * mode:	indicates if the tables should be generatred for the actual GTEx/TCGA
		#			networks, or for permuted versions. In the latter case, the network
		#			labels (GTEx or TCGA) are randomlhy shuffled.
		# * top:	indicates how many intereactions to include in the returned results.
		#			For each hub gene, interactions will be sorted in decreasing order
		#			of their support and the "top" ones will be included in the final
		#			sorted table. If a given hub has less that "top" interactions
		#			meeting the threshold filters listed below, only those are returned.
		#			
		# * mi:		MI threshold. Interactions with MI value < "mi" are not considered
		#			in the processing (i.e., they are treated as if they do not exist).
		# * pval:	P-value threshold. Interactions with MI value > "pval" are not considered
		#			in the processing (i.e., they are treated as if they do not exist).
		#
		# RETURN VLAUE:
		# A list with two elements:
		# * L_gtex: A named list with one entry per gene hub. Each entry contains an ordered 
		#			table of size "top" (or less, see description of "top" argument above). 
		#			The i-th entry contains a named table T listing interactions for hub gene 
		#			R = names(L_gtex)[[i]]. The j-th entry in T corresponds to a target 
		#			G = names(T)[j] such that the interaction (R,G) is seen only in GTEx
		#			networks but in none of the TCGA networks. The value T[j] is the support
		#			of the interaction (R,G) in GTEx, i.e., the number of networks containing it.
		# * L_tcga:	Same as L_gtex but containing interactions seen only in TCGA networks and 			
		# 			having 0 support in GTEx networks.
		# *****************************************************************************
		computeDiffinteractions <- function(mode=c("actual", "randomize"), nets, top, mi = 0, pval = 1){
			if (mode == "randomize")
				nets = nets[sample(1:length(nets))]
			common_hubs = intersect(unique(unlist(sapply(nets[gtex_indices], function(x){return(get(x)[[1]][,1])}))), 
					unique(unlist(sapply(nets[tcga_indices], function(x){return(get(x)[[1]][,1])}))))
			
			hubs = common_hubs
			# Get the unique GTEx interactions
			L_gtex = lapply(hubs, function(g){
						return(getTopDiffInteractions(g, nets = nets, comp="GTEx", top = top, mi = mi, pval = pval))
					})
			names(L_gtex) = hubs
			# Remove NULLs
			L_gtex = L_gtex[!sapply(L_gtex, is.null)]
			
			# Get the unique TCGA interactions
			L_tcga = lapply(hubs, function(g){
						return(getTopDiffInteractions(g, nets = nets, comp="TCGA", top = top, mi = mi, pval = pval))
					})
			names(L_tcga) = hubs
			# Remove NULLs
			L_tcga = L_tcga[!sapply(L_tcga, is.null)]
			
			res = list(L_gtex, L_tcga)
			names(res) = c("L_gtex", "L_tcga")
			return(res)	
		}
		
		if (is.null(params) | (length(params) != 5 & length(params) != 8))
			stop("params[[]] needs to be a list with five or eight members.")
		
		top = params[[1]]
		mi = params[[2]]
		pval = params[[3]]
		iterations = params[[4]]
		results_file = params[[5]]
		nets = varNames
		gtex_indices = 1:36
		tcga_indices = 37:length(nets)
		
		
		if (length(params) == 8){
			nets = params[[6]]
			gtex_indices = params[[7]]
			tcga_indices = params[[8]]
		}
		
		# Tabulate support for exlcusive interactions (eihter only in GTEx or only
		# in TCGA) for the actual GTEx and TCGA network collections.
		L_actual = computeDiffinteractions(mode = "actual", nets = nets, top = top, mi = mi, pval = pval)
		
		# Shuffle network assignements; for each shuffle, tabulate support for exclusive
		# interactions. These results will be used as a null, for assessing the significance
		# of coordinated regulation of exclusive interactions.
		if (iterations >0 ){
			for (i in 1:iterations)
				assign(paste("L_random_", i, sep = ""), 
						computeDiffinteractions(mode = "randomize", nets = nets, top = top, mi = mi, pval = pval), 
						envir = globalenv())
			vars = c("L_actual", paste("L_random_", seq(1,iterations), sep=""), "README")
		} else
			vars = c("L_actual", "README")
		
		# Save to file for later use, as the computations above are time consuming
		save(list = vars, file = results_file)
	}
	
	
	# ******************** which = exclusive_interactions_GO  **************************
	# Identifies genes participating in GTEx-exclusive interactions with unexpectedly 
	# high support (i.e., interactions that appear in many GTEx networks but in no TCGA
	# networks) and performs GO enrichment analysis on them. It then repeats the same
	# computation for TCGA-exlusive interactions.
	
	# The code below expects the following values to be passed in the params[[]] argument:
	# * params[[1]]:	Character string containing the full pathname of the file 
	#			produced by the call to oneOffs("exclusive_interactions", params), i.e. 
	#			to the file specified in params[[5]] in the parameter list passed to the
	#			call to oneOffs.
	# * params[[2]]:	significance threshold "alpha" for Bonferroni correction. We only 
	#			consider interactions with support at least N such that bonferroni-corrected 
	#			probability of seeing such a support in the random data is no more than
	#			this value. This is a number between 0 and 1.
	# * params[[3]]: name of CSV file where to store the results of the GO enrichment analysis
	# 		for the GTEx-exclusive interactions
	# * params[[4]]: name of CSV file where to store the results of the GO enrichment analysis
	# 		for the TCGA-exclusive interactions
	# * params[[5]]: name of CSV file where to store the results of the GO enrichment analysis
	# 		for the intersetion of genes found both in TCGA-exclusive and GTEX-exclusive
	# 		interactions.
	#
	# METHOD OPERATION
	# It examines the file generated by calling the method oneOffs("exclusive_interactions")
	# and counts the number M of GTEx-exlusive interactions. It then builds the null density 
	# of support values for the interactions in the randomized runs. It then uses uses 
	# the null density to identify the support level N such that the probability t one out
	# of M randomly drawn GTEx-exclusive interactions is larger than N is alpha. I.e., it
	# Bonferroni-corrects the type I error according to the number of actual GTEx-exclusive
	# interactions. Finally, it compiles all actual GTEx-exclusive interactions with suport
	# at least N+1, extracts all genes participating in these interactions, and runs GO
	# enrichment analysis on them. The same process is followed for the TCGA-exclusive 
	# interactions.
	#
	# Essentially, the method identifies interactions of high confidence via their high
	# support (and their high MI and low p-value used when calling oneOff("exclusive_interactions"))
	# that are mutually exlusive to the GTEx or TCGA netwoks. Given their high support, these 
	# should be ivolved in pan-cancer behavior and could represent dysregulated global tumor suppression
	# programs (in the case of GTEx-exclusive interactions) or tumor promotion programs (in the case
	# of TCGA-exclusive interactions).
	if (which == "exclusive_interactions_GO"){
		
		# *****************************************************************************
		# Return the list of all genes participating in interactions clearing a support
		# threshold.
		#
		# ARGUMENTS
		# * L:		a list of tables, in the format of L_actual[[1]] or L_actual[[2]],
		#			comprising support data for GTEx-exclusive ot TCGA-exclusive interactions.
		#			See documentation of method oneOffs("exclusive_interactions") for 
		#			the structure of these variables.
		# * min:	minimum support threshold for an interaction to clear
		#
		# RETURN VALUE
		# A vector comprising the entrez ids of all genes G such that G is a member of
		# a GTEx- ot TCGA-exclusive interaction with support at list "min". Both targets
		# and hubs are included.
		tabulateCommonRegulators <- function(L, min = 15){
			LL = unlist(sapply(names(L), function(gid){
								x = L[[gid]]
								if (x[1] < min)
									return(NA)
								else
									return(c(gid, names(x[x >= min])))}))
			LL = LL[!is.na(LL)]
			LL = unique(LL)
		}
		
		if (!exists("L_actual"))
			load(params[[1]])
		
		# Get the support of all GTEx-exclusive interactions from the randomized runs,
		# to build the null distribution of the support values
		random_runs = length(ls(pattern="L_random_"))		
		rg = unlist(sapply(1:random_runs, function(x){
							v = get(paste("L_random", x, sep="_"))
							r = v[[1]]
							return( unlist(sapply(r, function(x){return(x)})))
						}))
		
		# Next collect the support of the actual GTEx-exclusive interactions
		g = L_actual[[1]]
		gc = unlist(sapply(g, function(x){return(x)}))
		
		# Find the minimum support for the requested bonferonni p-value
		N = ceiling(quantile(rg, probs=(1-params[[2]]/length(gc))))
		# Select GTEx-exlusive interactions with support larger than N,
		# run GO enrichment analysis, and save the results to file
		tg = tabulateCommonRegulators(L_actual[[1]], N+1)
		resg = goEnrichment(tg)
		write.csv(resg[[1]], file = params[[3]], row.names=FALSE)
		
		# Now do same for the TCGA-exclusive interactions
		rt = unlist(sapply(1:random_runs, function(x){
							v = get(paste("L_random", x, sep="_"))
							r = v[[2]]
							return( unlist(sapply(r, function(x){return(x)})))
						}))
		t = L_actual[[2]]
		tc = unlist(sapply(t, function(x){return(x)}))
		
		N = ceiling(quantile(rt, probs=(1-params[[2]]/length(tc))))
		# Select TCGA-exlusive interactions with support larger than N
		tt = tabulateCommonRegulators(L_actual[[2]], N+1)
		rest = goEnrichment(tt)
		write.csv(rest[[1]], file = params[[4]], row.names=FALSE)
		
		# Finally, for good measure, run GO enrichement analysis on genes found
		# in both TCGA-exclusive and GTEx-exclusive interactions.
		resc = goEnrichment(intersect(tg, tt))
		write.csv(resc[[1]], file = params[[5]], row.names=FALSE)
	}
	
	
	# ******************** which = MakeViperRDA  *****************************
	# Load Viper matrices and save them as RDA file
	# The variable varNamesVP represents all the names of Viper matrices
	#
	# params[[1]]: the directory of Viper matrices in tab-separated files
	# params[[2]]: the output RDA file path
	#
	# e.g., params <- list("//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/viper/GTEx/", "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/GTEx_36_ViperMats.rda")
	
	if(which == "MakeViperRDA"){
		if(!is.null(params) && length(params) > 1) {
			
			### collect files from the fileNamePath
			f <- list.files(as.character(params[[1]]))
			f <- f[which(endsWith(f, ".txt"))]
			
			### set varNames
			varNamesVP <- sapply(f, function(x) paste("vmat", paste(strsplit(x, split = "_")[[1]][2:3], collapse = "_"), sep = "_"), USE.NAMES = FALSE)
			
			### assign viper matrices as R object
			sapply(1:length(f), function(i, f, varNamesVP) {
						assign(varNamesVP[i],
								as.matrix(read.table(paste0(params[[1]], f[i]),
												header = TRUE, sep = "\t", row.names = 1,
												stringsAsFactors = FALSE, check.names = FALSE)),
								envir = globalenv())
					}, f=f, varNamesVP=varNamesVP,
					USE.NAMES = FALSE)
			
			### set README function
			README <- function(){
				writeLines(paste(rep("#", 100), collapse = ""))
				writeLines("Each vmat object has Viper activity scores.")
				writeLines("varNamesVP contains all the variable names of the Viper objects.")
				writeLines("For each tissue, 50 samples from each of the remaining tissues were randomly collected and used as referenece.")
				writeLines("When generating Viper signature, z-score was chosen.")
				writeLines(paste(rep("#", 100), collapse = ""))
			}
			
			### save as RDA
			save(list = c(varNamesVP, "varNamesVP", "README"), file = as.character(params[[2]]))
			
		} else {
			writeLines("required params do not exist")
		}  
	}
	
	# ******************** which = DensityPerHubs  *****************************
	# Load Viper RDA file and draw a density line for each hub between similar tissues
	# one from GTEx and the other one from TCGA
	#
	# params[[1]]: one viper matrix variable name in GTEx
	# params[[2]]: one viper matrix variable name in TCGA
	# params[[3]]: the number of the top ranking hubs that are differentially activated
	# params[[4]]: output directory path
	#
	# e.g., params <- list("vmat_gtex_Breast", "vmat_tcga_brca", 1, "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/viper/density/")
	
	if(which == "DensityPerHubs"){
		if(!is.null(params) && length(params) > 3) {
			
			### load library
			if(!require(ggplot2)) {
				install.packages("ggplot2")
				library(ggplot2)
			}
			
			### get viper matrices      
			gtex_mat <- get(as.character(params[[1]]))
			tcga_mat <- get(as.character(params[[2]]))
			
			### since the viperMats may have different set of hubs, get the commons
			common_hubs <- intersect(rownames(gtex_mat), rownames(tcga_mat))
			gtex_mat <- gtex_mat[common_hubs,]
			tcga_mat <- tcga_mat[common_hubs,]
			
			### get differentially activated hubs between GTEx and TCGA
			t_gtex <- apply(gtex_mat, 1, function(x) mean(x)/sd(x))
			t_tcga <- apply(tcga_mat, 1, function(x) mean(x)/sd(x))
			t <- abs(t_gtex-t_tcga)
			t <- t[order(-t)]
			
			### select the top DA hubs 
			top_hubs <- names(t)[1:as.numeric(params[[3]])]
			
			### make a data for density plot
			v <- NULL
			for(i in 1:length(top_hubs)) {
				v <- c(v, as.numeric(gtex_mat[top_hubs[i],]))
				names(v)[((i-1)*(ncol(gtex_mat)+ncol(tcga_mat))+1):length(v)] <- paste("GTEx", top_hubs[i], sep = "_")
				tmpIdx <- length(v)+1
				v <- c(v, as.numeric(tcga_mat[top_hubs[i],]))
				names(v)[tmpIdx:length(v)] <- paste("TCGA", top_hubs[i], sep = "_")
			}
			
			### prepare data frame for ggplot
			dat <- data.frame(NES=v, Hub=names(v))
			
			### file name
			fName <- paste(c(strsplit(as.character(params[[1]]), "_")[[1]][2:3],
							strsplit(as.character(params[[2]]), "_")[[1]][2:3]), collapse = "_")
			
			### density plot
			ggplot(dat, aes(x=NES, col=Hub)) + geom_density(alpha = 0.5) +
					theme_classic(base_size = 16) +
					labs(title=fName, subtitle="Density per hub")
			
			### save the plot
			ggsave(filename = paste0(params[[4]], fName, "_density.png"),
					width = 15, height = 10)
			
		} else {
			writeLines("required params do not exist")
		}
	}
	
	# ******************** which = viper_correlation  *****************************
	# Load Viper RDA file and draw a correlation plot between GTEx and TCGA (specific tissue)
	# Each dot is a hub and NES if calculated as mean/std across all samples
	# Since the data is not from msVIPER, each sample has its own NES
	#
	# params[[1]]: one viper matrix variable name in GTEx
	# params[[2]]: one viper matrix variable name in TCGA
	# params[[3]]: the number of the top ranking hubs that are differentially activated
	# params[[4]]: output directory path
	#
	# e.g., params <- list("vmat_gtex_Breast", "vmat_tcga_brca", 10, "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/viper/correlation/")
	
	if(which == "viper_correlation"){
		if(!is.null(params) && length(params) > 3) {
			
			### load library
			if(!require(ggrepel)) {
				install.packages("ggrepel")
				library(ggrepel)
			}
			
			### get viper matrices      
			gtex_mat <- get(as.character(params[[1]]))
			tcga_mat <- get(as.character(params[[2]]))
			
			### since the viperMats may have different set of hubs, get the commons
			common_hubs <- intersect(rownames(gtex_mat), rownames(tcga_mat))
			gtex_mat <- gtex_mat[common_hubs,]
			tcga_mat <- tcga_mat[common_hubs,]
			
			### get mean/std across all samples for each hub
			t_gtex <- apply(gtex_mat, 1, function(x) mean(x)/sd(x))
			t_tcga <- apply(tcga_mat, 1, function(x) mean(x)/sd(x))
			
			### change names from Entrez ID to Gene symbol
			names(t_gtex) <- entrezIDtoSymbol(names(t_gtex))
			names(t_tcga) <- entrezIDtoSymbol(names(t_tcga))
			
			### create a data frame for the plot
			sig <- as.character(vector(mode = "character", length = length(t_gtex)))
			df <- data.frame(t_gtex, t_tcga, sig)
			
			### characterize the factors
			idx <- sapply(df, is.factor)
			df[idx] <- lapply(df[idx], function(x) as.character(x))
			
			### mark the top differentially activated hubs
			t <- abs(t_gtex-t_tcga)
			t <- t[order(-t)]
			top_hubs <- names(t)[1:as.numeric(params[[3]])]
			df[top_hubs,3] <- top_hubs
			
			### file name
			fName <- paste(c(strsplit(as.character(params[[1]]), "_")[[1]][2:3],
							strsplit(as.character(params[[2]]), "_")[[1]][2:3]), collapse = "_")
			
			### make a cor plot
			ggplot(data = df, aes(x=t_gtex, y=t_tcga)) +
					geom_point(color = "black", size = 1) +
					geom_label_repel(aes(t_gtex, t_tcga, label = sig), color = "red", box.padding = unit(0.45, "lines")) +
					labs(title=paste0(fName, "_correlation"),
							subtitle=sprintf("P.Cor = %s, p-value = %s", round(cor(df$t_gtex, df$t_tcga), 5),
									signif(cor.test(df$t_gtex, df$t_tcga)$p.value, 5)),
							t_gtex="GTEx", t_tcga="TCGA") +
					geom_smooth(method = lm, color="gray", se=FALSE) +
					theme_classic(base_size = 16)
			
			### save the plot
			ggsave(filename = paste0(params[[4]], fName, "_correlation.png"), width = 10, height = 10)
			
		} else {
			writeLines("required params do not exist")
		}
	}
	
  
  # ***************************** which = viper_densities ********************************
  # Make a density plot with viper activity scores
  # It would be useful to know difference of sample heterogeneity between GTEx and TCGA
  # params[[1]]: a character vector of variable names of viper matrices
  # params[[2]]: output directory path for density plots of all the tissues
  # e.g., params=list("varNamesVP", "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/viper/density/")
  
  if(which == "viper_densities"){
    if(!is.null(params && length(params) > 1)) {
      
      ### load library
      if(!require(ggplot2)) {
        install.packages("ggplot2")
        library(ggplot2)
      }
      
      ### get viper matrix names
      viperMatNames <- get(as.character(params[[1]]))
      
      ### since all the viperMats have different set of hubs, get the commons
      common_hubs <- rownames(get(viperMatNames[1]))
      for(i in 2:length(viperMatNames)) {
        common_hubs <- intersect(common_hubs, rownames(get(viperMatNames[i])))
      }
      
      ### print density plots of all the tissues
      for(i in 1:length(viperMatNames)) {
        ### load the viper matrix and select only the common hubs
        m <- get(viperMatNames[i])
        m <- m[common_hubs,]
        
        ### make the viper matrix to one concatenated vector
        v <- NULL
        for(j in 1:ncol(m)) {
          v <- c(v, m[,j])
          names(v)[((j-1)*nrow(m)+1):length(v)] <- colnames(m)[j]
        }
        
        ### prepare data frame for ggplot
        dat <- data.frame(NES=v, Sample=names(v))
        
        ### density plot
        ggplot(dat, aes(x=NES, col=Sample)) + geom_density(alpha = 0.5) +
          theme_classic(base_size = 16) +
          labs(title=viperMatNames[i], subtitle="") +
          theme(legend.position="none")
        
        ### save the plot
        ggsave(filename = paste0(params[[2]], viperMatNames[i], "_density.png"),
               width = 15, height = 10)
      }
      
      ### get median values for each hub for each tissue
      ### combine all the tissues into one
      v <- NULL
      for(i in 1:length(viperMatNames)) {
        m <- get(viperMatNames[i])
        m <- m[common_hubs,]
        v <- c(v, apply(m, 1, median))
        names(v)[((i-1)*nrow(m)+1):length(v)] <- viperMatNames[i]
      }
      
      ### prepare data frame for ggplot
      dat <- data.frame(NES=v, Tissue=names(v))
      
      ### density plot
      ggplot(dat, aes(x=NES, col=Tissue)) + geom_density(alpha = 0.5) +
        theme_classic(base_size = 16) +
        labs(title=paste("All", length(viperMatNames), "Density Plot"), subtitle="")
      
      ### save the plot
      ggsave(filename = paste0(params[[2]], "all_", length(viperMatNames), "_density.png"),
             width = 25, height = 10)
      
      ### GTEx vs TCGA
      names(v) <- sapply(strsplit(names(v), split = "_", fixed = TRUE), function(x) x[2])
      
      ### prepare data frame for ggplot
      dat <- data.frame(NES=v, Tissue=names(v))
      
      ### density plot
      ggplot(dat, aes(x=NES, col=Tissue)) + geom_density(alpha = 0.5) +
        theme_classic(base_size = 16) +
        labs(title=paste("All", length(viperMatNames), "GTEx-TCGA Density Plot"), subtitle="")
      
      ### save the plot
      ggsave(filename = paste0(params[[2]], "all_", length(viperMatNames), "gtex_tcga_density.png"),
             width = 15, height = 10)
      
    } else {
      stop("required params do not exist")
    }
  }
  
	
	# ******************** which = hub_investigation *****************************
	# Pick one hub of interest and investigate about it
	# 1. Present target genes of the hub in GTEx or in TCGA
	#    as well as intersections and exclusives
	# 2. Pathway analysis of the target genes
	# 3. Heatmap to show gene expression changes of the targets between GTEx and TCGA
	# Aracne and Gene expression needed
	#
	# Needs all_64_aracne_mi_fixed.rda, all_64_expmat.rda, and all_64_vipersigs.rda loaded
	#
	# params[[1]]: Hub Entrez (NCBI) ID
	# params[[2]]: GTEx Aracne object name
	# params[[3]]: TCGA Aracne object name
	# params[[4]]: GTEx expression matrix name
	# params[[5]]: TCGA expression matrix name
	# params[[6]]: GTEx Viper signature object name
	# params[[7]]: TCGA Viper signature object name
	# params[[8]]: Correlation threshold when building a network graph
	# params[[9]]: Output directory path
	#
	# e.g., params <- list("159119", "Breast", "tcga_brca",
	#                      "emat_gtex_Breast", "emat_tcga_brca",
	#                      "vipersig_gtex_Breast","vipersig_tcga_brca",
	#                      0.8, "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/results/cancer_hubs/")
	# e.g., params <- list("255738", "Prostate", "tcga_prad",
	#                      "emat_gtex_Prostate", "emat_tcga_prad",
	#                      "vipersig_gtex_Prostate","vipersig_tcga_prad",
	#                      0.85, "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/results/cancer_hubs/")
	# e.g., params <- list("5777", "Lung", "tcga_lusc",
	#                      "emat_gtex_Lung", "emat_tcga_lusc",
	#                      "vipersig_gtex_Lung","vipersig_tcga_lusc",
	#                      0.8, "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/results/cancer_hubs/")
	
	if(which == "hub_investigation"){
		if(!is.null(params) && length(params) > 8) {
			
			### get objects
			gtex_aracne <- get(as.character(params[[2]]))
			tcga_aracne <- get(as.character(params[[3]]))
			gtex_exp <- get(as.character(params[[4]]))
			tcga_exp <- get(as.character(params[[5]]))
			gtex_sig <- data.frame(get(as.character(params[[6]])), check.names = FALSE)
			tcga_sig <- data.frame(get(as.character(params[[7]])), check.names = FALSE)
			
			### get aracne info of the hub
			gtex_aracne <- gtex_aracne[[2]][[gtex_aracne[[1]][as.character(params[[1]]),3]]]
			tcga_aracne <- tcga_aracne[[2]][[tcga_aracne[[1]][as.character(params[[1]]),3]]]
			
			### get target gene info
			gtex_targets <- abs(gtex_aracne[,1])
			tcga_targets <- abs(tcga_aracne[,1])
			common_targets <- intersect(gtex_targets, tcga_targets)
			gtex_only_targets <- setdiff(gtex_targets, common_targets)
			tcga_only_targets <- setdiff(tcga_targets, common_targets)
			
			### make Excel data
			target_genes <- matrix(NA, max(length(gtex_targets), length(tcga_targets)), 10)
			colnames(target_genes) <- c("GTEx_Targets_Entrez", "GTEx_Targets_Symbol",
					"TCGA_Targets_Entrez", "TCGA_Targets_Symbol",
					"Common_Targets_Entrez", "Common_Targets_Symbol",
					"GTEx_Only_Targets_Entrez", "GTEx_Only_Targets_Symbol",
					"TCGA_Only_Targets_Entrez", "TCGA_Only_Targets_Symbol")
			
			target_genes[1:length(gtex_targets),1] <- gtex_targets
			target_genes[1:length(gtex_targets),2] <- entrezIDtoSymbol(gtex_targets)
			target_genes[1:length(tcga_targets),3] <- tcga_targets
			target_genes[1:length(tcga_targets),4] <- entrezIDtoSymbol(tcga_targets)
			if(length(common_targets) > 0) {
				target_genes[1:length(common_targets),5] <- common_targets
				target_genes[1:length(common_targets),6] <- entrezIDtoSymbol(common_targets)
			}
			if(length(gtex_only_targets > 0)) {
				target_genes[1:length(gtex_only_targets),7] <- gtex_only_targets
				target_genes[1:length(gtex_only_targets),8] <- entrezIDtoSymbol(gtex_only_targets)
			}
			if(length(tcga_only_targets > 0)) {
				target_genes[1:length(tcga_only_targets),9] <- tcga_only_targets
				target_genes[1:length(tcga_only_targets),10] <- entrezIDtoSymbol(tcga_only_targets)
			}
			
			
			### create output result directory
			dirName <- paste0(params[[2]], "_", entrezIDtoSymbol(as.character(params[[1]])))
			dir.create(paste0(params[[9]], dirName), showWarnings = FALSE)
			
			### print target gene list as Excel file
			write.xlsx2(target_genes, paste0(params[[9]], dirName, "/", dirName, "_target_genes.xlsx"), row.names = TRUE, col.names = TRUE)
			
			
			### pathway analysis
			# GTEx
			pathwayAnalysis_CP(de_thrsd = data.frame(Entrez_ID=gtex_targets), org = "human",
					lfcThreshold = NA, pThreshold = NA,
					title = paste("GTEx", params[[2]], params[[1]], "Target_Genes", sep = "_"),
					dir = paste0(params[[9]], dirName, "/"),
					suppl = FALSE)
			# TCGA
			pathwayAnalysis_CP(de_thrsd = data.frame(Entrez_ID=tcga_targets), org = "human",
					lfcThreshold = NA, pThreshold = NA,
					title = paste("TCGA", params[[2]], params[[1]], "Target_Genes", sep = "_"),
					dir = paste0(params[[9]], dirName, "/"),
					suppl = FALSE)
			# Exclusives
			pathwayAnalysis_CP(de_thrsd = data.frame(Entrez_ID=c(gtex_only_targets, tcga_only_targets)), org = "human",
					lfcThreshold = NA, pThreshold = NA,
					title = paste("Exclusive", params[[2]], params[[1]], "Target_Genes", sep = "_"),
					dir = paste0(params[[9]], dirName, "/"),
					suppl = FALSE)
			
			
			### heatmap
			### load library
			if(!require(gplots)) {
				install.packages("gplots")
				library(gplots)
			}
			
			### A function for scaling for heatmap
			scale_h <- function(data, type, na.rm=TRUE) {
				
				if(type == "row") {
					scaled <- t(scale(t(data)))
				} else if(type == "col") {
					scaled <- scale(data)
				} else {
					stop("Type is required: row or col")
				}
				
				if(na.rm == TRUE && (length(which(is.na(scaled))) > 0))  {
					scaled <- scaled[-unique(which(is.na(scaled), arr.ind = TRUE)[,1]),]
				}
				
				return(scaled)
			}
			
			### make heatmap data
			gtex_dat <- rbind(gtex_exp[as.character(common_targets),],
					gtex_exp[as.character(gtex_only_targets),],
					gtex_exp[as.character(tcga_only_targets),])
			rownames(gtex_dat) <- c(entrezIDtoSymbol(common_targets),
					entrezIDtoSymbol(gtex_only_targets),
					entrezIDtoSymbol(tcga_only_targets))
			tcga_dat <- rbind(tcga_exp[as.character(common_targets),],
					tcga_exp[as.character(gtex_only_targets),],
					tcga_exp[as.character(tcga_only_targets),])
			rownames(tcga_dat) <- c(entrezIDtoSymbol(common_targets),
					entrezIDtoSymbol(gtex_only_targets),
					entrezIDtoSymbol(tcga_only_targets))
			dat <- cbind(gtex_dat, tcga_dat)
			
			### scale the data
			dat <- scale_h(dat, type = "row", na.rm = FALSE)
			
			### set colside colors
			col_colors <- c(rep("gray", ncol(gtex_dat)), rep("black", ncol(tcga_dat)))
			names(col_colors) <- c(rep("GTEx", ncol(gtex_dat)), rep("TCGA", ncol(tcga_dat)))
			
			### set rowside colors
			row_colors <- c(rep("gold", length(common_targets)),
					rep("royalblue", length(gtex_only_targets)),
					rep("orange", length(tcga_only_targets)))
			names(row_colors) <- c(rep("Commons", length(common_targets)),
					rep("GTEx_Only", length(gtex_only_targets)),
					rep("TCGA_Only", length(tcga_only_targets)))
			
			### heatmap
			png(paste0(params[[9]], dirName, "/", dirName, "_target_genes_heatmap.png"),
					width = 1000, height = 1000)
			par(oma=c(0,0,2,6))
			breaks = seq(-5,5,length.out=101)
			heatmap.3(as.matrix(dat), main = paste0(dirName, "_Target_Genes"),
					xlab = "", ylab = "", col=greenred(100),
					scale="none", key=T, keysize=0.8, dendrogram = 'none', trace = 'none',
					labRow = rownames(dat), labCol = FALSE,
					Rowv = FALSE, Colv = FALSE, breaks = breaks,
					ColSideColors = cbind(as.vector(col_colors)),
					RowSideColors = t(cbind(as.vector(row_colors))),
					cexRow = 0.8, cexCol = 0.8, na.rm = TRUE)
			legend("topright", inset = 0.02, xpd = TRUE, title = "Sample Annotation", legend = unique(names(col_colors)), fill = unique(col_colors), cex = 1, box.lty = 0)
			legend("left", inset = -0.01, title = "Target Types", legend = unique(names(row_colors)), fill = unique(row_colors), cex = 1, box.lty = 0)
			dev.off()
			
			### heatmap with random targets for control
			common_genes <- intersect(rownames(gtex_exp), rownames(tcga_exp))
			set.seed(1234)
			random_genes <- common_genes[sample(length(common_genes), nrow(dat), replace = FALSE)]
			dat_random <- cbind(gtex_exp[random_genes,],
					tcga_exp[random_genes,])
			rownames(dat_random) <- entrezIDtoSymbol(random_genes)
			
			### scale the data
			dat_random <- scale_h(dat_random, type = "row", na.rm = FALSE)
			
			### heatmap - random targets
			png(paste0(params[[9]], dirName, "/", dirName, "_target_genes_heatmap_random.png"),
					width = 1000, height = 1000)
			par(oma=c(0,0,2,6))
			breaks = seq(-5,5,length.out=101)
			heatmap.3(as.matrix(dat_random), main = paste0(dirName, "_Random_Target_Genes"),
					xlab = "", ylab = "", col=greenred(100),
					scale="none", key=T, keysize=0.8, dendrogram = 'none', trace = 'none',
					labRow = rownames(dat_random), labCol = FALSE,
					Rowv = FALSE, Colv = FALSE, breaks = breaks,
					ColSideColors = cbind(as.vector(col_colors)),
					cexRow = 0.8, cexCol = 0.8, na.rm = TRUE)
			legend("topright", inset = 0.02, xpd = TRUE, title = "Sample Annotation", legend = unique(names(col_colors)), fill = unique(col_colors), cex = 1, box.lty = 0)
			dev.off()
			
			
			### network graph
			### need to save figures manually
			### load library
			if(!require(viper)) {
				source("https://bioconductor.org/biocLite.R")
				biocLite("viper")
				library(viper)
			}
			if(!require(RedeR)) {
				source("https://bioconductor.org/biocLite.R")
				biocLite("RedeR")
				library(RedeR)
			}
			if(!require(igraph)) {
				source("https://bioconductor.org/biocLite.R")
				biocLite("igraph")
				library(igraph)
			}
			
			### set edge threshold
			edgeThreshold <- as.numeric(params[[8]])
			
			### a function to remove nodes that have all the zero values
			removeZeroDegreeNodes<-function(mat) {
				zeroDegree <- 0
				cnt <- 1
				for(i in 1:nrow(mat)) {
					if((sum(mat[i,]) == 0) && (sum(mat[,i]) == 0)) {
						zeroDegree[cnt] <- i
						cnt <- cnt + 1
					}
				}
				
				return (mat[-zeroDegree, -zeroDegree])
			}
			
			### GTEx
			### correlation matrix using the target expressions
			c <- suppressWarnings(cor(t(gtex_dat), method = "spearman", use = "pairwise.complete.obs"))
			
			### diagonal -> 0
			diag(c) <- 0
			
			### abs and NA -> 0
			c <- abs(c)
			c[which(is.na(c), arr.ind = TRUE)] <- 0
			
			### filter edges with the threshold
			c[which(c < edgeThreshold, arr.ind = TRUE)] <- 0
			
			### remove nodes with all the zero values
			c <- removeZeroDegreeNodes(c)
			
			### exp and viper sig
			graph_exp <- apply(gtex_dat, 1, function(x) median(x, na.rm = TRUE))
			gtex_dat2 <- rbind(gtex_sig[as.character(common_targets),],
					gtex_sig[as.character(gtex_only_targets),],
					gtex_sig[as.character(tcga_only_targets),])
			rownames(gtex_dat2) <- rownames(gtex_dat)
			graph_sig <- apply(gtex_dat2, 1, function(x) median(x, na.rm = TRUE))
			
			### remove zero-degree nodes from exp and sig
			graph_exp <- graph_exp[rownames(c)]
			graph_sig <- graph_sig[rownames(c)]
			
			### igraph variables creation
			g <- graph.adjacency(c, mode = "undirected", weighted = TRUE)
			rescale <- function(x) (x-min(x))/(max(x) - min(x)) * 10
			E(g)$width <- rescale(E(g)$weight)
			E(g)$edgeColor <- "gray"
			g$legEdgeColor$scale <- "gray"
			V(g)$nodeSize <- (graph_exp / max(graph_exp, na.rm = TRUE)) * 30
			V(g)$size <- V(g)$nodeSize
			V(g)$legNodeSize <- V(g)$nodeSize
			V(g)$nodeLineColor <- "black"
			V(g)$VIPER_SIG <- graph_sig
			g <- att.setv(g, from = "VIPER_SIG", to = "nodeColor", breaks = seq(-3, 3, 0.5), pal = 2)
			
			### load RedeR screen and plot the graph
			rdp<-RedPort()
			calld(rdp)
			addGraph(rdp,g, layout.kamada.kawai(g))
			
			### add legends
			# color
			addLegend.color(rdp, g)
			# size
			circleLabel<-floor(seq(min(V(g)$nodeSize),max(V(g)$nodeSize),(max(V(g)$nodeSize) - min(V(g)$nodeSize))/4))
			circleSize<-(circleLabel / max(circleLabel)) * 30
			addLegend.size(rdp,sizevec=circleSize,labvec=circleLabel,title="Target_Expression")
			# title
			shape <- c("ELLIPSE", "DIAMOND")
			addLegend.shape(rdp,shape, title=paste0("GTEx_", dirName, "_Target_Genes Cor > ", edgeThreshold), position="topleft", ftsize=25, vertical=FALSE, dxtitle=50, dxborder=10, dyborder=-60)
			
			### TCGA
			### correlation matrix using the target expressions
			c <- suppressWarnings(cor(t(tcga_dat), method = "spearman", use = "pairwise.complete.obs"))
			
			### diagonal -> 0
			diag(c) <- 0
			
			### abs and NA -> 0
			c <- abs(c)
			c[which(is.na(c), arr.ind = TRUE)] <- 0
			
			### filter edges with the threshold
			c[which(c < edgeThreshold, arr.ind = TRUE)] <- 0
			
			### remove nodes with all the zero values
			c <- removeZeroDegreeNodes(c)
			
			### exp and viper sig
			graph_exp <- apply(tcga_dat, 1, function(x) median(x, na.rm = TRUE))
			tcga_dat2 <- rbind(tcga_sig[as.character(common_targets),],
					tcga_sig[as.character(gtex_only_targets),],
					tcga_sig[as.character(tcga_only_targets),])
			rownames(tcga_dat2) <- rownames(tcga_dat)
			graph_sig <- apply(tcga_dat2, 1, function(x) median(x, na.rm = TRUE))
			
			### remove zero-degree nodes from exp and sig
			graph_exp <- graph_exp[rownames(c)]
			graph_sig <- graph_sig[rownames(c)]
			
			### igraph variables creation
			g <- graph.adjacency(c, mode = "undirected", weighted = TRUE)
			rescale <- function(x) (x-min(x))/(max(x) - min(x)) * 10
			E(g)$width <- rescale(E(g)$weight)
			E(g)$edgeColor <- "gray"
			g$legEdgeColor$scale <- "gray"
			V(g)$nodeSize <- (graph_exp / max(graph_exp, na.rm = TRUE)) * 30
			V(g)$size <- V(g)$nodeSize
			V(g)$legNodeSize <- V(g)$nodeSize
			V(g)$nodeLineColor <- "black"
			V(g)$VIPER_SIG <- graph_sig
			g <- att.setv(g, from = "VIPER_SIG", to = "nodeColor", breaks = seq(-3, 3, 0.5), pal = 2)
			
			### load RedeR screen and plot the graph
			rdp<-RedPort()
			calld(rdp)
			addGraph(rdp,g, layout.kamada.kawai(g))
			
			### add legends
			# color
			addLegend.color(rdp, g)
			# size
			circleLabel<-floor(seq(min(V(g)$nodeSize),max(V(g)$nodeSize),(max(V(g)$nodeSize) - min(V(g)$nodeSize))/4))
			circleSize<-(circleLabel / max(circleLabel)) * 30
			addLegend.size(rdp,sizevec=circleSize,labvec=circleLabel,title="Target_Expression")
			# title
			shape <- c("ELLIPSE", "DIAMOND")
			addLegend.shape(rdp,shape, title=paste0("TCGA_", dirName, "_Target_Genes Cor > ", edgeThreshold), position="topleft", ftsize=25, vertical=FALSE, dxtitle=50, dxborder=10, dyborder=-60)
			
		} else {
			writeLines("required params do not exist")
		}
	}
	
	
	# ******************** which = viper_signature_all  *****************************
	# Generate Viper signatures using gene expression and Aracne network
	# Gene expression, Aracne network, and outputPath are required as "params"
	# Generates Viper signature matrix for all the tissues in params[[1]] and params[[2]]
	#
	# Gene expression: GTEx_36_EMat_Viper.rda or TCGA_28_EMat.rda
	# Aracne network: GTEx_36_Regulons.rda or TCGA_28_Regulons.rda
	#
	# There is one more option: number of samples for null model
	# Default value for the number is 15, if you don't specify it
	#
	# params[[1]]: Expression matrices variable name
	# params[[2]]: Aracne regulons variable name
	# params[[3]]: number of samples for random selection for reference (null) model
	# params[[4]]: method when making viper signature (should be either "zscore", "ttest", or "mean")
	# params[[5]]: permutation number for viper signature
	# params[[6]]: output path
	#
	# e.g., params <- list("emat_gtex_names", "gtexRegNames", 50, "zscore", 1000, "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/GTEx_36_ViperSigs.rda")
	# # e.g., params <- list("emat_tcga_names", "tcgaRegNames", 50, "zscore", 1000, "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/TCGA_28_ViperSigs.rda")
	
	if(which == "viper_signature_all"){
		if(!is.null(params) && length(params) > 5) {
			
			### load library
			if(!require(viper)) {
				source("https://bioconductor.org/biocLite.R")
				biocLite("viper")
				library(viper)
			}
			
			### a function to combine two matrices with different rows
			combineMats <- function(A, B) {
				if(is.null(dim(A))) {
					return(B)
				} else if(is.null(dim(B))) {
					return(A)
				} else {
					genes <- union(rownames(A), rownames(B))
					C <- data.frame(matrix(0, length(genes), (ncol(A)+ncol(B))))
					C <- cbind(data.frame(matrix(min(A), length(genes), ncol(A))), data.frame(matrix(min(B), length(genes), ncol(B))))
					rownames(C) <- genes
					colnames(C)[1:ncol(A)] <- colnames(A)
					colnames(C)[(ncol(A)+1):ncol(C)] <- colnames(B)
					C[rownames(A),1:ncol(A)] <- A
					C[rownames(B),(ncol(A)+1):ncol(C)] <- B
					return(C)
				}
			}
			
			### get emat and Aracne regulon names
			emat_variableNames <- get(params[[1]])
			reg_variableNames <- get(params[[2]])
			
			### make vipersig names
			assign("vipersig_names",
					sapply(strsplit(emat_variableNames, split = "_", fixed = TRUE),
							function(x) paste("vipersig", x[2], x[3], sep = "_")),
					envir = globalenv())
			
			### make temp directory
			dir.create("./temp/", showWarnings = FALSE)
			
			### iteratively makes vipersig for each tissue
			for(i in 1:length(emat_variableNames)) {
				### get expression
				exp <- get(emat_variableNames[i])
				
				### get Aracne regulon
				regulon <- get(reg_variableNames[i])
				
				### make a reference matrix
				N <- as.integer(params[[3]])
				set.seed(1234)
				refMat <- 0
				for(j in 1:length(emat_variableNames)) {
					if(emat_variableNames[j] != get(params[[1]])[i]) {
						tempMat <- get(emat_variableNames[j])
						refMat <- combineMats(refMat, tempMat[,sample(ncol(tempMat),N)])
					}
				}
				
				### make the same rows for the test matrix and the reference
				genes <- union(rownames(refMat), rownames(exp))
				
				# refMat
				refMat2 <- refMat
				refMat <- data.frame(matrix(min(refMat), length(genes), ncol(refMat2)))
				rownames(refMat) <- genes
				colnames(refMat) <- colnames(refMat2)
				refMat[rownames(refMat2),] <- refMat2
				
				# testMat
				testMat <- data.frame(matrix(min(exp), length(genes), ncol(exp)))
				rownames(testMat) <- genes
				colnames(testMat) <- colnames(exp)
				testMat[rownames(exp),] <- exp
				
				### viper
				sigMethod <- as.character(params[[4]])
				permutNum <- as.integer(params[[5]])
				vpsig <- viperSignature(eset = as.matrix(testMat), ref = as.matrix(refMat), method = sigMethod, per = permutNum, seed = 1, verbose = FALSE)
				
				### assign the signature
				assign(vipersig_names[i], vpsig$signature, envir = globalenv())
				
				### save the signature just in case
				save(list = c(vipersig_names[i]), file = paste0("./temp/", vipersig_names[i], ".rda"))
				
				### print the progress
				writeLines(paste(i, "/", length(emat_variableNames)))
			}
			
			### README function
			README <- function(){
				writeLines(paste(rep("#", 100), collapse = ""))
				writeLines("Each vipersig object has tissue-specific Viper signatures")
				writeLines("vipersig_names contains all the variable names of the Viper signature objects.")
				writeLines("For each tissue, 50 samples from each of the remaining tissues were randomly collected and used as referenece.")
				writeLines("When generating Viper signature, z-score was chosen.")
				writeLines(paste(rep("#", 100), collapse = ""))
			}
			
			### save the vipersig results as RDA file
			save(list = c(vipersig_names, "vipersig_names", "README"), file = as.character(params[[6]]))
			
		} else {
			writeLines("required params do not exist")
		}  
	}
	
	
	# ******************** which = systematic_hub_investigation  *****************************
	# Identifies hubs which have different Viper activity scores between GTEx and TCGA
	# There are 29 mapping info of GTEx vs TCGA
	# This will be done in the 29 comparisons
	# And based on Regulon info, target genes of the top hubs will be presented
	# The final results will be lists of differentially activated hubs,
	# Excel files of target genes of the hubs, Venn diagrams between GTEx and TCGA tissues,
	# and pathway analysis results of all the cases
	#
	# Needs all_64_aracne_mi_fixed.rda, all_64_vipermats.rda, all_64_vipersigs.rda,
	# and GTEx_TCGA_Map.rda loaded
	#
	# params[[1]]: Mapping info of GTEx vs TCGA tissues
	# params[[2]]: The number of the top hubs which will be investigated
	# params[[3]]: minimum MI threshold
	# params[[4]]: maximum p-value threshold
	# params[[5]]: measure to compare VIPER NES ["mean" or "ttest"]
	# params[[6]]: output directory
	#
	# e.g., params <- list("GTEx_TCGA_Map", 50, 0.1, 0.0001, "ttest", "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/viper/systematic_analysis/")
	
	if(which == "systematic_hub_investigation"){
		if(!is.null(params) && length(params) > 2) {
			
			### load library
			if(!require(VennDiagram)) {
				source("https://bioconductor.org/biocLite.R")
				biocLite("VennDiagram")
				library(VennDiagram)
			}
			if(!require(gridExtra)) {
				install.packages("gridExtra")
				library(gridExtra)
			}
			
			### get the mapping info
			mapping_info <- get(as.character(params[[1]]))
			
			### make empty pathway lists
			gtex_pathway <- list()
			tcga_pathway <- list()
			
			### iteratively perform the investigation for every comparison
			for(i in 1:nrow(mapping_info)) {
				
				### create output result directory
				dirName <- paste0(mapping_info[i,1], "_", mapping_info[i,2])
				dir.create(paste0(params[[6]], dirName), showWarnings = FALSE)
				
				### get all the data ready
				gtex_aracne <- get(as.character(mapping_info[i,1]))
				tcga_aracne <- get(paste("tcga", mapping_info[i,2], sep = "_"))
				gtex_viper <- get(paste0("vmat_gtex_", mapping_info[i,1]))
				tcga_viper <- get(paste0("vmat_tcga_", mapping_info[i,2]))
				gtex_viperSig <- get(paste0("vipersig_gtex_", mapping_info[i,1]))
				tcga_viperSig <- get(paste0("vipersig_tcga_", mapping_info[i,2]))
				
				### get common hubs
				common_hubs <- intersect(rownames(gtex_viper), rownames(tcga_viper))
				common_hubs_symbols <- entrezIDtoSymbol(common_hubs)
				
				### only have info of the common hubs
				gtex_viper <- gtex_viper[common_hubs,]
				tcga_viper <- tcga_viper[common_hubs,]
				
				if(as.character(params[[5]] == "mean")) {
					### VIPER NES mean difference of each hub
					mean_diff <- apply(gtex_viper, 1, mean) - apply(tcga_viper, 1, mean)
					
					### absoulute VIPER NES difference
					mean_diff_abs <- abs(mean_diff)
					mean_diff_abs <- mean_diff_abs[order(-mean_diff_abs)]
					
					### get top diff hubs
					top_diff_hubs <- names(mean_diff_abs)[1:as.integer(params[[2]])]
					names(top_diff_hubs) <- entrezIDtoSymbol(top_diff_hubs)
				} else if(as.character(params[[5]] == "ttest")) {
					### t-test
					t_diff <- t(sapply(1:length(common_hubs), function(x) {
										t <- t.test(gtex_viper[x,], tcga_viper[x,])
										return(c(t$statistic, t$p.value))
									}))
					rownames(t_diff) <- common_hubs
					
					### order t_diff based on p-value (ascending order)
					t_diff <- t_diff[order(t_diff[,2]),]
					
					### get top diff hubs based on t-test
					top_diff_hubs <- rownames(t_diff)[1:as.integer(params[[2]])]
					names(top_diff_hubs) <- entrezIDtoSymbol(top_diff_hubs)
				} else {
					stop("Wrong parameter [[5]]: it should be either mean or ttest")
				}
				
				### get target genes of the top diff hubs
				gtex_target_genes <- lapply(top_diff_hubs, function(x) {
							temp <- gtex_aracne[[2]][[gtex_aracne[[1]][x,2]]]
							return(rownames(temp[intersect(which(temp[,"MI"] > as.numeric(params[[3]])),
															which(temp[,"Pvalue"] < as.numeric(params[[4]]))),]))
						})
				tcga_target_genes <- lapply(top_diff_hubs, function(x) {
							temp <- tcga_aracne[[2]][[tcga_aracne[[1]][x,2]]]
							return(rownames(temp[intersect(which(temp[,"MI"] > as.numeric(params[[3]])),
															which(temp[,"Pvalue"] < as.numeric(params[[4]]))),]))
						})
				
				### order the target genes (common between GTEx and TCGA first)
				common <- lapply(1:length(gtex_target_genes), function(x) {
							return(intersect(gtex_target_genes[[x]], tcga_target_genes[[x]]))
						})
				gtex_target_genes <- lapply(1:length(gtex_target_genes), function(x) {
							return(c(common[[x]], setdiff(gtex_target_genes[[x]], common[[x]])))
						})
				tcga_target_genes <- lapply(1:length(tcga_target_genes), function(x) {
							return(c(common[[x]], setdiff(tcga_target_genes[[x]], common[[x]])))
						})
				
				### make the target gene format
				maxLen <- max(c(sapply(gtex_target_genes, length), sapply(tcga_target_genes, length)))
				target_genes <- matrix("", maxLen, (length(gtex_target_genes)*3))
				colnames(target_genes) <- rep("", ncol(target_genes))
				for(j in 1:ncol(target_genes)) {
					if(j %% 3 == 1) {
						if(length(gtex_target_genes[[j/3+1]]) > 0) {
							target_genes[1:length(gtex_target_genes[[j/3+1]]),j] <- entrezIDtoSymbol(gtex_target_genes[[j/3+1]])
							colnames(target_genes)[j] <- paste0("GTEx_Hub_", names(top_diff_hubs)[j/3+1])
						}
					} else if(j %% 3 == 2) {
						if(length(tcga_target_genes[[j/3+1]]) > 0) {
							target_genes[1:length(tcga_target_genes[[j/3+1]]),j] <- entrezIDtoSymbol(tcga_target_genes[[j/3+1]])
							colnames(target_genes)[j] <- paste0("TCGA_Hub_", names(top_diff_hubs)[j/3+1])
						}
					}
				}
				
				### save the target genes of the hubs
				write.table(target_genes, file = paste0(params[[6]], dirName, "/",
								params[[2]], "_top_diff_hubs_target_genes_mi_",
								params[[3]], "_pv_", params[[4]], "_",
								params[[5]], ".txt"),
						sep = "\t", row.names = FALSE)
				
				### add number of targets and p-value info to the top diff hubs
				# 2 x 2 contingency table for fisher's exact test
				# the test is to measure significance of relation of a hub between GTEx and TCGA
				x <- sapply(common, length)
				y <- sapply(gtex_target_genes, length) - x
				z <- sapply(tcga_target_genes, length) - x
				w <- total_geneNum
				
				### calculate p-values
				pv <- NULL
				for(j in 1:length(x)) {
					pv <- c(pv, fisher.test(matrix(c(x[j], z[j], y[j], w), 2, 2))$p.value)
				}
				
				if(as.character(params[[5]] == "mean")) {
					### make top diff hub result - data.frame + target #s & p-values
					hub_result <- cbind(TOP_DIFF_HUBS=names(top_diff_hubs),
							ENTREZ_ID=top_diff_hubs,
							MEAN_DIFF=mean_diff[1:as.integer(params[[2]])],
							GTEX_TARGET_NUM=(x+y),
							TCGA_TARGET_NUM=(x+z),
							COMMON_TARGET_NUM=x,
							INTERSECT_P_VALUE=pv)
				} else if(as.character(params[[5]] == "ttest")) {
					### make top diff hub result - data.frame + target #s & p-values
					hub_result <- cbind(TOP_DIFF_HUBS=names(top_diff_hubs),
							ENTREZ_ID=top_diff_hubs,
							T_STATISTIC=t_diff[1:as.integer(params[[2]]),1],
							DIFF_P_VALUE=t_diff[1:as.integer(params[[2]]),2],
							GTEX_TARGET_NUM=(x+y),
							TCGA_TARGET_NUM=(x+z),
							COMMON_TARGET_NUM=x,
							INTERSECT_P_VALUE=pv)
				}
				
				### save the top diff hubs
				write.table(hub_result, file = paste0(params[[6]], dirName, "/",
								params[[2]], "_top_diff_hubs_mi_",
								params[[3]], "_pv_", params[[4]], "_",
								params[[5]], ".txt"),
						sep = "\t", row.names = FALSE)
				
				### iteratively perform venn diagram and pathway analysis for each hub
				gtex_pathway[[i]] <- list()
				tcga_pathway[[i]] <- list()
				for(j in 1:length(top_diff_hubs)) {
					### create directory for each hub
					dir.create(paste0(params[[6]], dirName, "/", names(top_diff_hubs)[j]), showWarnings = FALSE)
					
					### venn diagram
					if((length(gtex_target_genes[[j]]) > 0) && (length(tcga_target_genes[[j]]) > 0)) {
						v <- venn.diagram(list(gtex_target_genes[[j]], tcga_target_genes[[j]]),
								category.names = c("GTEx", "TCGA"),
								cat.cex = 1.5, cex = 1.5,
								filename = NULL)
						
						### save the diagram as png
						png(paste0(params[[6]], dirName, "/", names(top_diff_hubs)[j], "/",
										names(top_diff_hubs)[j], "_Targets_Venn_diagram_mi_",
										params[[3]], "_pv_", params[[4]], ".png"),
								width = 1200, height = 1000, res = 150)
						grid.arrange(gTree(children=v),
								top=paste0(names(top_diff_hubs)[j], " Targets Venn Diagram"),
								bottom="")
						dev.off()
					} else if((length(gtex_target_genes[[j]]) > 0) && (length(tcga_target_genes[[j]]) == 0)) {
						v <- venn.diagram(list(gtex_target_genes[[j]]),
								category.names = c("GTEx"),
								cat.cex = 1.5, cex = 1.5,
								filename = NULL)
						
						### save the diagram as png
						png(paste0(params[[6]], dirName, "/", names(top_diff_hubs)[j], "/",
										names(top_diff_hubs)[j], "_Targets_Venn_diagram_mi_",
										params[[3]], "_pv_", params[[4]], ".png"),
								width = 1200, height = 1000, res = 150)
						grid.arrange(gTree(children=v),
								top=paste0(names(top_diff_hubs)[j], " Targets Venn Diagram"),
								bottom="")
						dev.off()
					} else if((length(gtex_target_genes[[j]]) == 0) && (length(tcga_target_genes[[j]]) > 0)) {
						v <- venn.diagram(list(tcga_target_genes[[j]]),
								category.names = c("TCGA"),
								cat.cex = 1.5, cex = 1.5,
								filename = NULL)
						
						### save the diagram as png
						png(paste0(params[[6]], dirName, "/", names(top_diff_hubs)[j], "/",
										names(top_diff_hubs)[j], "_Targets_Venn_diagram_mi_",
										params[[3]], "_pv_", params[[4]], ".png"),
								width = 1200, height = 1000, res = 150)
						grid.arrange(gTree(children=v),
								top=paste0(names(top_diff_hubs)[j], " Targets Venn Diagram"),
								bottom="")
						dev.off()
					} else {
						writeLines("No Venn diagram result")
					}
					
					### pathway analysis
					# GTEx
					gtex_pathway[[i]][[j]] <- pathwayAnalysis_CP(geneList = gtex_target_genes[[j]],
							org = "human",
							database = "GO",
							displayNum = 50,
							title = paste0("GTEx_", names(top_diff_hubs)[j], "_Targets_Pathways_mi_",
									params[[3]], "_pv_", params[[4]]),
							dir = paste0(params[[6]], dirName, "/", names(top_diff_hubs)[j], "/"))
					# TCGA
					tcga_pathway[[i]][[j]] <- pathwayAnalysis_CP(geneList = tcga_target_genes[[j]],
							org = "human",
							database = "GO",
							displayNum = 50,
							title = paste0("TCGA_", names(top_diff_hubs)[j], "_Targets_Pathways_mi_",
									params[[3]], "_pv_", params[[4]]),
							dir = paste0(params[[6]], dirName, "/", names(top_diff_hubs)[j], "/"))
				}
				
				### heatmap of changed NESs
				### load library
				if(!require(gplots)) {
					install.packages("gplots")
					library(gplots)
				}
				
				### A function for scaling for heatmap
				scale_h <- function(data, type, na.rm=TRUE) {
					
					if(type == "row") {
						scaled <- t(scale(t(data)))
					} else if(type == "col") {
						scaled <- scale(data)
					} else {
						stop("Type is required: row or col")
					}
					
					if(na.rm == TRUE && (length(which(is.na(scaled))) > 0))  {
						scaled <- scaled[-unique(which(is.na(scaled), arr.ind = TRUE)[,1]),]
					}
					
					return(scaled)
				}
				
				### make heatmap data
				gtex_dat <- lapply(common, function(x) gtex_viperSig[as.character(x),])
				tcga_dat <- lapply(common, function(x) tcga_viperSig[as.character(x),])
				
				### iteratively perform heatmap for each hub
				for(j in 1:length(gtex_dat)) {
					if(length(common[[j]]) > 0) {
						### cbind gtex and tcga data
						dat <- cbind(rbind(gtex_dat[[j]]), rbind(tcga_dat[[j]]))
						
						### change rownames to symbols from entrez id
						rownames(dat) <- entrezIDtoSymbol(common[[j]])
						
						### scale the data
						dat <- scale_h(dat, type = "row", na.rm = FALSE)
						
						### set colside colors
						col_colors <- c(rep("gray", ncol(rbind(gtex_dat[[j]]))), rep("black", ncol(rbind(tcga_dat[[j]]))))
						names(col_colors) <- c(rep("GTEx", ncol(rbind(gtex_dat[[j]]))), rep("TCGA", ncol(rbind(tcga_dat[[j]]))))
						
						### heatmap
						png(paste0(params[[6]], dirName, "/", names(top_diff_hubs)[j], "/",
										dirName, "_", names(top_diff_hubs)[j], "_common_genes_heatmap_mi_",
										params[[3]], "_pv_", params[[4]], ".png"),
								width = 1000, height = 1000)
						par(oma=c(0,0,2,6))
						if(nrow(as.matrix(dat)) == 1) {
							dat <- rbind(dat, dat)
							row_colors <- rep("blue", 2)
							names(row_colors) <- rownames(dat)
							heatmap.3(as.matrix(dat), main = paste0(dirName, "_", names(top_diff_hubs)[j], "_Common_Genes_ViperSigs"),
									xlab = "", ylab = "", col=greenred(100),
									scale="none", key=T, keysize=1, dendrogram = 'none', trace = 'none',
									labRow = FALSE, labCol = FALSE,
									Rowv = FALSE, Colv = FALSE,
									ColSideColors = cbind(as.vector(col_colors)),
									RowSideColors = t(cbind(as.vector(row_colors))),
									cexRow = 1.5, cexCol = 1.5, na.rm = TRUE)
							legend("left", inset = 0.01, title = "Target Name", legend = unique(names(row_colors)), fill = unique(row_colors), cex = 1, box.lty = 0)
						} else {
							heatmap.3(as.matrix(dat), main = paste0(dirName, "_", names(top_diff_hubs)[j], "_Common_Genes_ViperSigs"),
									xlab = "", ylab = "", col=greenred(100),
									scale="none", key=T, keysize=1, dendrogram = 'none', trace = 'none',
									labRow = rownames(dat), labCol = FALSE,
									Rowv = FALSE, Colv = FALSE,
									ColSideColors = cbind(as.vector(col_colors)),
									cexRow = 1.5, cexCol = 1.5, na.rm = TRUE)
						}
						legend("topright", inset = 0.02, xpd = TRUE, title = "Sample Annotation", legend = unique(names(col_colors)), fill = unique(col_colors), cex = 1.5, box.lty = 0)
						dev.off()
						
						### random heatmap
						set.seed(1234)
						common_targets <- intersect(rownames(gtex_viperSig), rownames(tcga_viperSig))
						random_targets <- common_targets[sample(length(common_targets), 10)]
						dat <- cbind(rbind(gtex_viperSig[random_targets,]),
								rbind(tcga_viperSig[random_targets,]))
						
						### change rownames to symbols from entrez id
						rownames(dat) <- entrezIDtoSymbol(random_targets)
						
						### scale the data
						dat <- scale_h(dat, type = "row", na.rm = FALSE)
						
						### heatmap
						png(paste0(params[[6]], dirName, "/", names(top_diff_hubs)[j], "/",
										dirName, "_", names(top_diff_hubs)[j], "_common_genes_heatmap_mi_",
										params[[3]], "_pv_", params[[4]], "_random.png"),
								width = 1000, height = 1000)
						par(oma=c(0,0,2,6))
						if(nrow(as.matrix(dat)) == 1) {
							dat <- rbind(dat, dat)
							row_colors <- rep("blue", 2)
							names(row_colors) <- rownames(dat)
							heatmap.3(as.matrix(dat), main = paste0(dirName, "_", names(top_diff_hubs)[j], "_Random_Genes_ViperSigs"),
									xlab = "", ylab = "", col=greenred(100),
									scale="none", key=T, keysize=1, dendrogram = 'none', trace = 'none',
									labRow = FALSE, labCol = FALSE,
									Rowv = FALSE, Colv = FALSE,
									ColSideColors = cbind(as.vector(col_colors)),
									RowSideColors = t(cbind(as.vector(row_colors))),
									cexRow = 1.5, cexCol = 1.5, na.rm = TRUE)
							legend("left", inset = 0.01, title = "Target Name", legend = unique(names(row_colors)), fill = unique(row_colors), cex = 1, box.lty = 0)
						} else {
							heatmap.3(as.matrix(dat), main = paste0(dirName, "_", names(top_diff_hubs)[j], "_Random_Genes_ViperSigs"),
									xlab = "", ylab = "", col=greenred(100),
									scale="none", key=T, keysize=1, dendrogram = 'none', trace = 'none',
									labRow = rownames(dat), labCol = FALSE,
									Rowv = FALSE, Colv = FALSE,
									ColSideColors = cbind(as.vector(col_colors)),
									cexRow = 1.5, cexCol = 1.5, na.rm = TRUE)
						}
						legend("topright", inset = 0.02, xpd = TRUE, title = "Sample Annotation", legend = unique(names(col_colors)), fill = unique(col_colors), cex = 1.5, box.lty = 0)
						dev.off()
					}
				}
			}
			
			### load hub info of every comparison
			total_hub_info <- lapply(1:nrow(mapping_info), function(i) {
						read.table(file = paste0(params[[6]], paste0(mapping_info[i,1], "_", mapping_info[i,2]), "/",
										params[[2]], "_top_diff_hubs_mi_",
										params[[3]], "_pv_", params[[4]], "_",
										params[[5]], ".txt"),
								header = TRUE, sep = "\t")
					})
			
			### make an empty union hub vector
			union_hub_names <- NULL
			for(i in 1:length(total_hub_info)) {
				union_hub_names <- union(union_hub_names, total_hub_info[[i]][,1])
			}
			union_hubs <- rep(0, length(union_hub_names))
			names(union_hubs) <- union_hub_names
			union_hub_comparison <- rep("", length(union_hub_names))
			names(union_hub_comparison) <- union_hub_names
			
			### count the number of hubs
			for(i in 1:length(total_hub_info)) {
				union_hubs[as.character(total_hub_info[[i]][,1])] <- union_hubs[as.character(total_hub_info[[i]][,1])]+1
				union_hub_comparison[as.character(total_hub_info[[i]][,1])] <- paste0(union_hub_comparison[as.character(total_hub_info[[i]][,1])], ",", mapping_info[i,1], "_", mapping_info[i,2])
			}
			union_hub_comparison <- sapply(union_hub_comparison, function(x) substr(x, 2, nchar(x)))
			
			### order the union hubs based on the count
			union_hub_comparison <- union_hub_comparison[order(-union_hubs)]
			union_hubs <- union_hubs[order(-union_hubs)]
			
			### add p-value
			total_hubs <- NULL
			for(i in 1:length(varNamesVP)) {
				total_hubs <- union(total_hubs, rownames(get(varNamesVP[i])))
			}
			union_hub_pv <- NULL
			for(i in 1:length(union_hubs)) {
				union_hub_pv[i] <- fisher.test(matrix(c(as.integer(params[[2]]),
										length(total_hubs)-as.integer(params[[2]]),
										union_hubs[i],
										length(total_hub_info)-union_hubs[i]), ncol=2))$p.value
			}
			
			### save the union hubs and the counts
			union_hub_result <- data.frame(cbind(Hub_Name=names(union_hubs),
							Count=union_hubs,
							Comparison_Name=union_hub_comparison,
							P_Value=union_hub_pv))
			union_hub_result$Count <- paste(union_hub_result$Count, "out of", length(total_hub_info))
			write.table(union_hub_result,
					file = paste0(params[[6]],
							params[[2]], "_top_diff_hubs_",
							params[[5]], "_counts.txt"),
					sep = "\t", row.names = FALSE)
			
			
			### load target info
			total_target_info <- lapply(1:nrow(mapping_info), function(i) {
						buf <- read.table(file = paste0(params[[6]], paste0(mapping_info[i,1], "_", mapping_info[i,2]), "/",
										params[[2]], "_top_diff_hubs_target_genes_mi_",
										params[[3]], "_pv_", params[[4]], "_",
										params[[5]], ".txt"),
								header = TRUE, sep = "\t")
						return(buf[,-3*(1:as.numeric(params[[2]]))])
					})
			
			### iteratively creates union targets and the counts
			for(i in 1:nrow(mapping_info)) {
				### empty union target names
				gtex_union_target_names <- NULL
				tcga_union_target_names <- NULL
				
				### get unique target names
				for(j in 2*(1:as.numeric(params[[2]]))-1) {
					gtex_union_target_names <- union(gtex_union_target_names, total_target_info[[i]][,j])
				}
				gtex_union_target_names <- gtex_union_target_names[-which(gtex_union_target_names == "")]
				for(j in 2*(1:as.numeric(params[[2]]))) {
					tcga_union_target_names <- union(tcga_union_target_names, total_target_info[[i]][,j])
				}
				tcga_union_target_names <- tcga_union_target_names[-which(tcga_union_target_names == "")]
				
				### remove NAs
				if(length(which(is.na(gtex_union_target_names))) > 0) {
					gtex_union_target_names <- gtex_union_target_names[-which(is.na(gtex_union_target_names))]
				}
				if(length(which(is.na(tcga_union_target_names))) > 0) {
					tcga_union_target_names <- tcga_union_target_names[-which(is.na(tcga_union_target_names))]
				}
				
				### set zero for all the union targets
				gtex_union_targets <- rep(0, length(gtex_union_target_names))
				names(gtex_union_targets) <- gtex_union_target_names
				tcga_union_targets <- rep(0, length(tcga_union_target_names))
				names(tcga_union_targets) <- tcga_union_target_names
				
				gtex_union_target_hubs <- rep("", length(gtex_union_target_names))
				names(gtex_union_target_hubs) <- gtex_union_target_names
				tcga_union_target_hubs <- rep("", length(tcga_union_target_names))
				names(tcga_union_target_hubs) <- tcga_union_target_names
				
				### count the number of targets
				for(j in 1:length(gtex_union_targets)) {
					gtex_union_targets[j] <- length(which(total_target_info[[i]][,2*(1:as.numeric(params[[2]]))-1] == gtex_union_target_names[j]))
					gtex_union_target_hubs[j] <- paste(as.character(total_hub_info[[i]]$TOP_DIFF_HUBS[which(total_target_info[[i]][,2*(1:as.numeric(params[[2]]))-1] == gtex_union_target_names[j], arr.ind = TRUE)[,2]]), collapse = ",")
				}
				for(j in 1:length(tcga_union_targets)) {
					tcga_union_targets[j] <- length(which(total_target_info[[i]][,2*(1:as.numeric(params[[2]]))] == tcga_union_target_names[j]))
					tcga_union_target_hubs[j] <- paste(as.character(total_hub_info[[i]]$TOP_DIFF_HUBS[which(total_target_info[[i]][,2*(1:as.numeric(params[[2]]))] == tcga_union_target_names[j], arr.ind = TRUE)[,2]]), collapse = ",")
				}
				
				### order the union targets based on the count
				gtex_union_target_hubs <- gtex_union_target_hubs[order(-gtex_union_targets)]
				gtex_union_targets <- gtex_union_targets[order(-gtex_union_targets)]
				tcga_union_target_hubs <- tcga_union_target_hubs[order(-tcga_union_targets)]
				tcga_union_targets <- tcga_union_targets[order(-tcga_union_targets)]
				
				### add p-value
				gtex_union_target_pv <- NULL
				gtex_avg_targetNum <- round(sum(gtex_union_targets) / as.integer(params[[2]]))
				for(j in 1:length(gtex_union_targets)) {
					gtex_union_target_pv[j] <- fisher.test(matrix(c(gtex_avg_targetNum,
											total_geneNum-gtex_avg_targetNum,
											gtex_union_targets[j],
											as.integer(params[[2]])-gtex_union_targets[j]), ncol=2))$p.value
				}
				tcga_union_target_pv <- NULL
				tcga_avg_targetNum <- round(sum(tcga_union_targets) / as.integer(params[[2]]))
				for(j in 1:length(tcga_union_targets)) {
					tcga_union_target_pv[j] <- fisher.test(matrix(c(tcga_avg_targetNum,
											total_geneNum-tcga_avg_targetNum,
											tcga_union_targets[j],
											as.integer(params[[2]])-tcga_union_targets[j]), ncol=2))$p.value
				}
				
				### save the target info
				gtex_union_target_result <- data.frame(cbind(Target_Name=names(gtex_union_targets),
								Count=gtex_union_targets,
								Hub_Info=gtex_union_target_hubs,
								P_Value=gtex_union_target_pv))
				gtex_union_target_result$Count <- paste(gtex_union_target_result$Count, "out of", as.numeric(params[[2]]))
				write.table(gtex_union_target_result,
						file = paste0(params[[6]], paste0(mapping_info[i,1], "_", mapping_info[i,2]), "/",
								params[[2]], "_top_diff_hubs_mi_",
								params[[3]], "_pv_", params[[4]], "_",
								params[[5]], "_gtex_target_counts.txt"),
						sep = "\t", row.names = FALSE)
				tcga_union_target_result <- data.frame(cbind(Target_Name=names(tcga_union_targets),
								Count=tcga_union_targets,
								Hub_Info=tcga_union_target_hubs,
								P_Value=tcga_union_target_pv))
				tcga_union_target_result$Count <- paste(tcga_union_target_result$Count, "/", as.numeric(params[[2]]))
				write.table(tcga_union_target_result,
						file = paste0(params[[6]], paste0(mapping_info[i,1], "_", mapping_info[i,2]), "/",
								params[[2]], "_top_diff_hubs_mi_",
								params[[3]], "_pv_", params[[4]], "_",
								params[[5]], "_tcga_target_counts.txt"),
						sep = "\t", row.names = FALSE)
			}
			
			
			### R object creation
			Systematic_Hub_Investigation_Results <- list()
			### info of all the hubs appeared in all the 29 comparisons
			Systematic_Hub_Investigation_Results[[1]] <- read.table(file = paste0(params[[6]],
							params[[2]], "_top_diff_hubs_",
							params[[5]], "_counts.txt"),
					sep = "\t", header = TRUE,
					row.names = 1)
			### load target genes of the hubs
			total_target_info <- lapply(1:nrow(mapping_info), function(i) {
						buf <- read.table(file = paste0(params[[6]], paste0(mapping_info[i,1], "_", mapping_info[i,2]), "/",
										params[[2]], "_top_diff_hubs_target_genes_mi_",
										params[[3]], "_pv_", params[[4]], "_",
										params[[5]], ".txt"),
								header = TRUE, sep = "\t")
						return(buf[,-3*(1:as.numeric(params[[2]]))])
					})
			Systematic_Hub_Investigation_Results[[2]] <- list()
			### iteratively add info of each comparison
			for(i in 1:nrow(mapping_info)) {
				Systematic_Hub_Investigation_Results[[2]][[i]] <- list()
				### info of the top 50 differentially activated hubs
				Systematic_Hub_Investigation_Results[[2]][[i]][[1]] <- read.table(file = paste0(params[[6]], paste0(mapping_info[i,1], "_", mapping_info[i,2]), "/",
								params[[2]], "_top_diff_hubs_mi_",
								params[[3]], "_pv_", params[[4]], "_",
								params[[5]], ".txt"),
						sep = "\t", header = TRUE,
						row.names = 1)
				### info of all the target genes from GTEx tissue
				Systematic_Hub_Investigation_Results[[2]][[i]][[2]] <- read.table(file = paste0(params[[6]], paste0(mapping_info[i,1], "_", mapping_info[i,2]), "/",
								params[[2]], "_top_diff_hubs_mi_",
								params[[3]], "_pv_", params[[4]], "_",
								params[[5]], "_gtex_target_counts.txt"),
						sep = "\t", header = TRUE,
						row.names = 1)
				### info of all the target genes from TCGA tissue
				Systematic_Hub_Investigation_Results[[2]][[i]][[3]] <- read.table(file = paste0(params[[6]], paste0(mapping_info[i,1], "_", mapping_info[i,2]), "/",
								params[[2]], "_top_diff_hubs_mi_",
								params[[3]], "_pv_", params[[4]], "_",
								params[[5]], "_tcga_target_counts.txt"),
						sep = "\t", header = TRUE,
						row.names = 1)
				Systematic_Hub_Investigation_Results[[2]][[i]][[4]] <- list()
				Systematic_Hub_Investigation_Results[[2]][[i]][[5]] <- list()
				Systematic_Hub_Investigation_Results[[2]][[i]][[6]] <- list()
				Systematic_Hub_Investigation_Results[[2]][[i]][[7]] <- list()
				### iteratively add target genes and pathway results
				for(j in 1:as.integer(params[[2]])) {
					### clean
					if((nrow(gtex_pathway[[i]][[j]]) == 0) || is.null(nrow(gtex_pathway[[i]][[j]]))) {
						gtex_pathway[[i]][[j]] <- NA
					}
					if((nrow(tcga_pathway[[i]][[j]]) == 0) || is.null(nrow(tcga_pathway[[i]][[j]]))) {
						tcga_pathway[[i]][[j]] <- NA
					}
					
					### GTEx
					Systematic_Hub_Investigation_Results[[2]][[i]][[4]][[j]] <- as.vector(total_target_info[[i]][,2*j-1])
					if(length(which(Systematic_Hub_Investigation_Results[[2]][[i]][[4]][[j]] == "")) > 0) {
						Systematic_Hub_Investigation_Results[[2]][[i]][[4]][[j]] <- Systematic_Hub_Investigation_Results[[2]][[i]][[4]][[j]][-which(Systematic_Hub_Investigation_Results[[2]][[i]][[4]][[j]] == "")]
					}
					Systematic_Hub_Investigation_Results[[2]][[i]][[6]][[j]] <- gtex_pathway[[i]][[j]]
					### TCGA
					Systematic_Hub_Investigation_Results[[2]][[i]][[5]][[j]] <- total_target_info[[i]][,2*j]
					if(length(which(Systematic_Hub_Investigation_Results[[2]][[i]][[5]][[j]] == "")) > 0) {
						Systematic_Hub_Investigation_Results[[2]][[i]][[5]][[j]] <- Systematic_Hub_Investigation_Results[[2]][[i]][[5]][[j]][-which(Systematic_Hub_Investigation_Results[[2]][[i]][[5]][[j]] == "")]
					}
					Systematic_Hub_Investigation_Results[[2]][[i]][[7]][[j]] <- tcga_pathway[[i]][[j]]
				}
				### name the hubs
				names(Systematic_Hub_Investigation_Results[[2]][[i]][[4]]) <- rownames(Systematic_Hub_Investigation_Results[[2]][[i]][[1]])
				names(Systematic_Hub_Investigation_Results[[2]][[i]][[5]]) <- rownames(Systematic_Hub_Investigation_Results[[2]][[i]][[1]])
				names(Systematic_Hub_Investigation_Results[[2]][[i]][[6]]) <- rownames(Systematic_Hub_Investigation_Results[[2]][[i]][[1]])
				names(Systematic_Hub_Investigation_Results[[2]][[i]][[7]]) <- rownames(Systematic_Hub_Investigation_Results[[2]][[i]][[1]])
			}
			### name the comparisons
			names(Systematic_Hub_Investigation_Results[[2]]) <- paste0(mapping_info[,1], "_", mapping_info[,2])
			
			### README function
			README = function() {
				writeLines(paste(rep("#", 100), collapse = ""))
				writeLines("The \"Systematic_Hub_Investigation_Results\" contains many hub and target info")
				writeLines(paste(rep("-", 100), collapse = ""))
				writeLines("Systematic_Hub_Investigation_Results[[1]]: It contains info of all the top 50 differentially activated hubs")
				writeLines("from all the 29 comparisons.")
				writeLines("The columns are Hub_Name, Count, Comparison_Name, and P_Value.")
				writeLines("The Hub_Name is gene symbol of a hub.")
				writeLines("The Count indicates how many times the hub is appeared among 29 comparisons.")
				writeLines("The Comparison_Name contains comparison info where the hub is appeared.")
				writeLines("Lastly, the P_Value has p-values from Fisher's exact test based on")
				writeLines("how significant the ratio (the count / the total number of top hubs) is when compared to the background ratio.")
				writeLines(paste(rep("-", 100), collapse = ""))
				writeLines("Systematic_Hub_Investigation_Results[[2]] is a list and has length of 29.")
				writeLines("length(Systematic_Hub_Investigation_Results[[2]]) == 29")
				writeLines("Each list represents one of the 29 mappings between GTEx and TCGA")
				writeLines("See names(Systematic_Hub_Investigation_Results[[2]])")
				writeLines(paste(rep("-", 100), collapse = ""))
				writeLines("Systematic_Hub_Investigation_Results[[2]][[x]][[1]] - x varies on 1:29")
				writeLines("It contains info of the top differentially activated hubs in the certain comparison")
				writeLines("There are 8 columns: TOP_DIFF_HUBS, ENTREZ_ID, T_STATISTIC, DIFF_P_VALUE,")
				writeLines("GTEX_TARGET_NUM, TCGA_TARGET_NUM, COMMON_TARGET_NUM, and INTERSECT_P_VALUE.")
				writeLines("The TOP_DIFF_HUBS means the hub name, and the ENTREZ_ID represents Entrez id of the hub.")
				writeLines("The T_STATISTIC is from t-test and the DIFF_P_VALUE contains p-values from the t-test.")
				writeLines("The GTEX_TARGET_NUM means the number of target genes of the hub in GTEx tissue")
				writeLines("while the TCGA_TARGET_NUM means that in TCGA tissue.")
				writeLines("The COMMON_TARGET_NUM has the number of target genes shared between GTEx and TCGA of the hub.")
				writeLines("The INTERSECT_P_VALUE indicates p-values of Fisher's exact test on the number of common target genes.")
				writeLines(paste(rep("-", 100), collapse = ""))
				writeLines("Systematic_Hub_Investigation_Results[[2]][[x]][[2]] - x varies on 1:29")
				writeLines("It contians info that how many times GTEx targets were appeared among the top hubs")
				writeLines("The first column contains target genes and the second column means")
				writeLines("how many times those targets genes were appeared among the top hubs in the comparison.")
				writeLines("(Hubs have their target genes connected based on MI)")
				writeLines("The third column means which top hubs have the target gene, and the fourth – p-values.")
				writeLines(paste(rep("-", 100), collapse = ""))
				writeLines("Systematic_Hub_Investigation_Results[[2]][[x]][[3]] - x varies on 1:29")
				writeLines("It contians info that how many times TCGA targets were appeared among the top hubs")
				writeLines("The first column contains target genes and the second column means")
				writeLines("how many times those targets genes were appeared among the top hubs in the comparison.")
				writeLines("(Hubs have their target genes connected based on MI)")
				writeLines("The third column means which top hubs have the target gene, and the fourth – p-values.")
				writeLines(paste(rep("-", 100), collapse = ""))
				writeLines(paste0("Systematic_Hub_Investigation_Results[[2]][[x]][[4]] is a list and has length of 50."))
				writeLines(paste0("length(Systematic_Hub_Investigation_Results[[2]][[x]][[4]]) == 50"))
				writeLines(paste0("Each list represents one of the top 50 hubs"))
				writeLines("which were differentially activated between GTEx and TCGA")
				writeLines("Each list has a vector of the hub's GTEx target genes")
				writeLines(paste(rep("-", 100), collapse = ""))
				writeLines(paste0("Systematic_Hub_Investigation_Results[[2]][[x]][[5]] is a list and has length of 50."))
				writeLines(paste0("length(Systematic_Hub_Investigation_Results[[2]][[x]][[5]]) == 50"))
				writeLines(paste0("Each list represents one of the top 50 hubs"))
				writeLines("which were differentially activated between GTEx and TCGA")
				writeLines("Each list has a vector of the hub's TCGA target genes")
				writeLines(paste(rep("#", 100), collapse = ""))
			}
			
			### save the R object
			save(list = c("Systematic_Hub_Investigation_Results", "README"),
					file = paste0(params[[6]],
							params[[2]], "_top_diff_hubs_",
							params[[5]], "_info.rda"))
			
		} else {
			writeLines("required params do not exist")
		}
	}
	
	# ******************** which = regulon_conservation  *****************************
	# For each hub gene X, it identifies the top interactome pairs where the regulon
	# of X is most highly cosnserved (using the FET p-values from thPairEnrich). It then
	# counts and reports how many of these pairs involve exlusively GTEx interactomes or
	# exclusively TCGA interactomes or one each from GTEx and TCGA. See README function 
	# below for details of how these counts are computed.
	#
	# ARGUMENTS:
	# The params[[]] list is expected to contain the following arguments:
	# * params[[1]]:	A vector with hub genes identifiers, in the form of character 
	#		gene symbols or integer/character Entrez IDs.
	# * params[[2]]:	The number of "top" pairwise interactome comparisons to use for
	#		calculating the counts (a value of 50 or 100 should be OK).
	# * params[[3]]: 	File name prefix to use for the results file, a ".rda" file containing
	#		the results object, and a ".csv" file containing a spreadsheet version of the 
	#		results object.
	# 
	# RETURN VALUE:
	# The metmod call returns the object "reg_conv". See the text in the README method below 
	# for a description of its contents.
	if(which == "regulon_conservation"){
		# README function, to be stored in the .rda output file
		README <- function(){
			writeLines("The object \"reg_conv\" contains the results of the analysis performed using")
			writeLines("the following method call from aracne.R:")
			writeLines("\toneOffs(which = \"regulon_conservation\", params = list(hubs_all, 100, \"regulon_conservation\"))")
			writeLines("where \"hubs_all\" is a vector containing the Entrez IDs of all hub genes seen in at")
			writeLines("least one of the GTEx and TCGA interactomes. The results object \"reg_conv\" is a list with")
			writeLines("one entry for each hub gene in hubs_all (and names(reg_con) = hubs_all). The entry")
			writeLines("corresponding to hub gene X provides data that quantify how well the regulon of X")
			writeLines("is conserved across GTEx and TCGA interactomes and further clarify if the regulons ")
			writeLines("of X in TCGA networks are generally different than those in GTEx networks, thus")
			writeLines("allowing the possibility to discover systematic regulatory rewiring. To generate")
			writeLines("the data in reg_conv[[X]] (which are described below) we use fhe following process.")
			writeLines("* STEP 1: Initialize counters named \"gtex\", \"tcga\", and \"both\" to the value of 0.")
			writeLines("* STEP 2: Use the object tfPairEnrich to extract the FET log(P-value) scores for all ")
			writeLines("  the N = choose(length(varNames), 2) pairwise interactome comparison, representing ")
			writeLines("  the conservation of the regulon of X across all N possible GTEx and TCGA itteractome. ")
			writeLines("  pairs. Order these N scores from the most significant to the least and keep the top")
			writeLines("  100 (this number comes from params[[2]], in the call above).")
			writeLines("* STEP 3: For each of these 100 scores, consider the two networks A, B in the ")
			writeLines("  interactome pair comparison. If A and B are both GTEx interactomes, increase")
			writeLines("  the counter \"gtex\" by one. If A and B are both TCGA interactomes, increase")
			writeLines("  the counter \"tcga\" by one. If A and B are one each GTEx and TCGA, increase ")
			writeLines("  the counter \"both\" by one.")
			writeLines("")
			writeLines("reg_conv[[X]][[1]] is the counter vector c(\"gtex\", \"both\", \"tcga\").")
			writeLines("reg_conv[[X]][[2]] is another length-3 vector containing the average FET log(pvalues).")
			writeLines("E.g.,  reg_conv[[X]][[2]][\"GTEx\"] contains the average FET log(pvalue) of all the ")
			writeLines("pairwise comparisons contributing counts to reg_conv[[X]][[1]][\"GTEx\"].")
		}
		
		
		# Save the results object to a CSV file, with one row per hub gene. Each row will contain
		# both the count and the average FET lot(pvalues).
		saveToCSV <- function(rc, csv_file){
			mat = t(sapply(rc, function(x){return(c(x[[1]], x[[2]]))}))
			sym = entrezIDtoSymbol(rownames(mat))
			sym[is.na(sym)] = "NA"
			mat = cbind(sym, rownames(mat), mat)
			colnames(mat) = c("Hub_gene_symbol", "Hub_gene_id", "GTEx_counts", "Both_counts", "TCGA_counts", "GTEx_FET", "Both_FET", "TCGA_FET")
			write.csv(mat, csv_file, row.names=FALSE)
		}
		
		if (is.null(params) | length(params) != 3)
			stop("Wrong number of parameters passed.")
		hubs = params[[1]]
		top = params[[2]]
		filePrefix = params[[3]]
		reg_conv = List()
		k = 1
		for (hub in hubs){
			fets = sort(getRegulonFETs(hub, replaceInf = TRUE))
			fets = fets[1:min(top, length(fets))]
			fets = fets[fets < 0]	# Don't look at meaningless entries
			hub_counts = avg_fets = rep(0, 3)
			names(hub_counts) = names(avg_fets) = c("GTEx", "Both", "TCGA")
			if (length(fets) > 0)
				for (i in 1:length(fets)){
					#writeLines(names(fets)[i])
					ind = length(grep(pattern="tcga", strsplit(names(fets)[i], "%")[[1]]))
					hub_counts[ind+1] = hub_counts[ind+1] + 1
					avg_fets[ind+1] = avg_fets[ind+1] + fets[i]
				}
			for (i in 1:3)
				if (avg_fets[i] != 0)
					avg_fets[i] = avg_fets[i] / hub_counts[i]
			reg_conv[[k]] = list(hub_counts, avg_fets)
			k = k + 1
		}
		names(reg_conv) = hubs
		save(reg_conv, README, file=paste(filePrefix, "rda", sep="."))
		saveToCSV(reg_conv, paste(filePrefix, "csv", sep="."))
		return(reg_conv)
	}
	
	# ******************** which = MakeChrRDA  *****************************
	# In Aracne network, each hub has its own targets
	# Store the location of targets on chromosomes in RDA
	# The ChrInfo is a list with length of the number of tissues
	# And each tissue list has list in it with length of the number of hubs in a tissue
	# Each hub list has a matrix that describes chromosome numbers of its targets
	# chrInfo[[tissue]][[hub]][target]
	# Needs All_64_Aracne_MI_Fixed.rda
	#
	# params[[1]]: a character vector of variable names of aracne networks
	# params[[2]]: the output RDA file path
	#
	# e.g., params <- list("varNames", "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/all_64_target_chrs.rda")
	
	if(which == "MakeChrRDA"){
		if(!is.null(params) && length(params) > 1) {
			
			### get variable names
			netNames <- get(as.character(params[[1]]))
			
			### iteratively getting chromosome info for each tissue
			chrInfo <- list()
			for(i in 1:length(netNames)) {
				### get a network
				net <- get(netNames[i])
				
				### get target names
				targets <- lapply(rownames(net[[1]]), function(x) {
							return(rownames(net[[2]][[net[[1]][as.character(x),2]]]))
						})
				names(targets) <- rownames(net[[1]])
				
				### get chromosome info of the targets
				chrInfo[[i]] <- list()
				for(j in 1:length(targets)) {
					chrInfo[[i]][[j]] <- getChromosome(targets[[j]])
				}
				names(chrInfo[[i]]) <- names(targets)
				
				### remove NAs
				for(j in 1:length(targets)) {
					if(length(which(is.na(chrInfo[[i]][[j]]))) > 0) {
						chrInfo[[i]][[j]][-which(is.na(chrInfo[[i]][[j]]))]
					}
				}
				
				### print progress
				writeLines(paste(i, "/", length(netNames)))
			}
			names(chrInfo) <- netNames
			
			### set README function
			README <- function(){
				writeLines(paste(rep("#", 100), collapse = ""))
				writeLines("The \"chrInfo\" contains chromosome info of target genes of every hub of every tissue.")
				writeLines("chrInfo has double list that the first one indicates tissue")
				writeLines("and the second one indicates hub.")
				writeLines("chrInfo[[tissue]][[hub]][target]")
				writeLines("For example, if you want to know \"AdiposeSub\" (tissue name) tissue's")
				writeLines("\"285855\" (hub Entrez ID) hub's \"23\" (target Entrez ID) target's chromosome number, ")
				writeLines("it is saved in chrInfo[[\"AdiposeSub\"]][[\"285855\"]][\"23\"].")
				writeLines(paste(rep("#", 100), collapse = ""))
			}
			
			### save as RDA
			save(list = c("chrInfo", "README"), file = as.character(params[[2]]))
			
		} else {
			writeLines("required params do not exist")
		}  
	}
	
	# ******************** which = hub_target_chromosomes  *****************************
	# The hub has its target genes and we would like to know which chromosome contains
	# those target genes
	# And also, statistical significance test (prop.test() and Fisher's exact test)
	# will be performed to show which chromosomes are statistically enriched with
	# target genes
	# The results will be a text file that has chromosome info and their statistics
	# Needs All_64_Aracne_MI_Fixed.rda and Gene_num_chr.rda
	# params[[1]]: a character vector of variable names of aracne networks
	# params[[2]]: a vector of input hub names in Entrez ID format
	# params[[3]]: "Gene_num_chr.rda" path - used for normalization
	# params[[4]]: target MI threshold
	# params[[5]]: target p-value threshold
	# params[[6]]: output file directory
	# e.g., params=list("varNames", 1081, "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/Gene_num_chr.rda", 0, 1, "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/chromosome_analysis/")
	
	if(which == "hub_target_chromosomes"){
		if(!is.null(params) && length(params) > 5) {
			
			### load library
			if(!require(scales)) {
				install.packages("scales")
				library(scales)
			}
			
			### get Aracne variable names
			netNames <- get(as.character(params[[1]]))
			
			### make an empty chromosome list
			chrs <- list()
			
			### iteratively getting chromosome info of targets in each tissue
			for(i in 1:length(netNames)) {
				### get a network
				net <- get(netNames[i])
				
				if(length(which(rownames(net[[1]]) == as.character(params[[2]]))) > 0) {
					### get target genes of the input hub
					temp <- net[[2]][[net[[1]][as.character(params[[2]]),2]]]
					temp <- temp[intersect(which(temp[,"MI"] >= as.numeric(params[[4]])),
									which(temp[,"Pvalue"] <= as.numeric(params[[5]]))),,drop=F]
					targets <- rownames(temp)
					
					### get chromosome info of the targets
					chrs[[i]] <- getChromosome(targets)
				} else {
					chrs[[i]] <- NA
				}
			}
			
			### load gene numbers of chromosomes
			load(as.character(params[[3]]))
			total_geneNum_chr <- sum(geneNum)
			
			### make an empty data frame for the result (chromosome significance)
			result <- matrix(NA, length(netNames)*length(geneNum), 7)
			colnames(result) <- c("Tissue_Name", "Chromosome", "Target_Num", "Normalized_Target_Num", "Scaled_Target_Num", "Prop_PV", "Fisher_PV")
			result <- as.data.frame(result)
			result$Target_Num <- 0
			result$Normalized_Target_Num <- 0
			result$Scaled_Target_Num <- 0
			
			### make an empty data frame for the result2 (distribution significance)
			result2 <- as.data.frame(matrix(NA, length(netNames), 26))
			colnames(result2)[1] <- "Tissue_Name"
			colnames(result2)[2:24] <- paste0(rep("chr", 23), 1:23)
			colnames(result2)[25:26] <- c("Total_Target_Num", "Chi-Square_PV")
			result2[,2:25] <- 0
			
			### put designated values
			for(i in 1:length(netNames)) {
				### for the result
				for(j in 1:length(geneNum)) {
					if(startsWith(netNames[i], "tcga")) {
						result$Tissue_Name[(i-1)*length(geneNum)+j] <- netNames[i]
					} else {
						result$Tissue_Name[(i-1)*length(geneNum)+j] <- paste0("gtex_", netNames[i])
					}
					
					result$Chromosome[(i-1)*length(geneNum)+j] <- j
				}
				
				### for the result2
				if(startsWith(netNames[i], "tcga")) {
					result2$Tissue_Name[i] <- netNames[i]
				} else {
					result2$Tissue_Name[i] <- paste0("gtex_", netNames[i])
				}
			}
			
			### get the number of targets
			for(i in 1:length(netNames)) {
				for(j in 1:length(chrs[[i]])) {
					if(!is.na(chrs[[i]][j])) {
						result$Target_Num[(i-1)*length(geneNum)+chrs[[i]][j]] <- result$Target_Num[(i-1)*length(geneNum)+chrs[[i]][j]] + 1
						result2[i,chrs[[i]][j]+1] <- result2[i,chrs[[i]][j]+1] + 1
					}
				}
				result2$Total_Target_Num[i] <- sum(result2[i,2:24])
			}
			
			### normalize the target numbers with "geneNum"
			result$Normalized_Target_Num <- result$Target_Num / geneNum[result$Chromosome]
			
			### scale the target numbers to 1-10
			unique_tissue_name <- unique(result$Tissue_Name)
			for(i in 1:length(unique_tissue_name)) {
				if((length(chrs[[i]]) > 1) || (!is.na(chrs[[i]]))) {
					idx <- which(result$Tissue_Name == unique_tissue_name[i])
					result$Scaled_Target_Num[idx] <- round(rescale(result$Normalized_Target_Num[idx], to=c(1,10)))
				}
			}
			
			### prepare the ratio of the number of genes located on every chromosome divided by the total number of genes across all chromosomes
			geneNum_ratio <- geneNum / total_geneNum_chr
			
			### get p-values
			for(i in 1:length(netNames)) {
				if((length(chrs[[i]]) > 1) || (!is.na(chrs[[i]]))) {
					target_num_sum <- sum(result$Target_Num[which(result$Tissue_Name == unique_tissue_name[i])])
					for(j in 1:length(geneNum)) {
						X <- geneNum[j]
						Y <- total_geneNum_chr
						Z <- result$Target_Num[(i-1)*length(geneNum)+j]
						W <- target_num_sum
						
						### prop.test()
						suppressWarnings(result$Prop_PV[(i-1)*length(geneNum)+j] <-
										prop.test(c(Z*10,
														X),
												c(W*10,
														Y))$p.value)
						
						### Fisher's exact test
						suppressWarnings(result$Fisher_PV[(i-1)*length(geneNum)+j] <-
										fisher.test(matrix(c(X, Y-X, Z, W-Z), ncol = 2))$p.value)
					}
					
					### chi-square test
					suppressWarnings(result2$`Chi-Square_PV`[i] <-
									chisq.test(x = result2[i,2:24], p = geneNum_ratio)$p.value)
				}
			}
			
			### write out the results
			write.table(result, file = paste0(params[[6]], entrezIDtoSymbol(params[[2]]), "_",
							params[[2]], "_chromosome_significance.txt"),
					sep = "\t", row.names = FALSE)
			write.table(result2, file = paste0(params[[6]], entrezIDtoSymbol(params[[2]]), "_",
							params[[2]], "_distribution_significance.txt"),
					sep = "\t", row.names = FALSE)
		} else {
			writeLines("required params do not exist")
		}  
	}
	
  
  # ***************************** which = hub_target_chromosomes ********************************
  # Make a beeswarm plot with input hub
  # The hub has its target genes and we would like to know which chromosome contains
  # those target genes
  # The plot would have target genes' chromosome info of all the 64 (GTEx + TCGA) tissues
  # Needs All_64_Aracne_MI_Fixed.rda
  # params[[1]]: a character vector of variable names of aracne networks
  # params[[2]]: a vector of input hub names in Entrez ID format
  # params[[3]]: target MI threshold
  # params[[4]]: target p-value threshold
  # params[[5]]: Normalization by chromosome length or by number of genes ["chrLen" or "geneNum"]
  # params[[6]]: output file path
  # params[[7]]: "Gene_num_chr.rda" path - if params[[5]] == geneNum otherwise, NULL
  # e.g., params=list("varNames", 1081, 0, 1, "geneNum", "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/chromosome_analysis/1081_target_chrs_norm.pdf", "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/Gene_num_chr.rda")
  
  if (which == "hub_target_chromosomes"){
    if(!is.null(params) && length(params) > 4) {
      
      ### load library
      if(!require(ggplot2)) {
        install.packages("ggplot2")
        library(ggplot2)
      }
      if(!require(ggbeeswarm)) {
        install.packages("ggbeeswarm")
        library(ggbeeswarm)
      }
      if(!require(ggpubr)) {
        install.packages("ggpubr")
        library(ggpubr)
      }
      
      ### get variable names
      netNames <- get(as.character(params[[1]]))
      
      ### make an empty chromosome list
      chrs <- list()
      
      ### iteratively getting chromosome info of targets in each tissue
      for(i in 1:length(netNames)) {
        ### get a network
        net <- get(netNames[i])
        
        if(length(which(rownames(net[[1]]) == as.character(params[[2]]))) > 0) {
          ### get target genes of the input hub
          temp <- net[[2]][[net[[1]][as.character(params[[2]]),2]]]
          temp <- temp[intersect(which(temp[,"MI"] >= as.numeric(params[[3]])),
                                 which(temp[,"Pvalue"] <= as.numeric(params[[4]]))),,drop=F]
          targets <- rownames(temp)
          
          ### get chromosome info of the targets
          chrs[[i]] <- getChromosome(targets)
        } else {
          chrs[[i]] <- NA
        }
      }
      
      ### make a data for a beeswarm plot
      df <- matrix(NA, sum(sapply(chrs, length)), 2)
      cnt <- 1
      for(i in 1:length(netNames)) {
        for(j in 1:length(chrs[[i]])) {
          if(startsWith(netNames[i], "tcga")) {
            df[cnt, 1] <- netNames[i]
          } else {
            df[cnt, 1] <- paste0("gtex_", netNames[i])
          }
          df[cnt, 2] <- chrs[[i]][j]
          cnt <- cnt + 1
        }
      }
      colnames(df) <- c("group", "chromosome")
      df <- data.frame(df)
      
      ### remove NAs
      if(length(which(is.na(df$chromosome))) > 0) {
        df <- df[-which(is.na(df$chromosome)),]
      }
      
      ### characterize and numerized the factor columns
      df[1] <- lapply(df[1], function(x) as.character(x))
      df[2] <- lapply(df[2], function(x) as.numeric(as.character(x)))
      
      ### order df
      df <- df[order(df$group, df$chromosome),]
      
      ### normalize the frequency by chromosome length
      if(params[[5]] == "chrLen") {
        
        ### load library
        if(!require(BSgenome.Hsapiens.UCSC.hg19)) {
          source("https://bioconductor.org/biocLite.R")
          biocLite("BSgenome.Hsapiens.UCSC.hg19")
          library(BSgenome.Hsapiens.UCSC.hg19)
        }
        if(!require(scales)) {
          install.packages("scales")
          library(scales)
        }
        
        ### get chr length
        chr.lengths <- seqlengths(Hsapiens)[1:23]
        chr.lengths[23] <- chr.lengths[23] + seqlengths(Hsapiens)[24]
        names(chr.lengths) <- 1:23
        
        ### get the frequencies
        unique_location <- which(!duplicated(df))
        frequency <- unique_location[2:length(unique_location)] - unique_location[1:(length(unique_location)-1)]
        frequency[length(unique_location)] <- nrow(df) - unique_location[length(unique_location)] + 1
        
        ### make frequency-combined data
        new_df <- df[unique_location,]
        new_df$frequency <- frequency
        new_df$length <- chr.lengths[as.character(new_df$chromosome)]
        
        ### normalize by chromosome length
        new_df$new_freq <- new_df$frequency / new_df$length
        
        ### scale the new frequency
        new_df$scaled_freq <- 0
        unique_group <- unique(new_df$group)
        for(i in 1:length(unique_group)) {
          idx <- which(new_df$group == unique_group[i])
          new_df$scaled_freq[idx] <- round(rescale(new_df$new_freq[idx], to=c(1,10)))
        }
        
        ### make data frame with new scaled frequency
        df <- matrix(NA, sum(new_df$scaled_freq), 2)
        cnt <- 1
        for(i in 1:nrow(new_df)) {
          for(j in 1:new_df$scaled_freq[i]) {
            df[cnt,1] <- new_df$group[i]
            df[cnt,2] <- new_df$chromosome[i]
            cnt <- cnt + 1
          }
        }
        colnames(df) <- c("group", "chromosome")
        df <- data.frame(df)
        
        ### characterize and numerized the factor columns
        df[1] <- lapply(df[1], function(x) as.character(x))
        df[2] <- lapply(df[2], function(x) as.numeric(as.character(x)))
        
        ### order df
        df <- df[order(df$group, df$chromosome),]
        
      } else if(params[[5]] == "geneNum") {
        
        ### load library
        if(!require(scales)) {
          install.packages("scales")
          library(scales)
        }
        
        ### load gene numbers of all the chromosomes
        load(params[[7]])
        
        ### get the frequencies
        unique_location <- which(!duplicated(df))
        frequency <- unique_location[2:length(unique_location)] - unique_location[1:(length(unique_location)-1)]
        frequency[length(unique_location)] <- nrow(df) - unique_location[length(unique_location)] + 1
        
        ### make frequency-combined data
        new_df <- df[unique_location,]
        new_df$frequency <- frequency
        new_df$geneNum <- geneNum[new_df$chromosome]
        
        ### normalize by gene numbers
        new_df$new_freq <- new_df$frequency / new_df$geneNum
        
        ### scale the new frequency
        new_df$scaled_freq <- 0
        unique_group <- unique(new_df$group)
        for(i in 1:length(unique_group)) {
          idx <- which(new_df$group == unique_group[i])
          new_df$scaled_freq[idx] <- round(rescale(new_df$new_freq[idx], to=c(1,10)))
        }
        
        ### make data frame with new scaled frequency
        df <- matrix(NA, sum(new_df$scaled_freq), 2)
        cnt <- 1
        for(i in 1:nrow(new_df)) {
          for(j in 1:new_df$scaled_freq[i]) {
            df[cnt,1] <- new_df$group[i]
            df[cnt,2] <- new_df$chromosome[i]
            cnt <- cnt + 1
          }
        }
        colnames(df) <- c("group", "chromosome")
        df <- data.frame(df)
        
        ### characterize and numerized the factor columns
        df[1] <- lapply(df[1], function(x) as.character(x))
        df[2] <- lapply(df[2], function(x) as.numeric(as.character(x)))
        
        ### order df
        df <- df[order(df$group, df$chromosome),]
        
      }
      
      ### get unique tissue names
      unique_tissue_names <- unique(df[,1])
      
      ### make a plot
      plot_rowNum <- 8
      beeswarmNum <- ceiling(length(unique_tissue_names) / plot_rowNum)
      p <- list()
      for(i in 1:plot_rowNum) {
        p[[i]] <- ggplot(df[which(df[,1] %in% unique_tissue_names[((i-1)*beeswarmNum+1):(i*beeswarmNum)]),],
                         aes(x=group, y=chromosome)) +
          labs(y="Chromosome") +
          theme_classic(base_size = 16) +
          geom_boxplot() +
          geom_beeswarm(aes(color=group)) +
          stat_compare_means()
      }
      
      ### save the plot
      pdf(as.character(params[[6]]), width = 20, height = 10)
      for(i in 1:plot_rowNum) {
        print(p[[i]])
      }
      dev.off()
    } else {
      stop("required params do not exist")
    }
  }
  
	
	# ******************** which = all_hub_target_chromosomes  *****************************
	# The hub has its target genes and we would like to know which chromosome contains
	# those target genes
	# And also, statistical significance test (Fisher's exact test)
	# will be performed to show which chromosomes are statistically enriched with
	# target genes
	# This will be done to all the existing hubs
	# The results will be text files that have chromosome statistics for every hub and for every tissue
	# Needs All_62_ARACNE.rda and All_62_chr_distribution.rda
	# params[[1]]: a character vector of variable names of aracne networks
	# params[[2]]: "all_64_chr_distribution.rda" path - used as background info or for normalization
	# params[[3]]: target MI threshold
	# params[[4]]: target p-value threshold
	# params[[5]]: output file directory
	# e.g., params=list("varNames", "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/All_62_chr_distribution.rda", 0, 1, "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/chromosome_analysis/")
	# e.g., params=list("varNames", "./All_62_chr_distribution.rda", 0, 1, "./results/chromosome/")
	
	if(which == "all_hub_target_chromosomes"){
		if(!is.null(params) && length(params) > 4) {
			### get Aracne variable names
			netNames <- get(as.character(params[[1]]))
			
			### load gene numbers of chromosomes
			load(as.character(params[[2]]))
			
			### get all the union set of hubs
			union_hub_names <- NULL
			for(i in 1:length(netNames)) {
				### get a network
				net <- get(netNames[i])
				
				### union all
				union_hub_names <- union(union_hub_names, rownames(net[[1]]))
			}
			
			### make a chromosome map
			all_genes <- getInteractomeGenes(netNames, count = FALSE)
			chromosome_map <- getChromosome(all_genes)
			
			### make a chromosome count data
			chr_count <- list()
			for(i in 1:length(netNames)) {
				### get a network
				net <- get(netNames[i])
				
				### make an empty matrix
				chr_count[[i]] <- matrix(0, length(union_hub_names), 23)
				rownames(chr_count[[i]]) <- union_hub_names
				colnames(chr_count[[i]]) <- paste0(rep("chr", 23), 1:23)
				
				### count the number of target genes for every chromosome
				for(j in 1:length(union_hub_names)) {
					if(length(which(rownames(net[[1]]) == as.character(union_hub_names[j]))) > 0) {
						### get target genes of a given hub
						temp <- net[[2]][[net[[1]][as.character(union_hub_names[j]),2]]]
						temp <- temp[intersect(which(temp[,"MI"] >= as.numeric(params[[3]])),
										which(temp[,"Pvalue"] <= as.numeric(params[[4]]))),,drop=F]
						targets <- rownames(temp)
						
						### +1 for the appropriate chromosome
						for(k in 1:length(targets)) {
							chrIdx <- chromosome_map[targets[k]]
							chr_count[[i]][j,chrIdx] <- chr_count[[i]][j,chrIdx] + 1
						}
					}
				}
				
				### name the list
				if(startsWith(netNames[i], "tcga")) {
					names(chr_count)[i] <- netNames[i]
				} else {
					names(chr_count)[i] <- paste0("gtex_", netNames[i])
				}
				
				### print progress
				writeLines(paste(i, "/", length(netNames)))
			}
			
			### a function to check if a given vector has all zero values
			isVectorAllZero <- function(v) {
				r = TRUE
				for(i in 1:length(v)) {
					if(v[i] != 0) {
						r = FALSE
						break
					}
				}
				return(r)
			}
			
			### calculate the distribution significance
			dist_sig <- matrix(NA, length(union_hub_names), length(netNames))
			rownames(dist_sig) <- union_hub_names
			colnames(dist_sig) <- netNames
			
			### chi-square test
			for(i in 1:length(netNames)) {
				### prepare the ratio of the number of genes located on every chromosome divided by the total number of genes across all chromosomes
				geneNum <- geneChrInfo[[i]]
				total_geneNum_chr <- sum(geneNum)
				geneNum_ratio <- geneNum / total_geneNum_chr
				
				for(j in 1:length(union_hub_names)) {
					if(!isVectorAllZero(chr_count[[i]][j,])) {
						suppressWarnings(dist_sig[j,i] <- chisq.test(x = chr_count[[i]][j,], p = geneNum_ratio)$p.value)
					}
				}
			}
			
			### add gene symbol column to the distribution significance matrix
			dist_sig <- cbind(Entrez_ID=rownames(dist_sig),
					Gene_Symbol=entrezIDtoSymbol(rownames(dist_sig)),
					dist_sig)
			
			### write out the distribution significance matrix
			write.table(dist_sig, file = paste0(params[[5]], "distribution_significance.txt"),
					sep = "\t", row.names = FALSE)
			
			### set variable for the RDA file
			dist_sig <- apply(dist_sig[,-c(1,2)], 2, as.numeric)
			rownames(dist_sig) <- union_hub_names
			colnames(dist_sig) <- netNames
			assign("chr_enrich_all62", dist_sig, envir = globalenv())
			
			### calculate the chromosome significance
			for(i in 1:length(netNames)) {
				### make an empty result matrix
				chr_sig <- matrix(NA, length(union_hub_names), ncol(chr_count[[i]]))
				rownames(chr_sig) <- union_hub_names
				colnames(chr_sig) <- colnames(chr_count[[i]])
				
				### prepare the ratio of the number of genes located on every chromosome divided by the total number of genes across all chromosomes
				geneNum <- geneChrInfo[[i]]
				total_geneNum_chr <- sum(geneNum)
				geneNum_ratio <- geneNum / total_geneNum_chr
				
				### Fisher's exact test
				for(j in 1:nrow(chr_sig)) {
					total_targetNum <- sum(chr_count[[i]][j,])
					for(k in 1:ncol(chr_sig)) {
						X <- geneNum[k]
						Y <- total_geneNum_chr
						Z <- chr_count[[i]][j,k]
						W <- total_targetNum
						
						suppressWarnings(chr_sig[j,k] <- fisher.test(matrix(c(X, Y-X, Z, W-Z), ncol = 2))$p.value)
					}
				}
				
				### add gene symbol column to the chromosome significance matrix
				chr_sig <- cbind(Entrez_ID=rownames(chr_sig),
						Gene_Symbol=entrezIDtoSymbol(rownames(chr_sig)),
						chr_sig)
				
				### write out the chromosome significance matrix
				write.table(chr_sig, file = paste0(params[[5]], names(chr_count)[i], "_chromosome_significance.txt"),
						sep = "\t", row.names = FALSE)
				
				### set variable for the RDA file
				chr_sig <- apply(chr_sig[,-c(1,2)], 2, as.numeric)
				rownames(chr_sig) <- union_hub_names
				colnames(chr_sig) <- colnames(chr_count[[i]])
				assign(paste0("chr_", names(chr_count)[i]), chr_sig, envir = globalenv())
				
				### print progress
				writeLines(paste(i, "/", length(netNames)))
			}
			
			### remove hubs which are not in the corresponding Aracne network
			for(i in 1:length(netNames)) {
				### get chr info object
				temp <- get(paste0("chr_", names(chr_count)[i]))
				
				### get Aracne object
				temp2 <- get(netNames[i])
				
				### keep hubs only in the Aracne network
				temp <- temp[rownames(temp2[[1]]),]
				
				### overwrite
				assign(paste0("chr_", names(chr_count)[i]), temp, envir = globalenv())
			}
			
			### make a vector of variable names for the RDA file
			chr_varNames <- c("chr_enrich_all62", paste0("chr_", names(chr_count)))
			
			### make a README function for the RDA file
			README <- function(){
				writeLines(paste(rep("#", 100), collapse = ""))
				writeLines("The \"chr_varNames\" has all the variable names of the chromosome info.")
				writeLines("The first one (chr_varNames[1]) has results of the Chi-square test.")
				writeLines("From the second to the last ones (chr_varNames[2:63]), they have results of the tissue-specific Fisher's exact test.")
				writeLines("The \"chr_enrich_all62\" has a 6033 (hubs) x 62 (tissues) matrix.")
				writeLines("A value in a cell in the matrix indicates the p-value of the Chi-square test of the corresponding hub (row) and the corresponding tissue (column).")
				writeLines("The \"chr_TISSUE_NAME\" has a 6033 (hubs) x 23 (chromosomes) matrix. ")
				writeLines("A value in a cell in the matrix indicates the p-value of the Fisher's exact test of the corresponding tissue (TISSUE_NAME in an object), the corresponding hub (row), and the corresponding chromosome (column).")
				writeLines("The Chi-square test p-values indicate how significant the distribution of target genes on the chromosomes of a hub is different from the background.")
				writeLines("And the Fisher's exact test p-values indicate how significant a specific chromosome has that number of target genes across all the chromosomes.")
				writeLines(paste(rep("#", 100), collapse = ""))
			}
			
			### make a RDA file with the results
			save(list = c(chr_varNames, "chr_varNames", "README"),
					file = paste0(params[[5]], "All_62_chrom_enrichment.rda"))
			
		}
	}
	
	
	# ******************** which = random_chromosome_distribution  *****************************
	# This is similar to "which = all_hub_target_chromosomes" but it generates
	# distribution significance 1000 times with 1000 random aracne network.
	# This is to make sure that the significance is still meaningful even if we do the same thing
	# with random networks.
	#
	# Needs All_64_Aracne_MI_Fixed.rda, Gene_num_chr.rda, and Random_Network_Gene_Map.rda
	#
	# params[[1]]: a character vector of variable names of aracne networks
	# params[[2]]: "Gene_num_chr.rda" path - used for normalization
	# params[[3]]: "Random_Network_Gene_Map.rda" path - used for calculating random networks
	# params[[4]]: target MI threshold
	# params[[5]]: target p-value threshold
	# params[[6]]: output file directory
	# e.g., params=list("varNames", "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/Gene_num_chr.rda", "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/Random_Network_Gene_Map.rda", 0, 1, "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/chromosome_analysis/")
	# e.g., params=list("varNames", "./Gene_num_chr.rda", "./Random_Network_Gene_Map.rda", 0, 1, "./results/chromosome/")
	
	if(which == "random_chromosome_distribution"){
		if(!is.null(params) && length(params) > 5) {
			### get Aracne variable names
			netNames <- get(as.character(params[[1]]))
			
			### load gene numbers of chromosomes
			load(as.character(params[[2]]))
			total_geneNum_chr <- sum(geneNum)
			
			### prepare the ratio of the number of genes located on every chromosome divided by the total number of genes across all chromosomes
			geneNum_ratio <- geneNum / total_geneNum_chr
			
			### a function to check if a given vector has all zero values
			isVectorAllZero <- function(v) {
				r = TRUE
				for(i in 1:length(v)) {
					if(v[i] != 0) {
						r = FALSE
						break
					}
				}
				return(r)
			}
			
			### get all the union set of hubs
			union_hub_names <- NULL
			for(i in 1:length(netNames)) {
				### get a network
				net <- get(netNames[i])
				
				### union all
				union_hub_names <- union(union_hub_names, rownames(net[[1]]))
			}
			
			### make a chromosome map
			all_genes <- getInteractomeGenes(varNames, count = FALSE)
			chromosome_map <- getChromosome(all_genes)
			
			### load random network map
			load(as.character(params[[3]]))
			
			### iteratively calculate distribution significance
			chr_enrich_varNames <- NULL
			for(p in 1:ncol(randomNetMap)) {
				
				### make a chromosome count data
				chr_count <- list()
				for(i in 1:length(netNames)) {
					### get a network
					net <- get(netNames[i])
					
					### make an empty matrix
					chr_count[[i]] <- matrix(0, length(union_hub_names), 23)
					rownames(chr_count[[i]]) <- union_hub_names
					colnames(chr_count[[i]]) <- paste0(rep("chr", 23), 1:23)
					
					### count the number of target genes for every chromosome
					for(j in 1:length(union_hub_names)) {
						if(length(which(rownames(net[[1]]) == as.character(union_hub_names[j]))) > 0) {
							### get target genes of a given hub
							temp <- net[[2]][[net[[1]][as.character(union_hub_names[j]),2]]]
							temp <- temp[intersect(which(temp[,"MI"] >= as.numeric(params[[4]])),
											which(temp[,"Pvalue"] <= as.numeric(params[[5]]))),,drop=F]
							targets <- rownames(temp)
							
							### map targets to random genes
							targets <- randomNetMap[targets,p]
							
							### +1 for the appropriate chromosome
							for(k in 1:length(targets)) {
								chrIdx <- chromosome_map[targets[k]]
								chr_count[[i]][j,chrIdx] <- chr_count[[i]][j,chrIdx] + 1
							}
						}
					}
					
					### name the list
					if(startsWith(netNames[i], "tcga")) {
						names(chr_count)[i] <- netNames[i]
					} else {
						names(chr_count)[i] <- paste0("gtex_", netNames[i])
					}
				}
				
				### calculate the distribution significance
				dist_sig <- matrix(NA, length(union_hub_names), length(netNames))
				rownames(dist_sig) <- union_hub_names
				colnames(dist_sig) <- varNames
				
				### chi-square test
				for(i in 1:length(netNames)) {
					for(j in 1:length(union_hub_names)) {
						if(!isVectorAllZero(chr_count[[i]][j,])) {
							suppressWarnings(dist_sig[j,i] <- chisq.test(x = chr_count[[i]][j,], p = geneNum_ratio)$p.value)
						}
					}
				}
				
				### set variable for the RDA file
				dist_sig <- apply(dist_sig, 2, as.numeric)
				rownames(dist_sig) <- union_hub_names
				colnames(dist_sig) <- varNames
				assign(paste0("chr_enrich_random", p), dist_sig, envir = globalenv())
				chr_enrich_varNames <- c(chr_enrich_varNames, paste0("chr_enrich_random", p))
				
				### print progress
				writeLines(paste(p, "/", ncol(randomNetMap)))
			}
			
			### make a README function for the RDA file
			README <- function(){
				writeLines(paste(rep("#", 100), collapse = ""))
				writeLines("The \"chr_enrich_varNames\" has all the distribution significance of 1000 random networks")
				writeLines("There are 1000 \"chr_enrich_random*\" objects")
				writeLines("A value in a cell in the matrix indicates the p-value of the Chi-square test of the corresponding hub (row) and the corresponding tissue (column).")
				writeLines("The Chi-square test p-values indicate how significant the distribution of target genes on the chromosomes of a hub is different from the background.")
				writeLines(paste(rep("#", 100), collapse = ""))
			}
			
			### make a RDA file with the results
			save(list = c(chr_enrich_varNames, "chr_enrich_varNames", "README"),
					file = paste0(params[[6]], "all_64_chrom_enrichment_random.rda"))
			
		}
	}
	
	
	# ******************** which = significant_chr_genes *****************************
	# The chromosome analysis was performed with "which = all_hub_target_chromosomes".
	# It produced p-values of each chromosome for each tissue based on
	# the number of same-chromosome-enriched target genes.
	# This function will identify which gene sets are significantly important
	# and what are their biological functions by pathway analysis
	# 
	# The results will be a text file that has significantly important gene set list,
	# and pathway analysis results of those gene sets
	# Needs All_62_chrom_enrichment.rda & All_62_target_chrs.rda
	# params[[1]]: "All_62_chrom_enrichment.rda" file path
	# params[[2]]: "All_62_target_chrs.rda" file path
	# params[[3]]: chromosome significance p-value threshold
	# params[[4]]: pathway analysis FDR threshold
	# params[[5]]: output directory
	# e.g., params=list("//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/All_62_chrom_enrichment.rda", "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/All_62_target_chrs.rda", 1e-50, 0.05, "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/chromosome_analysis/significant_genesets/")
	# e.g., params=list("./results/chromosome/mi_0.2_pv_0.0001/All_62_chrom_enrichment.rda", "./data/RDA_Files/All_62_target_chrs.rda", 1e-50, 0.05, "./results/chromosome/significant_genesets/")
	
	if(which == "significant_chr_genes"){
		if(!is.null(params) && length(params) > 4) {
			### load All_62_chrom_enrichment.rda
			load(as.character(params[[1]]))
			
			### load All_62_target_chrs.rda
			load(as.character(params[[2]]))
			
			### select significant tissue-hub-chromosome cases
			tissue_hub_chrom <- NULL
			for(i in 2:length(chr_varNames)) {
				### get chromosome significance info from each tissue
				chr_sig <- get(chr_varNames[i])
				
				### select indicies of hub-chromosome pairs with significant p-values
				idx <- which(chr_sig < as.numeric(params[[3]]), arr.ind = TRUE)
				if(nrow(idx) > 0) {
					### add tissue-hub-chromosome-pv
					tissue_hub_chrom <- rbind(tissue_hub_chrom,
							cbind(substring(chr_varNames[i], 5),
									rownames(idx),
									idx[,2],
									chr_sig[idx]))
				}
			}
			
			### add symbol of the hubs
			tissue_hub_chrom <- cbind(tissue_hub_chrom[,1:2],
					entrezIDtoSymbol(tissue_hub_chrom[,2]),
					tissue_hub_chrom[,3:4])
			
			### set row and column names
			colnames(tissue_hub_chrom) <- c("Tissue", "Hub", "Hub2", "Chromosome", "P-value")
			rownames(tissue_hub_chrom) <- NULL
			
			### write out the result
			write.table(tissue_hub_chrom,
					file = paste0(params[[5]], "significant_tissue_hub_chrom_", params[[3]], ".txt"),
					sep = "\t", row.names = FALSE, col.names = TRUE)
			
			### make a directory for each tissue
			unique_tissue <- unique(tissue_hub_chrom[,1])
			for(i in 1:length(unique_tissue)) {
				dir.create(paste0(params[[5]], unique_tissue[i]))
			}
			
			### get significant target gene sets
			for(i in 1:nrow(tissue_hub_chrom)) {
				### get target genes
				if(startsWith(tissue_hub_chrom[i,1], "gtex_")) {
					### target genes of the given tissue & the given hub
					temp <- chrInfo[[substring(tissue_hub_chrom[i,1], 6)]][[tissue_hub_chrom[i,2]]]
				} else {
					### target genes of the given tissue & the given hub
					temp <- chrInfo[[tissue_hub_chrom[i,1]]][[tissue_hub_chrom[i,2]]]
				}
				
				### filter the genes with the given chromosome 
				temp <- temp[which(temp == tissue_hub_chrom[i,4])]
				
				### get the target gene names - Entrez IDs
				target_genes <- names(temp)
				
				### pathway analysis
				### KEGG
				pKegg <- pathwayAnalysis_CP(target_genes,
						org = "human",
						database = "KEGG",
						title = paste0(tissue_hub_chrom[i,1], "_", tissue_hub_chrom[i,3], "_",
								tissue_hub_chrom[i,2], "_chr", tissue_hub_chrom[i,4],
								"_target_pathways"),
						pv_threshold = as.numeric(params[[4]]),
						displayNum = 50,
						imgPrint = TRUE,
						dir = paste0(params[[5]], tissue_hub_chrom[i,1], "/"))
				### GO
				pGo <- pathwayAnalysis_CP(target_genes,
						org = "human",
						database = "GO",
						title = paste0(tissue_hub_chrom[i,1], "_", tissue_hub_chrom[i,3], "_",
								tissue_hub_chrom[i,2], "_chr", tissue_hub_chrom[i,4],
								"_target_pathways"),
						pv_threshold = as.numeric(params[[4]]),
						displayNum = 50,
						imgPrint = TRUE,
						dir = paste0(params[[5]], tissue_hub_chrom[i,1], "/"))
				
				### save pathway analysis results in text
				write.table(pKegg, file = paste0(params[[5]], tissue_hub_chrom[i,1], "/",
								tissue_hub_chrom[i,1], "_", tissue_hub_chrom[i,3], "_",
								tissue_hub_chrom[i,2], "_chr", tissue_hub_chrom[i,4],
								"_target_pathways_kegg.txt"),
						sep = "\t", row.names = TRUE, col.names = TRUE)
				write.table(pGo, file = paste0(params[[5]], tissue_hub_chrom[i,1], "/",
								tissue_hub_chrom[i,1], "_", tissue_hub_chrom[i,3], "_",
								tissue_hub_chrom[i,2], "_chr", tissue_hub_chrom[i,4],
								"_target_pathways_go.txt"),
						sep = "\t", row.names = TRUE, col.names = TRUE)
			}
		}
	}
	
	
	# ******************** which = gene_chr_info  *****************************
	# We need to count how many genes lie in each chromsome for every interactome
	# "geneChrInfo" object will have the info by using org.eg.Hs.db package
	#
	# geneChrInfo[["Interactome_Name"]] will have an integer vector of length 23
	#
	# Need "All_62_ARACNE.rda" for getting gene names of all the genes
	#
	# params[[1]]: file path of "All_62_ARACNE.rda"
	# params[[2]]: the output RDA file path
	#
	# e.g., params <- list("//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/All_62_ARACNE.rda", "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/all_64_chr_distribution.rda")
	
	if(which == "gene_chr_info"){
		if(!is.null(params) && length(params) > 1) {
			
			### load gene names
			load(as.character(params[[1]]))
			
			### create an empty list
			geneChrInfo <- list()
			
			### get genes for each tissue
			geneList <- sapply(varNames, function(x) getInteractomeGenes(x, count = FALSE))
			
			### iteratively count the chromosome numbers for each tissue
			for(i in 1:length(varNames)) {
				### get chromosome info of the genes
				chrInfo <- getChromosome(geneList[[i]])
				
				### make an empty vector
				geneChrInfo[[i]] <- rep(0, 23)
				names(geneChrInfo[[i]]) <- paste0(rep("chr", length(geneChrInfo[[i]])), 1:length(geneChrInfo[[i]]))
				
				### count chromosome  numbers
				for(j in 1:length(chrInfo)) {
					geneChrInfo[[i]][chrInfo[j]] <- geneChrInfo[[i]][chrInfo[j]] + 1
				}
				
				### print progress
				writeLines(paste(i, "/", length(varNames)))
			}
			
			### name the list
			names(geneChrInfo) <- varNames
			
			### make a README function for the RDA file
			README <- function(){
				writeLines(paste(rep("#", 100), collapse = ""))
				writeLines("The \"geneChrInfo\" is a list and its length is 62 (based on the number of interactomes)")
				writeLines("Each element of the list has a vector of length 23 indicating the 23 chromosomes")
				writeLines("Each value in the vector represents the number of genes located in the corresponding chromosome")
				writeLines("E.g., geneChrInfo[[\"gtex_AdiposeSub\"]][2] = the number of genes from GTEx_AdiposeSub expression that are located in the Chromosome 2")
				writeLines(paste(rep("#", 100), collapse = ""))
			}
			
			### make a RDA file with the results
			save(list = c("geneChrInfo", "README"), file = as.character(params[[2]]))
			
		} else {
			writeLines("required params do not exist")
		}  
	}
	
	# ******************** which = make_raw_count_rda *****************************
	# Create a RDA file that contains raw count matrices of all the GTEx and
	# TCGA tissues. The Aracne-ready files are normalized, cleaned (removed genes
	# that have 0 or 1 across all samples), and even do not contain complete set
	# of samples (Because 100 <= the number of samples <= 200 is ideal for
	# running Aracne run). Therefore, the files are not suitable for DE analysis
	# between GTEx and TCGA. Here, this function intends to generate raw count
	# matrices for a tissue that already has an Aracne network. There will be
	# no cleaning, no normalization, and no sample selection. The row names
	# will be Entrez IDs.
	#
	# params[[1]]: The directory path of the GTEx raw count files
	#              (a character vector of length 1)
	#              
	# params[[2]]: The RDA file path of the TCGA raw counts
	#              (a character vector of length 1)
	#
	# params[[3]]: The result TCGA raw count RDA file path
	#              (a character vector of length 1)
	#
	# e.g., params = list("//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/ExpressionMatrices/GTEx_processed/separated_counts/not_cleaned/",
	#                      "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/TCGA_33_RAW_COUNTS.rda",
	#                      "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/All_62_raw_counts.rda")
	# e.g., params = list("./results/separated_counts/", "./data/RDA_Files/TCGA_33_RAW_COUNTS.rda", "./data/RDA_Files/All_62_raw_counts.rda")
	
	if(which == "make_raw_count_rda") {
		
		### argument checking
		assertString(params[[1]])
		assertString(params[[2]])
		assertString(params[[3]])
		
		### read GTEx raw counts
		f <- list.files(params[[1]], full.names = TRUE)
		f <- f[which(endsWith(f, ".txt"))]
		gtex_raw_counts <- readRNAseqData(file_names = f)
		
		### shorten the GTEx file names
		### E.g. "Brain-CerebellarHemisphere" -> "BrainCerHem" 
		getShortGTEx <- function(gtex_file_name) {
			fixed_name <- gsub("\\(", "-", gtex_file_name)
			fixed_name <- gsub("\\)", "-", fixed_name)
			fixed_name <- strsplit(fixed_name, "_", fixed = TRUE)[[1]][1]
			
			temp <- strsplit(fixed_name, "-", fixed = TRUE)[[1]]
			
			if(length(temp) > 1) {
				temp2 <- unlist(gregexpr("[A-Z]", temp[2]))
				
				temp3 <- ""
				if((length(temp2) > 1) && (abs(temp2[1] - temp2[2]) > 2)) {
					for(i in 1:2) {
						temp3 <- paste0(temp3, substr(temp[2], temp2[i], temp2[i]+2))
					}  
				} else {
					temp3 <- substr(temp[2], temp2[1], temp2[1]+2)
				}
				
				fixed_name <- paste0(temp[1], temp3)
			} else {
				fixed_name <- temp[1]
			}
			
			return(fixed_name)
		}
		
		### assign tissues that have equal or more than 100 samples to single variables
		gtex_rcntmat_names <- NULL
		for(i in 1:length(gtex_raw_counts)) {
			if(ncol(gtex_raw_counts[[i]]) >= 100) {
				gtex_rcntmat_names <- c(gtex_rcntmat_names, paste0("rcntmat_", getShortGTEx(names(gtex_raw_counts)[i])))
				gtex_raw_counts[[i]] <- gtex_raw_counts[[i]][order(as.numeric(rownames(gtex_raw_counts[[i]]))),]
				assign(gtex_rcntmat_names[length(gtex_rcntmat_names)], gtex_raw_counts[[i]], envir = globalenv())
			}
		}
		
		
		### load TCGA orignal raw counts
		load(params[[2]])
		
		### retain only primary tumor & new primary tumor
		tcga_sample_info <- tcga_sample_info[union(union(which(tcga_sample_info[,"Sample Type"] == "Primary Tumor"),
								which(tcga_sample_info[,"Sample Type"] == "Additional - New Primary")),
						which(tcga_sample_info[,"Sample Type"] == "Primary Blood Derived Cancer - Peripheral Blood")),]
		
		### remove FFPE samples
		tcga_sample_info <- tcga_sample_info[which(tcga_sample_info[,"is_derived_from_ffpe"] == "NO"),]
		
		### order the sample info based on Project ID and Case ID
		tcga_sample_info <- tcga_sample_info[order(tcga_sample_info[,"Project ID"],
						tcga_sample_info[,"Case ID"]),]
		
		### if there are multiple samples per one patient in each tissue, select one with the highest RIN
		unique_tissues <- unique(tcga_sample_info[,"Project ID"])
		rIdx <- NULL
		for(i in 1:length(unique_tissues)) {
			### get indicies for the given tissue
			tempIdx <- which(tcga_sample_info[,"Project ID"] == unique_tissues[i])
			
			### get duplicated indicies in the given tissue
			dupIdx <- tempIdx[which(duplicated(tcga_sample_info[tempIdx, "Case ID"]))]
			
			### if there are duplicates, select one with the highest RIN
			### tie breaker (multiple highest RIN) - select one with the highest lexical order
			if(length(dupIdx) > 0) {
				dups <- unique(tcga_sample_info[dupIdx, "Case ID"])
				
				### collect indicies except one that will remain
				### those indicies will be removed away later
				for(j in 1:length(dups)) {
					dIdx <- which(tcga_sample_info[,"Case ID"] == dups[j])
					rIdx <- c(rIdx, dIdx[order(tcga_sample_info[dIdx, "RIN"])])
					rIdx <- rIdx[-length(rIdx)]
				}
			}
		}
		tcga_sample_info <- tcga_sample_info[-rIdx,]
		
		### make all the raw count matrices have samples only appeared in the tcga_sample_info
		for(i in 1:length(rcnt_matNames)) {
			### get raw count matrix for the given tissue
			rcnt_mat <- get(rcnt_matNames[i])
			
			### retain samples only appeared in the tcga_sample_info (which means they are filtered)
			rcnt_mat <- rcnt_mat[,which(colnames(rcnt_mat) %in% rownames(tcga_sample_info))]
			
			### order the samples in lexical order
			rcnt_mat <- rcnt_mat[,order(colnames(rcnt_mat))]
			
			### change the colnames based on the first 15 characters
			colnames(rcnt_mat) <- substr(colnames(rcnt_mat), 1, 15)
			
			### save the result back to the variable
			assign(rcnt_matNames[i], rcnt_mat, envir = globalenv())
		}
		
		### change the row names based on the first 15 characters
		rownames(tcga_sample_info) <- substr(rownames(tcga_sample_info), 1, 15)
		
		### a function to transfrom Ensembl IDs to Gene symbols
		ensemblIDsToGeneSymbols <- function(ensembl_ids){
			
			### load library
			if(!require(biomaRt, quietly = TRUE)) {
				if(!requireNamespace("BiocManager", quietly = TRUE))
					install.packages("BiocManager")
				BiocManager::install("biomaRt", version = "3.8")
				require(biomaRt, quietly = TRUE)
			}
			
			### mapping information between Ensembl ID and Gene symbol
			mart <- useDataset("hsapiens_gene_ensembl", useEnsembl(biomart="ensembl", version=79))
			mart <- getBM( 
					mart = mart, 
					values = ensembl_ids, 
					filter = c("ensembl_gene_id"), 
					attributes = c("ensembl_gene_id", "external_gene_name"),
					verbose = FALSE
			) 
			
			### Create a dictionary from ensembl id to gene symbol
			ens_to_gene <- as.character(mart$external_gene_name) 
			names(ens_to_gene) <- as.character(mart$ensembl_gene_id) 
			
			### return corresponding gene symbols
			return(ens_to_gene[ensembl_ids])
			
		}
		
		### mapping information between Gene Symbol and Entrez ID (NCBI ID)
		map_symbol_eg <- mappedkeys(org.Hs.egSYMBOL2EG)
		list_symbol2eg <- as.list(org.Hs.egSYMBOL2EG[map_symbol_eg])
		
		### mapping information between Ensembl ID and Entrez ID (NCBI ID)
		map_ensembl_eg <- mappedkeys(org.Hs.egENSEMBL2EG)
		list_ensembl2eg <- as.list(org.Hs.egENSEMBL2EG[map_ensembl_eg])
		
		### change Ensembl IDs to Entrez IDs and retain tissues that have equal or more than 100 samples
		tcga_rcntmat_names <- NULL
		for(i in 1:length(rcnt_matNames)) {
			
			### get raw counts for the given tissue
			df <- get(rcnt_matNames[i])
			
			### change Ensembl IDs to Entrez IDs and retain tissues that have equal or more than 100 samples
			if(ncol(df) >= 100) {
				### store the rcntmat names that are needed
				tcga_rcntmat_names <- c(tcga_rcntmat_names, paste0("rcntmat_", substring(rcnt_matNames[i], 6)))
				
				### get Gencode ID transcript version cleaned
				ensemblIDs <- sapply(rownames(df), function(x) strsplit(x, ".", fixed = TRUE)[[1]][1])
				
				### annotate gene symbols to the raw counts
				df <- data.frame(Ensembl_ID=ensemblIDs,
						Gene_Symbol=ensemblIDsToGeneSymbols(ensemblIDs),
						df,
						stringsAsFactors = FALSE, check.names = FALSE)
				
				### remove rows that do not have gene symbols
				naIdx <- which(is.na(df$Gene_Symbol))
				if(length(naIdx) > 0) {
					df <- df[-naIdx,]
				}
				
				### duplicated gene symbol values
				dups <- unique(df$Gene_Symbol[duplicated(df$Gene_Symbol)])
				
				### initialize indices which should be retained
				retain <- rep(TRUE, nrow(df))
				
				### add up all the expressions for one gene symbol and remove the other rows
				for(dup in dups){
					ind <- which(df$Gene_Symbol == dup)
					df[ind[1],3:ncol(df)] <- apply(df[ind,3:ncol(df)], 2, sum)
					retain[setdiff(ind, ind[1])] <- FALSE
				}
				
				### remove the other rows
				df <- df[retain,]
				
				### get Entrez IDs correspond to the gene symbols
				df <- data.frame(Entrez_ID=as.character(list_symbol2eg[df$Gene_Symbol]),
						df,
						stringsAsFactors = FALSE, check.names = FALSE)
				
				### remove rows with NULL Entrez_ID
				nullIdx <- union(which(df$Entrez_ID == "NULL"), which(is.null(df$Entrez_ID)))
				if(length(nullIdx) > 0) {
					df <- df[-nullIdx,]
				}
				
				### get the accurate Entrez ID using Ensembl ID when there are multiple Entrez IDs mapped to one gene symbol
				entrez_dups <- grep("c", df$Entrez_ID)
				for(dup in entrez_dups) {
					df$Entrez_ID[dup] <- list_ensembl2eg[df$Ensembl_ID[dup]][[1]]
				}
				
				### set rownames with Entrez ID
				rownames(df) <- df$Entrez_ID
				
				### order based on Entrez ID
				df <- df[order(as.numeric(as.character(df$Entrez_ID))),]
				
				### assign processed tcga raw counts to a variable
				assign(tcga_rcntmat_names[length(tcga_rcntmat_names)], as.matrix(df[,-c(1,2,3)]), envir = globalenv())
				
			}
			
		}
		
		### combine the rcntmat names of GTEx and TCGA
		rcntmat_names <- c(gtex_rcntmat_names, tcga_rcntmat_names)
		
		### set README function
		README <- function(){
			writeLines(paste(rep("#", 100), collapse = ""))
			writeLines("A RDA file that contains raw count matrices of all the GTEx and")
			writeLines("TCGA tissues. The Aracne-ready files are normalized, cleaned (removed genes")
			writeLines("that have 0 or 1 across all samples), and even do not contain complete set")
			writeLines("of samples (Because 100 <= the number of samples <= 200 is ideal for")
			writeLines("running Aracne run). Therefore, the files are not suitable for DE analysis")
			writeLines("between GTEx and TCGA. Here, this function intends to generate raw count")
			writeLines("matrices for a tissue that already has an Aracne network. There will be")
			writeLines("no cleaning, no normalization, and no sample selection. The row names")
			writeLines("will be Entrez IDs.")
			writeLines("The \"rcntmat_names\" has all the variable nams of the raw count matrices.")
			writeLines(paste(rep("#", 100), collapse = ""))
		}
		
		### save the processed raw count matrices in a RDA file
		save(list = c("rcntmat_names", rcntmat_names, "README"), file = params[[3]])
		
	}
	
	# ******************************** which = MakeDEGRDA ********************************
	# Perform DE analysis on all the 29 tissue mappings between GTEx and TCGA.
	# The Aracne-ready files are normalized, cleaned (removed genes
	# that have 0 or 1 across all samples), and even do not contain complete set
	# of samples (Because 100 <= the number of samples <= 200 is ideal for
	# running Aracne run). Therefore, the files are not suitable for DE analysis
	# between GTEx and TCGA. This function uses raw count matrices for a tissue
	# that already has an Aracne network. There will be no cleaning, no normalization,
	# and no sample selection. And using the raw counts, DE analysis will be performed
	# between GTEx and TCGA samples. The DE results of 29 mappings will be saved as 
	# a RDA file.
	#
	# params[[1]]: The file path of the "All_62_raw_counts" file
	#              (a character vector of length 1)
	# params[[2]]: The file path of the "GTEx_TCGA_Map.rda" file
	#              (a character vector of length 1)
	# params[[3]]: The result RDA file path
	#              (a character vector of length 1)
	#
	# e.g., params = list("//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/All_62_raw_counts.rda",
	#                     "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/GTEx_TCGA_Map.rda",
	#                     "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/All_29_GTEx_vs_TCGA_DE_Results.rda")
	# e.g., params = list("./data/RDA_Files/All_62_raw_counts.rda", "./data/RDA_Files/GTEx_TCGA_Map.rda", "./data/RDA_Files/All_29_GTEx_vs_TCGA_DE_Results.rda")
	
	if(which == "MakeDEGRDA") {
		
		### argument checking
		assertString(params[[1]])
		assertString(params[[2]])
		assertString(params[[3]])
		
		### load the raw counts
		load(params[[1]])
		load(params[[2]])
		
		### create an empty DER_names
		der_names <- NULL
		
		### for each tissue mapping between GTEx and TCGA, perform DE analysis
		for(i in 1:nrow(GTEx_TCGA_Map)) {
			
			### get gene expressions
			gtex_gexp <- get(paste0("rcntmat_", GTEx_TCGA_Map[i,"GTEx"]))
			tcga_gexp <- get(paste0("rcntmat_tcga_", GTEx_TCGA_Map[i,"TCGA"]))
			
			### common shared genes between GTEx and TCGA
			common_genes <- intersect(rownames(gtex_gexp), rownames(tcga_gexp))
			
			### DE analysis
			gexp <- cbind(gtex_gexp[common_genes,], tcga_gexp[common_genes,])
			group <- c(rep("GTEx", ncol(gtex_gexp)), rep("TCGA", ncol(tcga_gexp)))
			deresult <- deseqWithComparisons(rCnt = gexp, grp = group, exp_class = "GTEx", ctrl_class = "TCGA")
			
			### save the DE result to a variable
			der_names <- c(der_names, paste0("DEG_GTEx_", GTEx_TCGA_Map[i,"GTEx"], "_vs_TCGA_", toupper(GTEx_TCGA_Map[i,"TCGA"])))
			assign(der_names[i], deresult, envir = globalenv())
			
		}
		
		### set README function
		README <- function() {
			writeLines(paste(rep("#", 100), collapse = ""))
			writeLines("A RDA file that contains DE analysis results using raw counts between")
			writeLines("GTEx and TCGA of 29 tissue mappings. The raw counts are not cleaned,")
			writeLines("not normalized, and include all the raw samples. And using the raw counts,")
			writeLines("DE analysis was performed between GTEx and TCGA samples.")
			writeLines("The \"DEG_GTEx_TISSUE_vs_TCGA_TISSUE\" object is a data frame that has")
			writeLines("a DE result of the corresponding GTEx-TCGA comparison of the TISSUE.")
			writeLines(paste(rep("#", 100), collapse = ""))
		}
		
		### save the DE result matrices in a RDA file
		save(list = c("der_names", der_names, "README"), file = params[[3]])
		
	}
	
  # ******************************** which = MakeDEGRDA2 ********************************
  # Perform DE analysis on all the 29 tissue mappings between GTEx and TCGA.
  # The Aracne-ready files are normalized, cleaned (removed genes
  # that have 0 or 1 across all samples), and even do not contain complete set
  # of samples (Because 100 <= the number of samples <= 200 is ideal for
  # running Aracne run). Therefore, the files are not suitable for DE analysis
  # between GTEx and TCGA. This function uses raw count matrices for a tissue
  # that already has an Aracne network. There will be no cleaning, no normalization,
  # and no sample selection. And using the raw counts, DE analysis will be performed
  # between GTEx and TCGA samples. The DE results of 29 mappings will be saved as 
  # a RDA file.
  #
  # params[[1]]: The file path of the "All_62_raw_counts" file
  #              (a character vector of length 1)
  # params[[2]]: The file path of the "GTEx_TCGA_Map.rda" file
  #              (a character vector of length 1)
  # params[[3]]: The result RDA file path
  #              (a character vector of length 1)
  #
  # e.g., params = list("//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/All_62_raw_counts.rda",
  #                     "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/GTEx_TCGA_Map.rda",
  #                     "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/All_29_GTEx_vs_TCGA_DE_Results.rda")
  # e.g., params = list("./data/RDA_Files/All_62_raw_counts.rda", "./data/RDA_Files/GTEx_TCGA_Map.rda", "./data/RDA_Files/All_29_GTEx_vs_TCGA_DE_Results.rda")
  
  if(which == "MakeDEGRDA2") {
    
    ### argument checking
    assertString(params[[1]])
    assertString(params[[2]])
    assertString(params[[3]])
    
    ### load the raw counts
    load(params[[1]])
    load(params[[2]])
    
    ### create an empty DER_names
    der_names <- NULL
    
    ### for each tissue mapping between GTEx and TCGA, perform DE analysis
    for(i in 1:nrow(GTEx_TCGA_Map)) {
      
      ### get gene expressions
      gtex_gexp <- get(paste0("rcntmat_", GTEx_TCGA_Map[i,"GTEx"]))
      tcga_gexp <- get(paste0("rcntmat_tcga_", GTEx_TCGA_Map[i,"TCGA"]))
      
      ### common shared genes between GTEx and TCGA
      common_genes <- intersect(rownames(gtex_gexp), rownames(tcga_gexp))
      
      ### DE analysis
      gexp <- cbind(gtex_gexp[common_genes,], tcga_gexp[common_genes,])
      group <- c(rep("GTEx", ncol(gtex_gexp)), rep("TCGA", ncol(tcga_gexp)))
      deresult <- deseqWithComparisons(rCnt = gexp, grp = group, exp_class = "GTEx", ctrl_class = "TCGA")
      
      ### save the DE result to a variable
      der_names <- c(der_names, paste0("DEG_GTEx_", GTEx_TCGA_Map[i,"GTEx"], "_vs_TCGA_", toupper(GTEx_TCGA_Map[i,"TCGA"])))
      assign(der_names[i], deresult, envir = globalenv())
      
    }
    
    ### set README function
    README <- function() {
      writeLines(paste(rep("#", 100), collapse = ""))
      writeLines("A RDA file that contains DE analysis results using raw counts between")
      writeLines("GTEx and TCGA of 29 tissue mappings. The raw counts are not cleaned,")
      writeLines("not normalized, and include all the raw samples. And using the raw counts,")
      writeLines("DE analysis was performed between GTEx and TCGA samples.")
      writeLines("The \"DEG_GTEx_TISSUE_vs_TCGA_TISSUE\" object is a data frame that has")
      writeLines("a DE result of the corresponding GTEx-TCGA comparison of the TISSUE.")
      writeLines(paste(rep("#", 100), collapse = ""))
    }
    
    ### save the DE result matrices in a RDA file
    save(list = c("der_names", der_names, "README"), file = params[[3]])
    
  }
  
	# ********************* which = generate_exclusivity_scores *************************
	# Computes and "exclusivity score" for every hub gene, a number that measures how strongly 
	# regulons are conserved within interactomes from the same sample collection (GTEx or TCGA) 
	# vs. across collections. A high exclusivity score means that a regulon is conserved more 
	# strongly between pairs of GTEx or TCGA interactomes rather than between pairs involving 
	# one interactome from each collection. Exclusivity scores are computed in several ways:
	# - by summing up log10(FET p-values) from the matrices in tfPairEnrich[[2]]. Specifically, 
	# 		for every hub gene we add up the hub's regulon conservation -log10(p-values) for all 
	#		interactome pairs involving either two GTEx or two TCGA interactomes and then subtract 
	#		-log10(p-values) for interactome pairs involving one interactome from each collection.
	# -	by using a significance threshold T. In this option, we identify all interactome pairs
	#		where the regulon conservation log10(FET p-values) is at most log10(T) and count
	#		the number A of such pairs involving either either two GTEx or two TCGA interactomes
	#		and the number B involving one GTEx and one TCGA interactome. The enrichment score 
	#		is then defined as (A-B).
	# In all score computations, scores are scaled to correct for differences in the numbers of 
	# pairs in the two categories and -Inf values are replaced by the smallest non-infinite value. 
	# Accordingly, 	# high positive are indicative of regulons strongly conserved within but not 
	# across collections.
	#
	# ARGUMENTS
	# * params[[1]]:	character strings, specifying which enrichment score computations option to 
	#		use. Available options are "FET_SUM" and "FET_COUNT"
	# * params[[2]]:	a positive number no more than 1 (e.g., 1e-10), specifying the significance
	#		threshold to use when params[[1]] = "FET_COUNT"
	#
	# The exclusivity scores are stored in the global variable "reg_exclusivity_scores".This is a
	# named vector with one entry per hub gene, containing the computed scores, ranked in decreasing 
	# score order. Vector names are entrez ids for the hub genes. The method also creates a couple
	# more global variables, used by the call makeGraphs(which="hub_FET_enrichment_boxplot").
	if(which == "generate_exclusivity_scores"){
		all_hubs = as.character(getInteractomeGenes(nets=varNames, count=FALSE, hubs_only = TRUE))
		
		# cons_mat (conservation matrix): a matrix with one row per hub and one column
		# per pair of interactomes, containing the log(FET p-value)s from the object
		# tfPairEnrich
		cons_mat = matrix(0, length(all_hubs), choose(length(varNames), 2))
		rownames(cons_mat) = all_hubs
		
		# Convenience object, map interactome pairs to one of three categories: "GTEX" (meaning
		# that both interactomes are from the GTEx collection), "TCGA" (meaning that both
		# interactomes are from the TCGA collection), or "BOTH" (meaning that the two interactomes
		# each come from the GTEx and TCGA collections).
		i_pair_map = rep("", choose(length(varNames), 2))
		k = 1
		for (i in 1:(length(varNames)-1))
			for (j in (i+1):length(varNames)){
				names(i_pair_map)[k] = paste(varNames[i], varNames[j], sep="%")
				i_pair_map[k] = switch(sum(grepl("tcga", c(varNames[i], varNames[j])))+1, "GTEX", "BOTH", "TCGA")
				k = k+1
			}
		colnames(cons_mat) = names(i_pair_map)
		
		for (i in 1:(length(varNames)-1))
			for (j in (i+1):length(varNames)){
				cname = paste(varNames[i], varNames[j], sep="%")
				fet_mat = tfPairEnrich[[2]][[tfPairEnrich[[1]][i ,j]]]
				cons_mat[rownames(fet_mat), cname] = fet_mat[, "log10_pval"]
			}
		# Replace -Inf with the smallest recorded non-infinite log(p-value).
		minp = min(cons_mat[cons_mat != -Inf])
		
		# Compute exclusivity scores. Begin by computing scaling factors to be used for
		# normalizing scores because the number of interactome pairs involving only GTEx
		# or only TCGA interactomes (these are the pairs that contribute positively to the
		# overall score) is different than the number of pairs ivolving one interactome
		# from each collection (these are the pairs that contribute negatively to the 
		# overall score).
		homog_num = choose(sum(grepl("tcga", varNames)),2)+choose(length(varNames)-sum(grepl("tcga", varNames)), 2)
		heter_num = choose(length(varNames), 2) - homog_num
		homog_scale = max(homog_num, heter_num) / homog_num
		heter_scale = max(homog_num, heter_num) / heter_num
		scores = rep(0, length(all_hubs))
		names(scores) = all_hubs
		scores2 = matrix(0, length(all_hubs), 3)
		rownames(scores2) <- all_hubs
		colnames(scores2) <- c("GTEX", "TCGA", "BOTH")
		if (params[[1]] == "FET_COUNT")
			thresh = params[[2]]
		for (hub in all_hubs){
			fets = cons_mat[hub, ]
			fets = replace(fets, fets==-Inf, minp)
			if (params[[1]] == "FET_SUM") {
				scores[hub] = round(-sum(fets[i_pair_map %in% c("TCGA", "GTEX")])*homog_scale + 
								sum(fets[i_pair_map %in% c("BOTH")])*heter_scale)
				scores2[hub, "GTEX"] = -sum(fets[i_pair_map %in% c("GTEX")])*homog_scale
				scores2[hub, "TCGA"] = -sum(fets[i_pair_map %in% c("TCGA")])*homog_scale
				scores2[hub, "BOTH"] = -sum(fets[i_pair_map %in% c("BOTH")])*heter_scale
			}
			else if (params[[1]] == "FET_COUNT") {
				scores[hub] = round(sum(fets[i_pair_map %in% c("TCGA", "GTEX")] <= log10(thresh))*homog_scale - 
								sum(fets[i_pair_map %in% c("BOTH")] <= log10(thresh))*heter_scale)
				scores2[hub, "GTEX"] = sum(fets[i_pair_map %in% c("GTEX")] <= log10(thresh))*homog_scale
				scores2[hub, "TCGA"] = sum(fets[i_pair_map %in% c("TCGA")] <= log10(thresh))*homog_scale
				scores2[hub, "BOTH"] = sum(fets[i_pair_map %in% c("BOTH")] <= log10(thresh))*heter_scale
			} else
				stop("params[[1]] needs to assume a value from c('FET_SUM', 'FET_COUNT')")
    }
		scores = sort(scores, decreasing = TRUE)	
		
		# Assign to global variables, to make accessible to other methods
		assign("reg_exclusivity_scores", scores, envir = globalenv())
		assign("reg_exclusivity_scores_collection", scores2, envir = globalenv())
		assign("reg_fet_cons_mat", cons_mat, envir = globalenv())
		assign("interactome_pair_map", i_pair_map, envir = globalenv())
	}
	
  # ******************************** which = Viper_activity_all2 ********************************
  # Create Viper activity score tables. This is similar to which = viper_acitivity_all, but
  # this is only for TCGA tissues and only for TCGA-GTEx mapping tissues.
  # Match each TCGA cancer A to its corresponding GTEx tissue B. Generate the viper profile
  # for each sample in A by using as reference the "best" 100 samples from B, based on RIN value.
  # Before this analysis we need to re-normalize together all A samples and the 100 samples
  # in the B reference set. Also experiment with the using Combat to code the collection
  # of origin (TCGA or GTEx). I.e., run viper with and without combat for a TCGA tumor and
  # compare the correlation of the viper signatures in each sample, to assess the impact
  # of using Combat.
  #
  # params[[1]]: The file path of the "All_62_raw_counts" file
  #              (a character vector of length 1)
  # params[[2]]: The file path of the "GTEx_TCGA_Map.rda" file
  #              (a character vector of length 1)
  # params[[3]]: The file path of the "TCGA_26_Regulons.rda" file
  #              (a character vector of length 1)
  # params[[4]]: The file path of the "GTEx_Data_V6_SampleData.csv" file
  #              (a character vector of length 1)
  # params[[5]]: Number of samples for reference (null) model
  #              (an integer value)
  # params[[6]]: Method when making viper signature (should be either "zscore", "ttest", or "mean")
  #              (a character vector of length 1)
  # params[[7]]: Permutation number for viper signature
  #              (an integer value)
  # params[[8]]: Batch effect correction between GTEx and TCGA (should be either TRUE or FALSE)
  #              (a boolean value)
  # params[[9]]: Output directory
  #              (A character vector of length 1)
  #
  # e.g., params = list("//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/All_62_raw_counts.rda",
  #                     "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/GTEx_TCGA_Map.rda",
  #                     "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/TCGA_26_Regulons.rda",
  #                     "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/Annotations/GTEx_Data_V6_SampleData.csv",
  #                     100, "zscore", 1000, FALSE,
  #                     "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/results/viper/TCGA/gtex_associated/")
  # e.g., params = list("./data/RDA_Files/All_62_raw_counts.rda", "./data/RDA_Files/GTEx_TCGA_Map.rda", "./data/RDA_Files/TCGA_26_Regulons.rda", "./data/GTEx_Data_V6_SampleData.csv", 100, "zscore", 1000, FALSE, "./results/viper/TCGA/gtex_associated/")
  
  if(which == "Viper_activity_all2") {
    
    ### argument checking
    assertString(params[[1]])
    assertString(params[[2]])
    assertString(params[[3]])
    assertString(params[[4]])
    assertIntegerish(params[[5]])
    assertString(params[[6]])
    assertIntegerish(params[[7]])
    assertLogical(params[[8]])
    assertString(params[[9]])
    
    ### load required data
    load(params[[1]])
    load(params[[2]])
    load(params[[3]])
    
    ### remove unnecessary objects
    rm(list = setdiff(rcntmat_names, union(paste0("rcntmat_", unique(GTEx_TCGA_Map[,1])),
                                           paste0("rcntmat_tcga_", unique(GTEx_TCGA_Map[,2])))))
    rcntmat_names <- union(paste0("rcntmat_", unique(GTEx_TCGA_Map[,1])),
                           paste0("rcntmat_tcga_", unique(GTEx_TCGA_Map[,2])))
    rm(list = setdiff(tcgaRegNames, paste0("regulon_tcga_", unique(GTEx_TCGA_Map[,2]))))
    tcgaRegNames <- paste0("regulon_tcga_", unique(GTEx_TCGA_Map[,2]))
    gc()
    
    ### load sample info
    sampleInfo <- read.csv(params[[4]], check.names = FALSE, stringsAsFactors = FALSE)
    rownames(sampleInfo) <- sampleInfo$SAMPID
    
    ### make RIN info
    rin <- sampleInfo$SMRIN
    names(rin) <- sampleInfo$SAMPID
    
    ### make an empty list for the result RDA
    vipermat <- vector("list", length = nrow(GTEx_TCGA_Map))
    
    ### for each tissue mapping between GTEx and TCGA, perform the analysis
    for(i in 1:nrow(GTEx_TCGA_Map)) {
      
      ### file name for the output
      fileName <- paste("TCGA", toupper(GTEx_TCGA_Map[i,2]), "GTEX", toupper(GTEx_TCGA_Map[i,1]), "ViperMat", params[[5]], params[[6]], params[[7]], sep = "_")
      
      ### get raw counts
      gtex_rawcnt <- get(paste0("rcntmat_", GTEx_TCGA_Map[i,1]))
      tcga_rawcnt <- get(paste0("rcntmat_tcga_", GTEx_TCGA_Map[i,2]))
      
      ### only retain the top params[[5]] samples in GTEx raw counts based on RIN
      target_samples <- names(rin[colnames(gtex_rawcnt)][order(-rin[colnames(gtex_rawcnt)])])
      if(length(target_samples) > as.integer(params[[5]])) {
        target_samples <- target_samples[1:as.integer(params[[5]])]
      }
      gtex_rawcnt <- gtex_rawcnt[,target_samples, drop = FALSE]
      
      ### combine the gtex and tcga raw counts
      combined_rawcnt <- merge(gtex_rawcnt, tcga_rawcnt, by = "row.names")
      rownames(combined_rawcnt) <- combined_rawcnt[,1]
      combined_rawcnt <- combined_rawcnt[,-1]
      combined_rawcnt <- combined_rawcnt[order(as.numeric(rownames(combined_rawcnt))),]
      
      ### normalization
      norm_cnt <- normalizeRNASEQwithVST(combined_rawcnt, filter_thresh = 0)
      
      ### batch effect correction
      if(params[[8]] == TRUE) {
        ### load library
        if(!require(sva, quietly = TRUE)) {
          if(!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager")
          BiocManager::install("sva")
          require(sva, quietly = TRUE)
        }
        
        norm_cnt <- ComBat(dat = as.matrix(norm_cnt),
                           batch = substr(colnames(norm_cnt), 1, 4))
        
        fileName <- paste0(fileName, "_Batch_Corrected")
      }
      
      ### viper signature - an input for viper activity
      vpsig <- viperSignature(eset = as.matrix(norm_cnt[,which(startsWith(colnames(norm_cnt), "TCGA"))]),
                              ref = as.matrix(norm_cnt[,which(startsWith(colnames(norm_cnt), "GTEX"))]),
                              method = as.character(params[[6]]),
                              per = as.integer(params[[7]]),
                              seed = 1,
                              verbose = FALSE)
      
      ### get TCGA regulon object for viper activity
      tcga_regulon <- get(paste0("regulon_tcga_", GTEx_TCGA_Map[i,2]))
      
      # ### remove non-existing genes from the regulon
      # existing_genes <- intersect(rownames(gtex_rawcnt), rownames(tcga_rawcnt))
      # tcga_regulon <- tcga_regulon[which(names(tcga_regulon) %in% existing_genes)]
      # tcga_regulon <- lapply(tcga_regulon, function(x) {
      #   retain_idx <- which(names(x[[1]]) %in% existing_genes)
      #   y <- list(x[[1]][retain_idx], x[[2]][retain_idx])
      #   names(y) <- c("tfmode", "likelihood")
      #   return(y)
      # })
      
      ### compute viper activity table
      vpres <- viper(vpsig, tcga_regulon, verbose = FALSE)
      
      ### write out the viper activity table
      write.table(vpres, file = paste0(params[[9]], fileName, ".txt"), sep = "\t", row.names = TRUE, quote = FALSE)
      
      ### save the result
      vipermat[[i]] <- vpres
      
      ### garbage collection
      gc()
      
    }
    
    ### name the result
    names(vipermat) <- apply(GTEx_TCGA_Map, 1, function(x) paste("TCGA", toupper(x[2]), "GTEX", toupper(x[1]), sep = "_"))
    
    ### make a README
    README <- function() {
      writeLines(paste(rep("#", 100), collapse = ""))
      writeLines("The \"vipermat\" is a list that contains 29 matrices.")
      writeLines("Each matrix is a viper signature (viper activity table) of a TCGA tissue that")
      writeLines("used the reference model as the corresponding GTEx tissue.")
      writeLines("The \"GTEx_TCGA_Map\" has the tissue mapping information between GTEx and TCGA.")
      writeLines("The matrices were generated as the following:")
      writeLines("Match each TCGA cancer A to its corresponding GTEx tissue B. Generate the viper profile")
      writeLines("for each sample in A by using as reference the best 100 samples from B, based on RIN value.")
      writeLines("Before this analysis we need to re-normalize together all A samples and the 100 samples")
      writeLines("in the B reference set. Any existing batch effect between GTEx and TCGA can be")
      writeLines("also considered if the params[[8]] is set as TRUE.")
      writeLines("The \"names(vipermat)\" indicates which viper matrices are for which tissues.")
      writeLines("The format is \"TCGA_[A tissue]_GTEX_[B tissue].")
      writeLines(paste(rep("#", 100), collapse = ""))
    }
    
    ### make RDA file with all the objects
    if(params[[8]] == TRUE) {
      save(list = c("vipermat", "GTEx_TCGA_Map", "README"), file = paste0(params[[9]], paste("All", nrow(GTEx_TCGA_Map), "GTEx", "vs", "TCGA", "ViperMats", "Batch", "Corrected.rda", sep = "_")))
    } else {
      save(list = c("vipermat", "GTEx_TCGA_Map", "README"), file = paste0(params[[9]], paste("All", nrow(GTEx_TCGA_Map), "GTEx", "vs", "TCGA", "ViperMats.rda", sep = "_")))
    }
    
  }
  
  # ******************************** which = Viper_activity_all3 ********************************
  # Create Viper activity score tables. This is similar to which = viper_acitivity_all, but
  # this is only for TCGA tissues and only for TCGA-GTEx mapping tissues.
  # Match each TCGA cancer A to its corresponding GTEx tissue B. Generate the viper profile
  # for each sample in A by using as reference the "best" 100 samples from B, based on RIN value.
  # For the null model, we employ gene expression data of GTEx and TCGA from Schultz's group.
  # They preprocessed GTEx and TCGA data with the same pipeline and batch corrected.
  #
  # params[[1]]: The file path of the "All_21_Schultz_Gene_Expressions.rda" file
  #              (a character vector of length 1)
  # params[[2]]: The file path of the "TCGA_26_Regulons.rda" file
  #              (a character vector of length 1)
  # params[[3]]: The file path of the "GTEx_Data_V6_SampleData.csv" file
  #              (a character vector of length 1)
  # params[[4]]: Number of samples for reference (null) model
  #              (an integer value)
  # params[[5]]: Method when making viper signature (should be either "zscore", "ttest", or "mean")
  #              (a character vector of length 1)
  # params[[6]]: Permutation number for viper signature
  #              (an integer value)
  # params[[7]]: Output directory
  #              (A character vector of length 1)
  #
  # e.g., params = list("//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/All_21_Schultz_Gene_Expressions.rda",
  #                     "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/TCGA_26_Regulons.rda",
  #                     "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/Annotations/GTEx_Data_V6_SampleData.csv",
  #                     100, "zscore", 1000,
  #                     "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/results/viper/TCGA/gtex_associated/")
  # e.g., params = list("./data/RDA_Files/All_21_Schultz_Gene_Expressions.rda", "./data/RDA_Files/TCGA_26_Regulons.rda", "./data/GTEx_Data_V6_SampleData.csv", 100, "zscore", 1000, "./results/viper/TCGA/gtex_associated/")
  
  if(which == "Viper_activity_all3") {
    
    ### argument checking
    assertString(params[[1]])
    assertString(params[[2]])
    assertString(params[[3]])
    assertIntegerish(params[[4]])
    assertString(params[[5]])
    assertIntegerish(params[[6]])
    assertString(params[[7]])
    
    ### load required data
    load(params[[1]])
    load(params[[2]])
    
    ### remove unnecessary objects
    GTEx_TCGA_Map <- GTEx_TCGA_Map[intersect(which(!is.na(GTEx_TCGA_Map[,3])), which(!is.na(GTEx_TCGA_Map[,4]))),]
    rm(list = setdiff(tcgaRegNames, paste0("regulon_tcga_", unique(GTEx_TCGA_Map[,2]))))
    tcgaRegNames <- paste0("regulon_tcga_", unique(GTEx_TCGA_Map[,2]))
    gc()
    
    ### load sample info
    sampleInfo <- read.csv(params[[3]], check.names = FALSE, stringsAsFactors = FALSE)
    rownames(sampleInfo) <- sampleInfo$SAMPID
    
    ### make RIN info
    rin <- sampleInfo$SMRIN
    names(rin) <- sampleInfo$SAMPID
    
    ### make an empty list for the result RDA
    vipermat <- vector("list", length = nrow(GTEx_TCGA_Map))
    
    ### for each tissue mapping between GTEx and TCGA, perform the analysis
    for(i in 1:nrow(GTEx_TCGA_Map)) {
      
      ### file name for the output
      fileName <- paste("TCGA", toupper(GTEx_TCGA_Map[i,2]), "GTEX", toupper(GTEx_TCGA_Map[i,1]), "ViperMat", params[[4]], params[[5]], params[[6]], sep = "_")
      
      ### get raw counts
      gtex_gexp <- get(paste0("Schultz_Gene_Expression_GTEx_", GTEx_TCGA_Map[i,1]))
      tcga_gexp <- get(paste0("Schultz_Gene_Expression_TCGA_", GTEx_TCGA_Map[i,2]))
      
      ### only retain the top params[[4]] samples in GTEx raw counts based on RIN
      target_samples <- names(rin[colnames(gtex_gexp)][order(-rin[colnames(gtex_gexp)])])
      if(length(target_samples) > as.integer(params[[4]])) {
        target_samples <- target_samples[1:as.integer(params[[4]])]
      }
      gtex_gexp <- gtex_gexp[,target_samples, drop = FALSE]
      
      ### combine the gtex and tcga raw counts
      combined_gexp <- merge(gtex_gexp, tcga_gexp, by = "row.names")
      rownames(combined_gexp) <- combined_gexp[,1]
      combined_gexp <- combined_gexp[,-1]
      combined_gexp <- combined_gexp[order(as.numeric(rownames(combined_gexp))),]
      
      ### viper signature - an input for viper activity
      vpsig <- viperSignature(eset = as.matrix(combined_gexp[,which(startsWith(colnames(combined_gexp), "TCGA"))]),
                              ref = as.matrix(combined_gexp[,which(startsWith(colnames(combined_gexp), "GTEX"))]),
                              method = as.character(params[[5]]),
                              per = as.integer(params[[6]]),
                              seed = 1,
                              verbose = FALSE)
      
      ### get TCGA regulon object for viper activity
      tcga_regulon <- get(paste0("regulon_tcga_", GTEx_TCGA_Map[i,2]))
      
      # ### remove non-existing genes from the regulon
      # existing_genes <- intersect(rownames(gtex_gexp), rownames(tcga_gexp))
      # tcga_regulon <- tcga_regulon[which(names(tcga_regulon) %in% existing_genes)]
      # tcga_regulon <- lapply(tcga_regulon, function(x) {
      #   retain_idx <- which(names(x[[1]]) %in% existing_genes)
      #   y <- list(x[[1]][retain_idx], x[[2]][retain_idx])
      #   names(y) <- c("tfmode", "likelihood")
      #   return(y)
      # })
      
      ### compute viper activity table
      vpres <- viper(vpsig, tcga_regulon, verbose = FALSE)
      
      ### write out the viper activity table
      write.table(vpres, file = paste0(params[[7]], fileName, ".txt"), sep = "\t", row.names = TRUE, quote = FALSE)
      
      ### save the result
      vipermat[[i]] <- vpres
      
      ### garbage collection
      gc()
      
    }
    
    ### name the result
    names(vipermat) <- apply(GTEx_TCGA_Map, 1, function(x) paste("TCGA", toupper(x[2]), "GTEX", toupper(x[1]), sep = "_"))
    
    ### make a README
    README <- function() {
      writeLines(paste(rep("#", 100), collapse = ""))
      writeLines("The \"vipermat\" is a list that contains 12 matrices.")
      writeLines("For viper profiles, gene expression data of Schultz's group have been used.")
      writeLines("The gene expression data of both GTEx and TCGA were from the same pipeline and")
      writeLines("batch effect between the two was also corrected.")
      writeLines("Each matrix is a viper signature (viper activity table) of a TCGA tissue that")
      writeLines("used the reference model as the corresponding GTEx tissue.")
      writeLines("The \"GTEx_TCGA_Map\" has the tissue mapping information between GTEx and TCGA.")
      writeLines("The matrices were generated as the following:")
      writeLines("Match each TCGA cancer A to its corresponding GTEx tissue B. Generate the viper profile")
      writeLines("for each sample in A by using as reference the best 100 samples from B, based on RIN value.")
      writeLines("The \"names(vipermat)\" indicates which viper matrices are for which tissues.")
      writeLines("The format is \"TCGA_[A tissue]_GTEX_[B tissue].")
      writeLines(paste(rep("#", 100), collapse = ""))
    }
    
    ### make RDA file with all the objects
    save(list = c("vipermat", "GTEx_TCGA_Map", "README"), file = paste0(params[[7]], paste("All", nrow(GTEx_TCGA_Map), "GTEx", "vs", "TCGA", "ViperMats.rda", sep = "_")))
    
  }
  
  # ************************** which = exclusive_conservation_analysis **************************
  # After running "regulonConservationMode()", we have a matrix with one row per hub gene
  # and three columns titled "GTEX", "TCGA", and "BOTH", listing the number of interactome pairs
  # where the hub gene regulon is conserved at a FET-pvalue above a given p-value threshold.
  # Now starting from the matrix, we would like to give answers to five questions below:
  #
  # 1. How the top activated hubs based on Viper activity (GTEx vs TCGA) are enriched
  #    with the significance score?
  # 2. How the driver mutation genes (that can be derived by TCGA MAF file) are enriched
  #    with the significance score?
  # 3. How the top exclusively conserved hubs are enriched with the Viper NES? 
  # 4. How the DE genes (GTEx vs TCGA) are enriched with the score?
  #
  # For those questions, this function will generate tables that have summaries of
  # [enrichment score, correlation, and p-values], and GSEA figures, correlation plots. 
  #
  # params[[1]]: 1. The name of returned object from regulonConservationMode()
  #              2. or the name of returned object from oneOffs("generate_exclusivity_scores")
  #              if you want to run this function with the second option above,
  #              params[[2]] should be "NONE".
  #              (a character vector of length 1)
  # params[[2]]: The sample collection name of interest (GTEX/TCGA/BOTH/NONE)
  #              The "None" indicates that the there will be no filtering and params[[1]] is
  #              a returned result from oneOffs("generate_exclusivity_scores").
  #              i.e., params[[1]] = "reg_exclusivity_scores" 
  #              (a character vector of length 1)
  # params[[3]]: The cut-off number of non-interesting collections
  #              (e.g., if  params[[2]] == "GTEX" and params[[3]] == 10,
  #               then remove hubs that have the "TCGA" count < 10 and the "BOTH" count < 10
  #              (an integer value)
  # params[[4]]: The file path of the "All_12_GTEx_vs_TCGA_ViperMats.rda" file
  #              (a character vector of length 1)
  # params[[5]]: The number of top active hubs that will be used for Q2 & Q3
  #              (an integer value)
  # params[[6]]: The file path of the "TCGA_33_Driver_Mutation_Genes.rda" file
  #              (a character vector of length 1)
  # params[[7]]: The cut-off FDR of top driver mutation genes that will be used
  #              (a number)
  # params[[8]]: The file path of the "All_12_GTEx_TCGA_DE_Results.rda" file
  #              (a character vector of length 1)
  # params[[9]]: The cut-off adjusted p-value of the DE genes that will be used
  #              (a number)
  # params[[10]]: Output directory
  #              (A character vector of length 1)
  #
  # e.g., params = list("t", "TCGA", 10,
  #                     "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/All_12_GTEx_vs_TCGA_ViperMats.rda",
  #                     100,
  #                     "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/TCGA_33_Driver_Mutation_Genes.rda",
  #                     0.1,
  #                     "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/All_12_GTEx_TCGA_DE_Results.rda",
  #                     0.01,
  #                     "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/results/exclusive_conservation/")
  # e.g., params = list("t", "TCGA", 10,
  #                     "./data/RDA_Files/All_12_GTEx_vs_TCGA_ViperMats.rda",
  #                     100,
  #                     "./data/RDA_Files/TCGA_33_Driver_Mutation_Genes.rda",
  #                     0.1,
  #                     "./data/RDA_Files/All_12_GTEx_TCGA_DE_Results.rda",
  #                     0.01,
  #                     "./results/exclusive_conservation/")
  # e.g., params = list("reg_exclusivity_scores", "NONE", NULL,
  #                     "./data/RDA_Files/All_12_GTEx_vs_TCGA_ViperMats.rda",
  #                     100,
  #                     "./data/RDA_Files/TCGA_33_Driver_Mutation_Genes.rda",
  #                     0.1,
  #                     "./data/RDA_Files/All_12_GTEx_TCGA_DE_Results.rda",
  #                     0.01,
  #                     "./results/exclusive_conservation/")
  #
  # load("./data/RDA_Files/All_62_ARACNE.rda")
  # oneOffs("generate_exclusivity_scores", params = list("FET_SUM"))
  # oneOffs("generate_exclusivity_scores", params = list("FET_COUNT", 1e-20))
  # t <- regulonConservationMode(0.01)
  
  if(which == "exclusive_conservation_analysis") {
    
    ### argument checking
    assertString(params[[1]])
    assertString(params[[2]])
    assertIntegerish(params[[3]])
    assertString(params[[4]])
    assertIntegerish(params[[5]])
    assertString(params[[6]])
    assertNumeric(params[[7]])
    assertString(params[[8]])
    assertNumeric(params[[9]])
    assertString(params[[10]])
    
    ### load required libraries
    if(!require("fgsea", quietly = TRUE)) {
      source("https://bioconductor.org/biocLite.R")
      biocLite("fgsea")
      require("fgsea", quietly = TRUE)
    }
    if(!require(ggplot2, quietly = TRUE)) {
      install.packages("ggplot2")
      require(ggplot2, quietly = TRUE)
    }
    if(!require(xlsx, quietly = TRUE)) {
      install.packages("xlsx")
      require(xlsx, quietly = TRUE)
    }
    
    ### load Viper activity scores
    load(params[[4]])
    
    ### load driver mutation gene info
    load(params[[6]])
    
    if(params[[2]] != "None") {
      ### get exclusivity counts
      exclusivity_cnt <- get(params[[1]])
      
      ### filter the counts with the parameters
      exclusivity_cnt <- exclusivity_cnt[order(exclusivity_cnt[,params[[2]]], decreasing = TRUE),]
      filterIdx <- which(colnames(exclusivity_cnt) != params[[2]])
      for(idx in filterIdx) {
        exclusivity_cnt <- exclusivity_cnt[exclusivity_cnt[,idx] < params[[3]],]
      }
      
      #
      ### Question #1 for each TCGA tissue (Top Active Hubs from Viper)
      #
      ### create a directory for the Q1 results
      dir.create(paste0(params[[10]], "viper"))
      
      ### make an empty list
      top_active_hubs <- vector("list", length = length(vipermat))
      names(top_active_hubs) <- names(vipermat)
      
      ### get top active hubs
      for(i in 1:length(vipermat)) {
        means <- apply(vipermat[[i]], 1, mean)
        means <- means[order(-abs(means))]
        top_active_hubs[[i]] <- entrezIDtoSymbol(names(means)[1:params[[5]]])
      }
      
      ### run GSEA
      gsea_result <- fgsea(pathways = top_active_hubs, stats = exclusivity_cnt[,params[[2]]], nperm = 10000)
      
      ### remove rows with negative NES
      ### because they are highly enriched with 0 counts
      ### and we have no interests on them
      retainIdx <- which(gsea_result[,"NES"] >= 0)
      gsea_result <- gsea_result[retainIdx,]
      top_active_hubs <- top_active_hubs[retainIdx]
      
      ### write the result table
      result_table <- data.frame(Viper_Tissue=gsea_result$pathway,
                                 PVal=gsea_result$pval,
                                 Adj_PVal=gsea_result$padj,
                                 ES=gsea_result$ES,
                                 NES=gsea_result$NES,
                                 Hub_Set_Size=gsea_result$size,
                                 stringsAsFactors = FALSE, check.names = FALSE)
      write.table(result_table, file = paste0(params[[10]], "viper/Viper_Hubs_Enrichment_With_Conservation_Counts.txt"),
                  sep = "\t", row.names = FALSE)
      
      ### plot GSEA results
      for(i in 1:length(top_active_hubs)) {
        png(paste0(params[[10]], "viper/GSEA_", names(top_active_hubs)[i], ".png"),
            width = 1200, height = 1000, res = 130)
        print(plotEnrichment(top_active_hubs[[i]],exclusivity_cnt[,params[[2]]]) +
                labs(title = paste0("GSEA_", names(top_active_hubs)[i])))
        dev.off()
      }
      
      ### plot correlation results
      for(i in 1:length(top_active_hubs)) {
        ### get viper scores for the tissue
        means <- apply(vipermat[[i]], 1, mean)
        names(means) <- entrezIDtoSymbol(names(means))
        
        ### common shared genes
        common_genes <- intersect(names(means), rownames(exclusivity_cnt))
        
        ### draw a scatter plot
        x <- exclusivity_cnt[common_genes,params[[2]]]
        y <- means[common_genes]
        png(paste0(params[[10]], "viper/Cor_", names(top_active_hubs)[i], ".png"),
            width = 1200, height = 1000, res = 130)
        plot(x = x, y = y, pch = 19, col = "black",
             main = paste0("Correlation_", names(top_active_hubs)[i], "\n",
                          "P.Cor = ", round(cor(as.numeric(x), as.numeric(y), use = "pairwise.complete.obs"), 5)),
             xlab = "Exclusive Conservation Counts",
             ylab = "Average Viper NES")
        abline(lm(as.numeric(y)~as.numeric(x)), col="blue", lwd=2)
        dev.off()
      }
      
      #
      ### Question #2 for each TCGA tissue (Driver Mutation Genes)
      #
      ### create a directory for the Q1 results
      dir.create(paste0(params[[10]], "DMG"))
      
      ### make an empty list
      dm_genes <- vector("list", length = length(tcga_driver_mutation_genes))
      names(dm_genes) <- names(tcga_driver_mutation_genes)
      
      ### get mutation driver genes with the cut-off FDR
      for(i in 1:length(tcga_driver_mutation_genes)) {
        tIdx <- which(tcga_driver_mutation_genes[[i]][,"fdr"] < params[[7]])
        if(length(tIdx) == 0) {
          dm_genes[[i]] <- NA
        } else {
          dm_genes[[i]] <- as.character(tcga_driver_mutation_genes[[i]]$Hugo_Symbol[tIdx])
        }
      }
      
      ### run GSEA
      gsea_result <- fgsea(pathways = dm_genes, stats = exclusivity_cnt[,params[[2]]], nperm = 10000)
      
      ### remove rows with negative NES
      ### because they are highly enriched with 0 counts
      ### and we have no interests on them
      retainIdx <- which(gsea_result[,"NES"] >= 0)
      gsea_result <- gsea_result[retainIdx,]
      dm_genes <- dm_genes[which(names(dm_genes) %in% gsea_result$pathway)]
      
      ### write the result table
      result_table <- data.frame(DMG_Tissue=gsea_result$pathway,
                                 PVal=gsea_result$pval,
                                 Adj_PVal=gsea_result$padj,
                                 ES=gsea_result$ES,
                                 NES=gsea_result$NES,
                                 Hub_Set_Size=gsea_result$size,
                                 stringsAsFactors = FALSE, check.names = FALSE)
      write.table(result_table, file = paste0(params[[10]], "DMG/DMG_Enrichment_With_Conservation_Counts.txt"),
                  sep = "\t", row.names = FALSE)
      
      ### plot GSEA results
      for(i in 1:length(dm_genes)) {
        png(paste0(params[[10]], "DMG/GSEA_", names(dm_genes)[i], ".png"),
            width = 1200, height = 1000, res = 130)
        print(plotEnrichment(dm_genes[[i]],exclusivity_cnt[,params[[2]]]) +
                labs(title = paste0("GSEA_", names(dm_genes)[i])))
        dev.off()
      }
    }
    
    #
    ### Question #3 for each TCGA tissue (Top exclusively conservative hubs)
    #
    ### create a directory for the Q3 results
    dir.create(paste0(params[[10]], "ECH"))
    
    ### get top exclusively conservative hubs
    if(params[[2]] != "None") {
      top_ECHs <- rownames(exclusivity_cnt)[1:params[[5]]]
    } else {
      exclusivity_cnt <- get(params[[1]])
      top_ECHs <- entrezIDtoSymbol(names(exclusivity_cnt)[1:params[[5]]])
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
    
    ### iteratively perform GSEA per each case per each sample
    for(case in names(vipermat)) {
      ### create a directory for the given case
      dir.create(paste0(params[[10]], "ECH/", case))
      
      ### GSEA per each sample
      for(i in 1:ncol(vipermat[[case]])) {
        ### get stats ready
        stats <- vipermat[[case]][,i]
        stats <- stats[order(stats)]
        names(stats) <- entrezIDtoSymbol(names(stats))
        
        ### run GSEA
        temp_list <- list(top_ECHs)
        names(temp_list) <- colnames(vipermat[[case]])[i]
        temp <- fgsea(pathways = temp_list, stats = stats, nperm = 10000)
        if(i == 1) {
          gsea_result <- temp
        } else {
          gsea_result <- rbind(gsea_result, temp)
        }
        
        ### progress print
        if(i %% 10 == 0) {
          writeLines(paste(case, "PROGRESS", i, "/", ncol(vipermat[[case]])))
        }
      }
      
      ### calculate FDRs
      gsea_result <- gsea_result[order(gsea_result$pval),]
      gsea_result$padj <- correct_bh(gsea_result$pval)
      
      ### draw GSEA plot if the case's FDR < 0.01
      mIdx <- which(gsea_result$padj < 0.01)
      if(length(mIdx) > 0) {
        for(i in 1:length(mIdx)) {
          ### get stats ready
          stats <- vipermat[[case]][,gsea_result$pathway[mIdx[i]]]
          stats <- stats[order(stats)]
          names(stats) <- entrezIDtoSymbol(names(stats))
          
          ### draw the plot
          png(paste0(params[[10]], "ECH/", case, "/GSEA_", case, "_", gsea_result$pathway[mIdx[i]], ".png"),
              width = 1200, height = 1000, res = 130)
          print(plotEnrichment(top_ECHs, stats) +
                  labs(title = paste0("GSEA_", case, "_", gsea_result$pathway[mIdx[i]])))
          dev.off()
        }
      }
      
      ### write out the GSEA table
      result_table <- data.frame(Viper_Sample=gsea_result$pathway,
                                 PVal=gsea_result$pval,
                                 Adj_PVal=gsea_result$padj,
                                 ES=gsea_result$ES,
                                 NES=gsea_result$NES,
                                 Hub_Set_Size=gsea_result$size,
                                 stringsAsFactors = FALSE, check.names = FALSE)
      write.table(result_table, file = paste0(params[[10]], "ECH/", case, "/ECH_Enrichment_With_Viper_Profiles.txt"),
                  sep = "\t", row.names = FALSE)
    }
    
    ### GSEA of random hubs on the Viper profiles
    ### create a directory for the random results
    dir.create(paste0(params[[10]], "random"))
    
    ### random hubs
    set.seed(1234)
    top_ECHs <- entrezIDtoSymbol(names(exclusivity_cnt)[sample(length(exclusivity_cnt), params[[5]])])
    
    ### GSEA per tissue case
    for(case in names(vipermat)) {
      ### GSEA per each sample
      for(i in 1:ncol(vipermat[[case]])) {
        ### get stats ready
        stats <- vipermat[[case]][,i]
        stats <- stats[order(stats)]
        names(stats) <- entrezIDtoSymbol(names(stats))
        
        ### run GSEA
        temp_list <- list(top_ECHs)
        names(temp_list) <- colnames(vipermat[[case]])[i]
        temp <- fgsea(pathways = temp_list, stats = stats, nperm = 10000)
        if(i == 1) {
          gsea_result <- temp
        } else {
          gsea_result <- rbind(gsea_result, temp)
        }
        
        ### progress print
        if(i %% 10 == 0) {
          writeLines(paste(case, "PROGRESS", i, "/", ncol(vipermat[[case]])))
        }
      }
      
      ### calculate FDRs
      gsea_result <- gsea_result[order(gsea_result$pval),]
      gsea_result$padj <- correct_bh(gsea_result$pval)
      
      ### write out the GSEA table
      result_table <- data.frame(Viper_Sample=gsea_result$pathway,
                                 PVal=gsea_result$pval,
                                 Adj_PVal=gsea_result$padj,
                                 ES=gsea_result$ES,
                                 NES=gsea_result$NES,
                                 Hub_Set_Size=gsea_result$size,
                                 stringsAsFactors = FALSE, check.names = FALSE)
      write.table(result_table, file = paste0(params[[10]], "random/", case, "_RH_Enrichment_With_Viper_Profiles.txt"),
                  sep = "\t", row.names = FALSE)
    }
    
  }
  
  # ***************************** which = two_vipers_ma_plot *****************************
  # Draw a MA plot with two Viper matrices to examine difference between them.
  # params[[1]]: The name of the first Viper matrix
  # params[[2]]: The name of the second Viper matrix
  # params[[3]]: The tissue name
  # params[[4]]: The output directory
  # e.g., params=list("vmat1", "vmat2", "TCGA_BRCA_GTEX_BREAST", "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/results/viper/TCGA/ma_plots/")
  # e.g., params=list("vmat1", "vmat2", "TCGA_BRCA_GTEX_BREAST", "./results/viper/TCGA/ma_plots/")
  
  if (which == "two_vipers_ma_plot"){
    
    ### argument checking
    assertString(params[[1]])
    assertString(params[[2]])
    assertString(params[[3]])
    assertString(params[[4]])
    
    ### get two viper matrices
    vmat1 <- get(params[[1]])
    vmat2 <- get(params[[2]])
    
    ### get shared hubs and samples
    shared_hubs <- intersect(rownames(vmat1), rownames(vmat2))
    shared_samples <- intersect(colnames(vmat1), colnames(vmat2))
    
    ### only retain shared rows and columns
    vmat1 <- vmat1[shared_hubs, shared_samples]
    vmat2 <- vmat2[shared_hubs, shared_samples]
    
    ### create a directory for the tissue
    dir.create(paste0(params[[4]], params[[3]]))
    
    ### draw MA plots
    for(i in 1:length(shared_samples)) {
      x <- apply(cbind(vmat1[,i], vmat2[,i]), 1, mean)
      y <- vmat1[,i] - vmat2[,i]
      
      png(paste0(params[[4]], params[[3]], "/", "MA_Plot_", shared_samples[i], ".png"),
          width = 1200, height = 1000, res = 130)
      plot(x = x, y = y,
           xlab = "Average NES", ylab = "NES Difference",
           main = paste0(params[[3]], "_", shared_samples[i]))
      abline(lm(y~x), col="red")
      dev.off()
    }
    
  }
  
  # ***************************** which = ech_group_difference *****************************
  # From which = exclusive_conservation_analysis of oneOffs(), we have GSEA result tables that
  # represent how the selected hubs are enriched with the input signatures (i.e., Viper profiles).
  # The table contains a column with sample names and a column of their GSEA FDRs.
  # Here, we distinguish the samples into two groups: 1. signicantly enriched, 2. otherwise
  # and then see survival rate, age, gender, subgroup, etc. between the two groups.
  #
  # params[[1]]: The file path of the GSEA result table file
  #              (a character vector of length 1)
  # params[[2]]: FDR cut-off for determining the significant group
  #              (a number)
  # e.g., params=list("//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/results/exclusive_conservation/ECH/reg_exclusivity/top_100_hubs/TCGA_BRCA_GTEX_BREAST/ECH_Enrichment_With_Viper_Profiles.txt",
  #                   0.01)
  # e.g., params=list("./results/exclusive_conservation/ECH/reg_exclusivity/top_100_hubs/TCGA_BRCA_GTEX_BREAST/ECH_Enrichment_With_Viper_Profiles.txt",
  #                   0.01)
  
  if (which == "ech_group_difference") {
    
    ### argument checking
    assertString(params[[1]])
    assertNumeric(params[[2]])
    
    ### load libraries
    if(!require(ggbeeswarm, quietly = TRUE)) {
      install.packages("ggbeeswarm")
      require(ggbeeswarm, quietly = TRUE)
    }
    if(!require(ggpubr, quietly = TRUE)) {
      install.packages("ggpubr")
      require(ggpubr, quietly = TRUE)
    }
    if(!require(gridExtra, quietly = TRUE)) {
      install.packages("gridExtra")
      require(gridExtra, quietly = TRUE)
    }
    if(!require(survminer, quietly = TRUE)) {
      install.packages("survminer")
      require(survminer, quietly = TRUE)
    }
    if(!require(survival, quietly = TRUE)) {
      install.packages("survival")
      require(survival, quietly = TRUE)
    }
    if(!require(TCGAbiolinks, quietly = TRUE)) {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install("TCGAbiolinks")
      require(TCGAbiolinks, quietly = TRUE)
    }
    
    ### load GSEA result table
    gsea_result <- read.table(file = params[[1]], header = TRUE, sep = "\t",
                              stringsAsFactors = FALSE, check.names = FALSE)
    
    ### get the TCGA tissue name from params[[1]]
    tissue <- strsplit(params[[1]], split = "/", fixed = TRUE)[[1]]
    tissue <- tissue[length(tissue)-1]
    tissue <- strsplit(tissue, split = "_", fixed = TRUE)[[1]][2]
    tissue <- paste0("TCGA-", tissue)
    
    ### load TCGA clinical info
    tcga_clinical_info <- GDCquery_clinic(project = tissue, type = "clinical")
    rownames(tcga_clinical_info) <- tcga_clinical_info$submitter_id
    
    ### annotate patient ID to the result table
    gsea_result$Patient_ID <- sapply(gsea_result$Viper_Sample, function(x) substr(x, 1, 12))
    
    ### annotate significance to the result table
    gsea_result$Significance <- "NO"
    gsea_result$Significance[which(gsea_result$Adj_PVal < params[[2]])] <- "YES"
    
    ### annotate survival to the result table
    gsea_result$Day_To_Death <- tcga_clinical_info[gsea_result$Patient_ID, "days_to_death"]
    gsea_result$Vital <- tcga_clinical_info[gsea_result$Patient_ID, "vital_status"]
    
    ### annotate primary_diagnosis to the result table
    gsea_result$Primary_Diagnosis <- tcga_clinical_info[gsea_result$Patient_ID, "primary_diagnosis"]
    
    ### annotate tumor_stage to the result table
    gsea_result$Tumor_Stage <- tcga_clinical_info[gsea_result$Patient_ID, "tumor_stage"]
    
    ### annotate age_at_diagnosis to the result table
    gsea_result$Age_At_Diagnosis <- tcga_clinical_info[gsea_result$Patient_ID, "age_at_diagnosis"]
    
    ### annotate prior_malignancy to the result table
    gsea_result$Prior_Malignancy <- tcga_clinical_info[gsea_result$Patient_ID, "prior_malignancy"]
    
    ### annotate gender to the result table
    gsea_result$Gender <- tcga_clinical_info[gsea_result$Patient_ID, "gender"]
    
    ### annotate race to the result table
    gsea_result$Race <- tcga_clinical_info[gsea_result$Patient_ID, "race"]
    
    ### annotate age_at_index to the result table
    gsea_result$Age_At_Index <- tcga_clinical_info[gsea_result$Patient_ID, "age_at_index"]
    
    ### annotate treatment_type to the result table
    gsea_result$Treatment_Type <- tcga_clinical_info[gsea_result$Patient_ID, "treatment_type"]
    
    ### annotate treatment_or_therapy to the result table
    gsea_result$Treatment_Or_Therapy <- tcga_clinical_info[gsea_result$Patient_ID, "treatment_or_therapy"]
    
    ### write out the result table
    write.table(gsea_result, file = paste0(dirname(params[[1]]), "/ECH_Enrichment_With_Info.txt"),
                sep = "\t", row.names = FALSE)
    
    
    ### create a directory for the survival results
    result_dir <- paste0(dirname(params[[1]]), "/survival/")
    dir.create(result_dir)
    
    ### pie chart with vital status
    pie_data <- data.frame(Vital = c("Alive", "Dead", "Alive", "Dead"),
                           Significance = c("YES", "YES", "NO", "NO"),
                           Number = c(length(intersect(which(gsea_result$Vital == "Alive"),
                                                       which(gsea_result$Significance == "YES"))),
                                      length(intersect(which(gsea_result$Vital == "Dead"),
                                                       which(gsea_result$Significance == "YES"))),
                                      length(intersect(which(gsea_result$Vital == "Alive"),
                                                       which(gsea_result$Significance == "NO"))),
                                      length(intersect(which(gsea_result$Vital == "Dead"),
                                                       which(gsea_result$Significance == "NO")))),
                           stringsAsFactors = FALSE, check.names = FALSE)
    p1 <- ggplot(data = pie_data[which(pie_data$Significance == "YES"),],
                 aes(x = "", y = Number, fill = Vital)) +
            theme_minimal() +
            theme(plot.title = element_text(hjust = 0.5, color = "black")) +
            geom_bar(stat = "identity", width = 1) +
            coord_polar(theta="y") +
            geom_text(aes(label = pie_data$Number[which(pie_data$Significance == "YES")]),
                      position = position_stack(vjust = 0.5)) +
            labs(x = NULL, y = NULL, title = "GSEA FDR < 0.01")
    p2 <- ggplot(data = pie_data[which(pie_data$Significance == "NO"),],
                 aes(x = "", y = Number, fill = Vital)) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, color = "black")) +
      geom_bar(stat = "identity", width = 1) +
      coord_polar(theta="y") +
      geom_text(aes(label = pie_data$Number[which(pie_data$Significance == "NO")]),
                position = position_stack(vjust = 0.5)) +
      labs(x = NULL, y = NULL, title = "GSEA FDR >= 0.01")
    
    ### arrange the plots and print out
    g1 <- arrangeGrob(p1, p2, layout_matrix = rbind(c(1, 2), c(1, 2)),
                      top = paste0("Vital Status Difference Between Two Sample Groups"))
    ggsave(file = paste0(result_dir, "pie_chart_vital_status_", tissue, ".png"), g1, width = 20, height = 12)
    
    if(length(unique(gsea_result$Vital)) > 1 && length(unique(gsea_result$Significance)) > 1) {
      ### beeswarm plot with the day_to_death
      ggplot(gsea_result[which(!is.na(gsea_result$Day_To_Death)),], aes_string(x="Significance", y="Day_To_Death")) +
        theme_classic(base_size = 16) +
        geom_boxplot() +
        geom_beeswarm(aes_string(color="Significance"), na.rm = TRUE) +
        stat_compare_means() +
        labs(x = paste0("GSEA FDR < ", params[[2]]), y = "Survival (Days)") +
        theme(legend.position = "None")
      ggsave(filename = paste0(result_dir, "beeswarm_plot_survival(days)_", tissue, ".png"), width = 12, height = 10)
      
      ### survival plot with the day_to_death
      isSurvPlot <- TRUE
      unique_vital <- unique(gsea_result$Vital)
      unique_sig <- unique(gsea_result$Significance)
      for(i in 1:length(unique_vital)) {
        for(j in 1:length(unique_sig)) {
          if(length(intersect(which(gsea_result$Vital == unique_vital[i]),
                              which(gsea_result$Significance == unique_sig[j]))) < 2) {
            isSurvPlot <- FALSE
            break
          }
        }
      }
      if(isSurvPlot) {
        gsea_result$Day_To_Death <- as.numeric(gsea_result$Day_To_Death)
        gsea_result$Vital[gsea_result$Vital == "Alive"] <- 0
        gsea_result$Vital[gsea_result$Vital == "Dead"] <- 1
        gsea_result$Vital <- as.numeric(gsea_result$Vital)
        
        gsea_result$Significance[which(gsea_result$Significance == "YES")] <- paste0("GSEA FDR < ", params[[2]])
        gsea_result$Significance[which(gsea_result$Significance == "NO")] <- paste0("GSEA FDR >= ", params[[2]])
        
        fit <- survfit(as.formula(paste("Surv(Day_To_Death, Vital)", "~", "Significance")), data = gsea_result)
        p3 <- ggsurvplot(
                fit,
                data = gsea_result,
                title = paste0("Survival Differences Between Two Groups In ", tissue),
                legend.labs = levels(as.factor(gsea_result[,"Significance"])),
                risk.table = TRUE,
                tables.col = "strata",
                pval = TRUE,
                conf.int = TRUE,
                conf.int.style = "ribbon",
                xlab = "Time in Days",
                break.time.by = round(max(gsea_result$Day_To_Death, na.rm = TRUE)/5),
                ggtheme = theme_classic(),
                risk.table.y.text.col = TRUE,
                risk.table.height = 0.25,
                risk.table.y.text = FALSE,
                ncensor.plot = FALSE,
                ncensor.plot.height = 0.25
              )
        ggsave(filename = paste0(result_dir, "survival_plot_", tissue, ".png"),
               plot = print(p3), width = 12, height = 10)
      }
    }
    
    
    ### for every newly added column, compare the values
    for(column in colnames(gsea_result)[11:19]) {
      ### create a directory for the additional results
      result_dir <- paste0(dirname(params[[1]]), "/", column, "/")
      dir.create(result_dir) 
      
      ### if it's integer column, compare the values between [FDR < 0.01] vs [FDR >= 0.01]
      ### if it's categorical charater column, compare FDRs of the unique values of the column
      if(class(gsea_result[,column]) == "integer") {
        if(length(unique(gsea_result$Significance)) > 1) {
          ### beeswarm plot
          ggplot(gsea_result[which(!is.na(gsea_result[,column])),], aes_string(x="Significance", y=column)) +
            theme_classic(base_size = 16) +
            geom_boxplot() +
            geom_beeswarm(aes_string(color="Significance"), na.rm = TRUE) +
            stat_compare_means() +
            labs(x = paste0("GSEA FDR < ", params[[2]]), y = column) +
            theme(legend.position = "None")
          ggsave(filename = paste0(result_dir, "beeswarm_plot_", column, "_", tissue, ".png"), width = 12, height = 10)
        }
      } else if(class(gsea_result[,column]) == "character") {
        if(length(unique(gsea_result[,column])) > 1) {
          ### only retain samples that have the same value more than 5 in total
          retain_idx <- NULL
          for(unique_value in unique(gsea_result[,column])) {
            tempIdx <- which(gsea_result[,column] == unique_value)
            if(length(tempIdx) > 5) {
              retain_idx <- c(retain_idx, tempIdx)
            }
          }
          retain_idx <- intersect(retain_idx, which(!is.na(gsea_result[,"Adj_PVal"])))
          
          ### get median for each group
          medians <- aggregate(as.formula(paste0("Adj_PVal ~ ", column)), gsea_result[retain_idx,], median)
          medians$Adj_PVal <- round(medians$Adj_PVal, 5)
          
          ### beeswarm plot
          ggplot(gsea_result[retain_idx,], aes_string(x=column, y="Adj_PVal")) +
            theme_classic(base_size = 16) +
            geom_boxplot(outlier.shape = NA) +
            geom_text(data = medians, aes(label = Adj_PVal, y = Adj_PVal + 0.5)) +
            stat_compare_means() +
            labs(x = column, y = "GSEA FDR") +
            theme(legend.position = "None", axis.text.x = element_text(angle = 90, hjust = 1))
          ggsave(filename = paste0(result_dir, "beeswarm_plot_", column, "_", tissue, ".png"), width = 12, height = 10)
        }
      }
    }
    
  }
  
  # ***************************** which = survival_ech_difference *****************************
  # From which = exclusive_conservation_analysis of oneOffs(), we have GSEA result tables that
  # represent how the selected hubs are enriched with the input signatures (i.e., Viper profiles).
  # The table contains a column with sample names and a column of their GSEA FDRs.
  # Here, we will divide the samples into two groups: good and bad survival groups and see
  # their GSEA FDR differences.
  #
  # params[[1]]: The file path of the GSEA result table file
  #              (a character vector of length 1)
  # params[[2]]: The survival percentage threshold for determining good or bad survival
  #              e.g., 30 means it will be comparing top 30% vs bottom 30% based on survival
  #              (an integer: should be 1-50)
  # e.g., params=list("//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/results/exclusive_conservation/ECH/reg_exclusivity/top_100_hubs/TCGA_BRCA_GTEX_BREAST/ECH_Enrichment_With_Viper_Profiles.txt",
  #                   30)
  # e.g., params=list("./results/exclusive_conservation/ECH/reg_exclusivity/top_100_hubs/TCGA_BRCA_GTEX_BREAST/ECH_Enrichment_With_Viper_Profiles.txt",
  #                   30)
  
  if (which == "survival_ech_difference") {
    
    ### argument checking
    assertString(params[[1]])
    assertIntegerish(params[[2]])
    
    ### load libraries
    if(!require(ggbeeswarm, quietly = TRUE)) {
      install.packages("ggbeeswarm")
      require(ggbeeswarm, quietly = TRUE)
    }
    if(!require(ggpubr, quietly = TRUE)) {
      install.packages("ggpubr")
      require(ggpubr, quietly = TRUE)
    }
    if(!require(gridExtra, quietly = TRUE)) {
      install.packages("gridExtra")
      require(gridExtra, quietly = TRUE)
    }
    if(!require(TCGAbiolinks, quietly = TRUE)) {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install("TCGAbiolinks")
      require(TCGAbiolinks, quietly = TRUE)
    }
    
    ### load GSEA result table
    gsea_result <- read.table(file = params[[1]], header = TRUE, sep = "\t",
                              stringsAsFactors = FALSE, check.names = FALSE)
    
    ### get the TCGA tissue name from params[[1]]
    tissue <- strsplit(params[[1]], split = "/", fixed = TRUE)[[1]]
    tissue <- tissue[length(tissue)-1]
    tissue <- strsplit(tissue, split = "_", fixed = TRUE)[[1]][2]
    tissue <- paste0("TCGA-", tissue)
    
    ### load TCGA clinical info
    tcga_clinical_info <- GDCquery_clinic(project = tissue, type = "clinical")
    rownames(tcga_clinical_info) <- tcga_clinical_info$submitter_id
    
    ### annotate patient ID to the result table
    gsea_result$Patient_ID <- sapply(gsea_result$Viper_Sample, function(x) substr(x, 1, 12))
    
    ### annotate survival to the result table
    gsea_result$Day_To_Death <- tcga_clinical_info[gsea_result$Patient_ID, "days_to_death"]
    gsea_result$Vital <- tcga_clinical_info[gsea_result$Patient_ID, "vital_status"]
    
    ### annotate primary_diagnosis to the result table
    gsea_result$Primary_Diagnosis <- tcga_clinical_info[gsea_result$Patient_ID, "primary_diagnosis"]
    
    ### annotate tumor_stage to the result table
    gsea_result$Tumor_Stage <- tcga_clinical_info[gsea_result$Patient_ID, "tumor_stage"]
    
    ### annotate age_at_diagnosis to the result table
    gsea_result$Age_At_Diagnosis <- tcga_clinical_info[gsea_result$Patient_ID, "age_at_diagnosis"]
    
    ### annotate prior_malignancy to the result table
    gsea_result$Prior_Malignancy <- tcga_clinical_info[gsea_result$Patient_ID, "prior_malignancy"]
    
    ### annotate gender to the result table
    gsea_result$Gender <- tcga_clinical_info[gsea_result$Patient_ID, "gender"]
    
    ### annotate race to the result table
    gsea_result$Race <- tcga_clinical_info[gsea_result$Patient_ID, "race"]
    
    ### annotate age_at_index to the result table
    gsea_result$Age_At_Index <- tcga_clinical_info[gsea_result$Patient_ID, "age_at_index"]
    
    ### annotate treatment_type to the result table
    gsea_result$Treatment_Type <- tcga_clinical_info[gsea_result$Patient_ID, "treatment_type"]
    
    ### annotate treatment_or_therapy to the result table
    gsea_result$Treatment_Or_Therapy <- tcga_clinical_info[gsea_result$Patient_ID, "treatment_or_therapy"]
    
    ### only keep samples that have survival data
    gsea_result <- gsea_result[!is.na(gsea_result$Day_To_Death),]
    
    ### annotate -log10(FDR)
    gsea_result$logFDR <- -log10(gsea_result$Adj_PVal)
    
    ### annotate survival status
    gsea_result$Survival_Status <- "Middle"
    gsea_result <- gsea_result[order(gsea_result$Day_To_Death),]
    thresholdNum <- floor(nrow(gsea_result) * as.numeric(params[[2]]) / 100)
    if(thresholdNum > 0 && thresholdNum < (nrow(gsea_result)/2)) {
      gsea_result$Survival_Status[1:thresholdNum] <- "Bad"
      gsea_result$Survival_Status[(nrow(gsea_result)-thresholdNum+1):nrow(gsea_result)] <- "Good"
    } else {
      gsea_result <- NULL
    }
    
    ### the directory for the survival results
    result_dir <- paste0(dirname(params[[1]]), "/survival/")
    dir.create(result_dir, showWarnings = FALSE)
    
    ### beeswarm plot with the day_to_death
    if(!is.null(gsea_result)) {
      max_y <- max(gsea_result$logFDR)
      ggplot(gsea_result, aes_string(x="Survival_Status", y="logFDR")) +
        theme_classic(base_size = 16) +
        geom_boxplot() +
        geom_beeswarm(aes_string(color="Survival_Status"), na.rm = TRUE) +
        stat_compare_means(label.y = max_y*1.1) +
        ggtitle(paste0("ECHs on Viper(", tissue, ") -log10(GSEA FDR) between Good and Bad Survival")) +
        labs(x = paste0("Survival Status Top-Bottom ", params[[2]], "%"), y = "-log10(GSEA FDR)") +
        theme(legend.position = "None")
      ggsave(filename = paste0(result_dir, "beeswarm_plot_FDR_", tissue, ".png"), width = 12, height = 10)
    }
    
  }
  
  # ***************************** which = ech_cosmic_analysis *****************************
  # Examine significancy of target genes of the exclusively conserved hubs based on
  # enrichment with Cosmic Cancer Gene Census (CGC). In each tissue, for all the regulons,
  # compute how many of those regulons are enriched with the CGC. Using Fisher's exact test,
  # calculate p-value of the top N ECHs (similar to pahtway analysis in each regulon x N).
  # Plus, calculate p-values of all the regulons in each network and plot their p-values.
  # The result will be written as TXT and PNG under the params[[4]] directory.
  #
  # params[[1]]: The file path of the Aracne RDA file (All_62_ARACNE.rda)
  #              (a character vector of length 1)
  # params[[2]]: The file path of Cosmic Cancer Gene Census
  #              (a character vector of length 1)
  # params[[3]]: The number of top ECHs that will be tested
  #              (a number)
  # params[[4]]: The file path of the GTEx-TCGA tissue mapping info (GTEx_TCGA_Map.rda)
  #              (a character vector of length 1)
  # params[[5]]: Directory path for the results
  #              (a character vector of length 1)
  #
  # e.g., params=list("//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/All_62_ARACNE.rda",
  #                   "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/Cosmic/Cosmic_Census_100419_all.tsv",
  #                   100,
  #                   "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/GTEx_TCGA_Map.rda",
  #                   "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/exclusive_conservation/ECH/Cosmic/")
  # e.g., params=list("./data/RDA_Files/All_62_ARACNE.rda",
  #                   "./data/Cosmic/Cosmic_Census_100419_all.tsv",
  #                   100,
  #                   "./data/RDA_Files/GTEx_TCGA_Map.rda",
  #                   "./results/exclusive_conservation/ECH/Cosmic/")
  
  if (which == "ech_cosmic_analysis") {
    
    ### argument checking
    assertString(params[[1]])
    assertString(params[[2]])
    assertNumeric(params[[3]])
    assertString(params[[4]])
    assertString(params[[5]])
    
    ### load data
    load(params[[1]], envir = globalenv())
    cgc <- read.table(file = params[[2]], header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE, check.names = FALSE)
    load(params[[4]])
    Sys.sleep(3)
    
    ### get top ECHs
    oneOffs("generate_exclusivity_scores", params = list("FET_SUM"))
    echs <- names(reg_exclusivity_scores)[1:params[[3]]]
    
    ### get TCGA Aracne names
    tcga_aracne_names <- varNames[grep("tcga", varNames)]
    
    ### perform analysis for each TCGA tissue
    for(aracne_name in tcga_aracne_names) {
      ### get Aracne network
      aracne <- get(aracne_name)
      
      ### get total genes in the network
      interactome_total_genes <- as.character(getInteractomeGenes(aracne_name, count = FALSE))
      
      ### get cancer genes in the network
      interactome_cancer_genes <- intersect(interactome_total_genes, as.character(cgc$`Entrez GeneId`))
      
      ### make an empty data frame
      result_table <- data.frame(matrix(NA, params[[3]], 7))
      colnames(result_table) <- c("Hub_Gene_Symbol", "Hub_Entrez_ID", "Enriched_Regulons_Gene_Symbol",
                                  "Enriched_Regulons_Entrez_ID", "PVal", "Enrichment_Count", "Background")
      
      ### hub names
      result_table$Hub_Entrez_ID <- echs
      result_table$Hub_Gene_Symbol <- entrezIDtoSymbol(echs)
      rownames(result_table) <- echs
      
      ### Fisher's exact test for each hub
      for(hub in echs) {
        ### get target genes
        target_genes <- rownames(aracne[[2]][[hub]])
        
        ### compute enriched genes
        enriched_genes <- intersect(target_genes, interactome_cancer_genes)
        if(length(enriched_genes) > 0) {
          result_table[hub,"Enriched_Regulons_Entrez_ID"] <- paste(enriched_genes, collapse = "/")
          result_table[hub,"Enriched_Regulons_Gene_Symbol"] <- paste(entrezIDtoSymbol(enriched_genes), collapse = "/")
          result_table[hub,"Enrichment_Count"] <- paste0(length(enriched_genes), "/", length(target_genes))
        } else {
          result_table[hub,"Enrichment_Count"] <- paste0("0", "/", length(target_genes))
        }
        result_table[hub,"Background"] <- paste0(length(interactome_cancer_genes), "/", length(interactome_total_genes))
        
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
        result_table[hub,"PVal"] <- fisher.test(matrix(c(x, z, y, w), 2, 2), alternative = "greater")$p.value
      }
      
      ### print out the result table
      write.table(result_table, file = paste0(params[[5]], aracne_name, "_ech_cosmic_enrichment.txt"),
                  sep = "\t", row.names = FALSE)
      
      ### compute Fisher's exact test p-values for all the hubs in the TCGA Aracne network
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
      
      ### plot all the hubs
      barplot_data <- all_enrichment_pvs[order(all_enrichment_pvs)]
      barplot_data <- log10(barplot_data)
      barplot_data[echs] <- -barplot_data[echs]
      colors <- rep("black", length(barplot_data))
      names(colors) <- names(barplot_data)
      colors[echs] <- "red"
      png(paste0(params[[5]], aracne_name, "_top_", params[[3]], "_echs_cosmic_enrichment_pvals.png"),
          width = 2000, height = 1400, res = 130)
      barplot(barplot_data, col = colors, border = colors,
              main = paste(toupper(aracne_name), "Cosmic Cancer Gene Census Enrichment"),
              xaxt = "n", ylim = c(min(barplot_data, na.rm=TRUE)*1.3, max(barplot_data, na.rm=TRUE)*1.3),
              xlab = "All the hubs in the Aracne network",
              ylab = "log10(Enrichment p-values)")
      legend("topright", legend = c(paste0("Top ", params[[3]], " ECHs"), "Others"), col=c("red", "black"), lty = 1)
      dev.off()
      
      ### plot enrichment
      input <- list(echs)
      names(input) <- paste0(aracne_name, "_top_", params[[3]], "_echs_cosmic_enrichment_GSEA")
      temp <- run_gsea(gene_list = input, signature = list(abs(barplot_data)),
                       printPlot = TRUE, fdr_cutoff = 1,
                       printPath = paste0(params[[5]]))
      
      ### plot all the hubs (ordering based on regulon size)
      barplot_data <- barplot_data[all_hubs]
      barplot_data <- -abs(barplot_data)
      barplot_data[echs] <- -barplot_data[echs]
      colors <- rep("black", length(barplot_data))
      names(colors) <- names(barplot_data)
      colors[echs] <- "red"
      png(paste0(params[[5]], aracne_name, "_top_", params[[3]], "_echs_cosmic_enrichment_pvals_RS.png"),
          width = 2000, height = 1400, res = 130)
      barplot(barplot_data, col = colors, border = colors,
              main = paste(toupper(aracne_name), "Cosmic Cancer Gene Census Enrichment"),
              xaxt = "n", ylim = c(min(barplot_data, na.rm=TRUE)*1.3, max(barplot_data, na.rm=TRUE)*1.3),
              xlab = "All the hubs in the Aracne network (Ordered by regulon size)",
              ylab = "log10(Enrichment p-values)")
      legend("topright", legend = c(paste0("Top ", params[[3]], " ECHs"), "Others"), col=c("red", "black"), lty = 1)
      dev.off()
      
      ### get GTEx Aracne name
      gtex_aracne_name <- GTEx_TCGA_Map[which(GTEx_TCGA_Map[,"TCGA"] == strsplit(aracne_name, "_", TRUE)[[1]][2])[1],"GTEx"]
      
      if(!is.na(gtex_aracne_name)) {
        ### get GTEx Aracne network
        gtex_aracne <- get(gtex_aracne_name)
        
        ### get total genes in the GTEx network
        gtex_interactome_total_genes <- as.character(getInteractomeGenes(gtex_aracne_name, count = FALSE))
        
        ### get cancer genes in the  GTEx network
        gtex_interactome_cancer_genes <- intersect(gtex_interactome_total_genes, as.character(cgc$`Entrez GeneId`))
        
        ### compute Fisher's exact test p-values for all the hubs in the GTEx Aracne network
        gtex_all_hubs <- rownames(gtex_aracne[[1]])
        gtex_all_enrichment_pvs <- rep(0, length(gtex_all_hubs))
        names(gtex_all_enrichment_pvs) <- gtex_all_hubs
        for(hub in gtex_all_hubs) {
          ### get target genes
          target_genes <- rownames(gtex_aracne[[2]][[hub]])
          
          ### compute enriched genes
          enriched_genes <- intersect(target_genes, gtex_interactome_cancer_genes)
          
          ### calculate p-value
          ### Fisher's exact test
          ###
          ###                 regulon   no-regulon
          ###               -----------------------
          ### cancer gene   |   X           Y
          ### no-cancer gene|   Z           W
          x <- length(enriched_genes)
          y <- length(gtex_interactome_cancer_genes) - x
          z <- length(target_genes) - x
          w <- length(gtex_interactome_total_genes) - x - y - z
          
          ### Fisher's exact test p-value
          gtex_all_enrichment_pvs[hub] <- fisher.test(matrix(c(x, z, y, w), 2, 2), alternative = "greater")$p.value
        }
        
        ### get common hubs between GTEx and TCGA networks
        common_hubs <- intersect(all_hubs, gtex_all_hubs)
        
        ### plot all the hubs (top: TCGA - ordered, bottom: GTEx)
        tcga_barplot_data <- all_enrichment_pvs[common_hubs]
        tcga_barplot_data <- tcga_barplot_data[order(tcga_barplot_data)]
        tcga_barplot_data <- -log10(tcga_barplot_data)
        
        gtex_barplot_data <- gtex_all_enrichment_pvs[names(tcga_barplot_data)]
        gtex_barplot_data <- log10(gtex_barplot_data)
        
        png(paste0(params[[5]], aracne_name, "_top_", params[[3]], "_echs_cosmic_enrichment_pvals_tcga_ordered.png"),
            width = 2000, height = 1400, res = 100)
        par(mai = c(7, 0.5, 0.5, 0))
        barplot(tcga_barplot_data, axes = FALSE, xaxt = "n")
        barplot(gtex_barplot_data, add = TRUE, axes = FALSE, xaxt = "n")
        abline(h = 0, col = "white")
        mtext("TCGA ORDERED -log10(Enrichment p-values)", side = 2, at = 5)
        mtext("GTEx log10(Enrichment p-values)", side = 2, at = -5)
        mtext(paste(toupper(aracne_name), "Cosmic Cancer Gene Census Enrichment with Corresponding GTEx Hubs"), side = 3)
        dev.off()
        
        ### plot all the hubs (top: TCGA, bottom: GTEx - ordered)
        gtex_barplot_data <- gtex_barplot_data[order(gtex_barplot_data)]
        tcga_barplot_data <- tcga_barplot_data[names(gtex_barplot_data)]
        png(paste0(params[[5]], aracne_name, "_top_", params[[3]], "_echs_cosmic_enrichment_pvals_gtex_ordered.png"),
            width = 2000, height = 1400, res = 100)
        par(mai = c(7, 0.5, 0.5, 0))
        barplot(tcga_barplot_data, axes = FALSE, xaxt = "n")
        barplot(gtex_barplot_data, add = TRUE, axes = FALSE, xaxt = "n")
        abline(h = 0, col = "white")
        mtext("TCGA -log10(Enrichment p-values)", side = 2, at = 5)
        mtext("GTEx ORDERED log10(Enrichment p-values)", side = 2, at = -6)
        mtext(paste(toupper(aracne_name), "Cosmic Cancer Gene Census Enrichment with Corresponding GTEx Hubs"), side = 3)
        dev.off()
        
        ### plot all the hubs (top: TCGA - ordered, bottom: GTEx - ordered)
        tcga_barplot_data <- tcga_barplot_data[order(-tcga_barplot_data)]
        png(paste0(params[[5]], aracne_name, "_top_", params[[3]], "_echs_cosmic_enrichment_pvals_not_matched.png"),
            width = 2000, height = 1400, res = 100)
        par(mai = c(7, 0.5, 0.5, 0))
        barplot(tcga_barplot_data, axes = FALSE, xaxt = "n")
        barplot(gtex_barplot_data, add = TRUE, axes = FALSE, xaxt = "n")
        abline(h = 0, col = "white")
        mtext("TCGA ORDERED -log10(Enrichment p-values)", side = 2, at = 5)
        mtext("GTEx ORDERED log10(Enrichment p-values)", side = 2, at = -6)
        mtext(paste(toupper(aracne_name), "Cosmic Cancer Gene Census Enrichment with GTEx Hubs - NOT MATCHED"), side = 3)
        dev.off()
      }
    }
    
  }
  
  # ***************************** which = viper_cosmic_analysis *****************************
  # Examine significancy of target genes of the high Viper activity hubs based on
  # enrichment with Cosmic Cancer Gene Census (CGC). In each tissue, for all the regulons,
  # compute how many of those regulons are enriched with the CGC. Using Fisher's exact test,
  # calculate p-value of the top N hubs (similar to pahtway analysis in each regulon x N).
  # Plus, calculate p-values of all the regulons in each network and plot their p-values.
  # The result will be written as TXT and PNG under the params[[4]] directory.
  # This is similar to "which = ech_cosmic_analysis" but this time, high viper activity
  # hubs will be used instead of highly exclusively conserved hubs. 
  #
  # params[[1]]: The file path of the Aracne RDA file (All_62_ARACNE.rda)
  #              (a character vector of length 1)
  # params[[2]]: The file path of the Viper profile RDA file (All_12_GTEx_vs_TCGA_ViperMats.rda)
  #              (a character vector of length 1)
  # params[[3]]: The file path of Cosmic Cancer Gene Census
  #              (a character vector of length 1)
  # params[[4]]: The number of top Viper activity hubs that will be tested
  #              (a number)
  # params[[5]]: Directory path for the results
  #              (a character vector of length 1)
  #
  # e.g., params=list("//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/All_62_ARACNE.rda",
  #                   "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/All_12_GTEx_vs_TCGA_ViperMats.rda",
  #                   "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/Cosmic/Cosmic_Census_100419_all.tsv",
  #                   100,
  #                   "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/viper/TCGA/Cosmic/")
  # e.g., params=list("./data/RDA_Files/All_62_ARACNE.rda",
  #                   "./data/RDA_Files/All_12_GTEx_vs_TCGA_ViperMats.rda",
  #                   "./data/Cosmic/Cosmic_Census_100419_all.tsv",
  #                   100,
  #                   "./results/viper/TCGA/Cosmic/")
  
  if (which == "viper_cosmic_analysis") {
    
    ### argument checking
    assertString(params[[1]])
    assertString(params[[2]])
    assertString(params[[3]])
    assertNumeric(params[[4]])
    assertString(params[[5]])
    
    ### load data
    load(params[[1]])
    load(params[[2]])
    cgc <- read.table(file = params[[3]], header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE, check.names = FALSE)
    
    ### for every available tissue that has Viper profile
    for(viper_name in names(vipermat)) {
      ### get viper profile
      viper_profile <- vipermat[[viper_name]]
      
      ### get Aracne network
      aracne_name <- tolower(paste(strsplit(viper_name, split = "_", fixed = TRUE)[[1]][1:2], collapse = "_"))
      aracne <- get(aracne_name)
      
      ### get total genes in the network
      interactome_total_genes <- as.character(getInteractomeGenes(aracne_name, count = FALSE))
      
      ### get cancer genes in the network
      interactome_cancer_genes <- intersect(interactome_total_genes, as.character(cgc$`Entrez GeneId`))
      
      ### highest Viper NES hub matrix
      viper_hub_mat <- matrix(NA, nrow(viper_profile), ncol(viper_profile))
      colnames(viper_hub_mat) <- colnames(viper_profile)
      for(colname in colnames(viper_hub_mat)) {
        viper_hub_mat[,colname] <- names(viper_profile[,colname][order(-abs(viper_profile[,colname]))])
      }
      
      ### give ranks to the hubs based on Viper NES across all samples
      hub_ranks <- rep(0, nrow(viper_profile))
      names(hub_ranks) <- rownames(viper_profile)
      for(hub in names(hub_ranks)) {
        for(colname in colnames(viper_hub_mat)) {
          hub_ranks[hub] <- hub_ranks[hub] + which(viper_hub_mat[,colname] == hub)
        }
      }
      hub_ranks <- hub_ranks[order(hub_ranks)]
      
      ### get top Viper activity hubs
      top_viper_hubs <- names(hub_ranks)[1:params[[4]]]
      
      ### make an empty data frame
      result_table <- data.frame(matrix(NA, params[[4]], 7))
      colnames(result_table) <- c("Hub_Gene_Symbol", "Hub_Entrez_ID", "Enriched_Regulons_Gene_Symbol",
                                  "Enriched_Regulons_Entrez_ID", "PVal", "Enrichment_Count", "Background")
      
      ### hub names
      result_table$Hub_Entrez_ID <- top_viper_hubs
      result_table$Hub_Gene_Symbol <- entrezIDtoSymbol(top_viper_hubs)
      rownames(result_table) <- top_viper_hubs
      
      ### Fisher's exact test for each hub
      for(hub in top_viper_hubs) {
        ### get target genes
        target_genes <- rownames(aracne[[2]][[hub]])
        
        ### compute enriched genes
        enriched_genes <- intersect(target_genes, interactome_cancer_genes)
        if(length(enriched_genes) > 0) {
          result_table[hub,"Enriched_Regulons_Entrez_ID"] <- paste(enriched_genes, collapse = "/")
          result_table[hub,"Enriched_Regulons_Gene_Symbol"] <- paste(entrezIDtoSymbol(enriched_genes), collapse = "/")
          result_table[hub,"Enrichment_Count"] <- paste0(length(enriched_genes), "/", length(target_genes))
        } else {
          result_table[hub,"Enrichment_Count"] <- paste0("0", "/", length(target_genes))
        }
        result_table[hub,"Background"] <- paste0(length(interactome_cancer_genes), "/", length(interactome_total_genes))
        
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
        result_table[hub,"PVal"] <- fisher.test(matrix(c(x, z, y, w), 2, 2), alternative = "greater")$p.value
      }
      
      ### print out the result table
      write.table(result_table, file = paste0(params[[5]], aracne_name, "_tvh_cosmic_enrichment.txt"),
                  sep = "\t", row.names = FALSE)
      
      ### compute Fisher's exact test p-values for all the hubs in the Aracne network
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
      
      ### plot all the hubs
      barplot_data <- all_enrichment_pvs[order(all_enrichment_pvs)]
      barplot_data <- log10(barplot_data)
      barplot_data[top_viper_hubs] <- -barplot_data[top_viper_hubs]
      colors <- rep("black", length(barplot_data))
      names(colors) <- names(barplot_data)
      colors[top_viper_hubs] <- "red"
      png(paste0(params[[5]], aracne_name, "_top_", params[[4]], "_tvhs_cosmic_enrichment_pvals.png"),
          width = 2000, height = 1400, res = 130)
      barplot(barplot_data, col = colors, border = colors,
              main = paste(toupper(aracne_name), "Cosmic Cancer Gene Census Enrichment"),
              xaxt = "n", ylim = c(min(barplot_data, na.rm=TRUE)*1.3, max(barplot_data, na.rm=TRUE)*1.3),
              xlab = "All the hubs in the Aracne network",
              ylab = "log10(Enrichment p-values)")
      legend("topright", legend = c(paste0("Top ", params[[4]], " Viper Hubs"), "Others"), col=c("red", "black"), lty = 1)
      dev.off()
      
      ### plot enrichment
      input <- list(top_viper_hubs)
      names(input) <- paste0(aracne_name, "_top_", params[[4]], "_tvhs_cosmic_enrichment_GSEA")
      temp <- run_gsea(gene_list = input, signature = list(abs(barplot_data)),
                printPlot = TRUE, fdr_cutoff = 1,
                printPath = paste0(params[[5]]))
      
      ### plot all the hubs (ordering based on regulon size)
      barplot_data <- barplot_data[all_hubs]
      barplot_data <- -abs(barplot_data)
      barplot_data[top_viper_hubs] <- -barplot_data[top_viper_hubs]
      colors <- rep("black", length(barplot_data))
      names(colors) <- names(barplot_data)
      colors[top_viper_hubs] <- "red"
      png(paste0(params[[5]], aracne_name, "_top_", params[[4]], "_tvhs_cosmic_enrichment_pvals_RS.png"),
          width = 2000, height = 1400, res = 130)
      barplot(barplot_data, col = colors, border = colors,
              main = paste(toupper(aracne_name), "Cosmic Cancer Gene Census Enrichment"),
              xaxt = "n", ylim = c(min(barplot_data, na.rm=TRUE)*1.3, max(barplot_data, na.rm=TRUE)*1.3),
              xlab = "All the hubs in the Aracne network (Ordered by regulon size)",
              ylab = "log10(Enrichment p-values)")
      legend("topright", legend = c(paste0("Top ", params[[4]], " Viper Hubs"), "Others"), col=c("red", "black"), lty = 1)
      dev.off()
    }
    
  }
  
  # ******************* which = cosmic_hubs_enrichment_with_exclusivity_score *******************
  # This function examines exclusivity scores of the top cosmic enriched hubs.
  # The top cosmic enriched hubs are determined based on statistical significance (p-value)
  # that how many cosmic cancer genes are in regulons of a hub.
  # The results can tell whether the regulons enriched with cancer genes also have
  # large exclusivity scores or not.
  #
  # params[[1]]: The file path of the Aracne RDA file (All_62_ARACNE.rda)
  #              (a character vector of length 1)
  # params[[2]]: The file path of Cosmic Cancer Gene Census
  #              (a character vector of length 1)
  # params[[3]]: The number of top cosmic enriched hubs that will be tested
  #              (a number)
  # params[[4]]: The number of permutation test for computing p-value
  #              (an integer)
  # params[[5]]: The file path of the GTEx-TCGA tissue mapping info (GTEx_TCGA_Map.rda)
  #              (a character vector of length 1)
  # params[[6]]: Character string, specifying which enrichment score computations option
  #		           to use when computing exclusivity scores. Should be either "FET_SUM" or "FET_COUNT".
  #              For more details, please refer to oneOffs(which="generate_exclusivity_scores").
  #              (a character vector of length 1)
  # params[[7]]: Directory path for the results
  #              (a character vector of length 1)
  #
  # e.g., params=list("//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/All_62_ARACNE.rda",
  #                   "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/Cosmic/Cosmic_Census_100419_all.tsv",
  #                   100, 10000,
  #                   "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/GTEx_TCGA_Map.rda",
  #                   "FET_SUM",
  #                   "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/exclusive_conservation/ECH/Cosmic/cosmic_hubs_enrichment_with_exclusivity_score/")
  # e.g., params=list("./data/RDA_Files/All_62_ARACNE.rda",
  #                   "./data/Cosmic/Cosmic_Census_100419_all.tsv",
  #                   100, 10000,
  #                   "./data/RDA_Files/GTEx_TCGA_Map.rda",
  #                   "FET_SUM",
  #                   "./results/exclusive_conservation/ECH/Cosmic/cosmic_hubs_enrichment_with_exclusivity_score/")
  
  if (which == "cosmic_hubs_enrichment_with_exclusivity_score") {
    
    ### argument checking
    assertString(params[[1]])
    assertString(params[[2]])
    assertNumeric(params[[3]])
    assertIntegerish(params[[4]])
    assertString(params[[5]])
    assertChoice(params[[6]], c("FET_SUM", "FET_COUNT"))
    assertString(params[[7]])
    
    ### load library
    if(!require(xlsx, quietly = TRUE)) {
      install.packages("xlsx")
      require(xlsx, quietly = TRUE)
    }
    
    ### load data
    load(params[[1]], envir = globalenv())
    cgc <- read.table(file = params[[2]], header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE, check.names = FALSE)
    load(params[[5]])
    Sys.sleep(3)
    
    ### compute exclusivity scores of all the hubs
    if(params[[6]] == "FET_SUM") {
      oneOffs("generate_exclusivity_scores", params = list("FET_SUM"))
    } else {
      oneOffs("generate_exclusivity_scores", params = list("FET_COUNT", 1e-20))
    }
    
    ### make an empty list for saving top cosmic enriched regulons
    cosmic_enriched_regulons <- vector("list", length = length(varNames))
    names(cosmic_enriched_regulons) <- varNames
    
    ### get TCGA Aracne names
    tcga_aracne_names <- varNames[grep("tcga", varNames)]
    
    ### perform analysis for each TCGA tissue
    for(aracne_name in tcga_aracne_names) {
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
      top_cosmic_hubs <- intersect(names(all_enrichment_pvs),
                                   names(reg_exclusivity_scores))[1:params[[3]]]
      
      ### save the cosmic enriched regulons
      cosmic_enriched_regulons[[aracne_name]] <- aracne[[2]][top_cosmic_hubs]
      
      ### permutation test for getting a p-value (1-tail)
      set.seed(1234)
      permu_result <- sapply(1:(params[[4]]-1), function(x) {
        total_hubs <- intersect(names(all_enrichment_pvs), names(reg_exclusivity_scores))
        random_hubs <- total_hubs[sample(length(all_enrichment_pvs), params[[3]])]
        return(sum(reg_exclusivity_scores[random_hubs]))
      })
      pVal <- (length(which(permu_result > sum(reg_exclusivity_scores[top_cosmic_hubs])))+1) / params[[4]]
      
      ### draw a plot
      barplot_data <- -reg_exclusivity_scores
      barplot_data[top_cosmic_hubs] <- -barplot_data[top_cosmic_hubs]
      colors <- rep("black", length(barplot_data))
      names(colors) <- names(barplot_data)
      colors[top_cosmic_hubs] <- "red"
      png(paste0(params[[7]], aracne_name, "_top_", params[[3]], "_cosmic_hubs_exclusivity_scores.png"),
          width = 2000, height = 1400, res = 130)
      barplot(barplot_data, col = colors, border = colors,
              main = paste(toupper(aracne_name), "Cosmic Enriched Hubs on Exclusivity Scores", "\n",
                           params[[4]], " Permutation P-Value = ", pVal),
              xaxt = "n", ylim = c(min(barplot_data, na.rm=TRUE)*1.3, max(barplot_data, na.rm=TRUE)*1.3),
              xlab = "All the hubs in the Aracne network",
              ylab = "Regulon Exclusivity Scores")
      legend("topright", legend = c(paste0("Top ", params[[3]], " Cosmic-enriched Hubs"), "Others"), col=c("red", "black"), lty = 1)
      dev.off()
    }
    
    ### get GTEx Aracne names
    gtex_aracne_names <- varNames[!grepl("tcga", varNames)]
    
    ### perform analysis for each GTEx tissue
    for(aracne_name in gtex_aracne_names) {
      ### get Aracne network
      aracne <- get(aracne_name)
      
      ### get total genes in the network
      interactome_total_genes <- as.character(getInteractomeGenes(aracne_name, count = FALSE))
      
      ### get cancer genes in the network
      interactome_cancer_genes <- intersect(interactome_total_genes, as.character(cgc$`Entrez GeneId`))
      
      ### compute cosmic enrichment (Fisher's exact test p-values) for all the hubs in the GTEx Aracne network
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
      top_cosmic_hubs <- intersect(names(all_enrichment_pvs),
                                   names(reg_exclusivity_scores))[1:params[[3]]]
      
      ### save the cosmic enriched regulons
      cosmic_enriched_regulons[[aracne_name]] <- aracne[[2]][top_cosmic_hubs]
      
      ### permutation test for getting a p-value (1-tail)
      set.seed(1234)
      permu_result <- sapply(1:(params[[4]]-1), function(x) {
        total_hubs <- intersect(names(all_enrichment_pvs), names(reg_exclusivity_scores))
        random_hubs <- total_hubs[sample(length(all_enrichment_pvs), params[[3]])]
        return(sum(reg_exclusivity_scores[random_hubs]))
      })
      pVal <- (length(which(permu_result > sum(reg_exclusivity_scores[top_cosmic_hubs])))+1) / params[[4]]
      
      ### draw a plot
      barplot_data <- -reg_exclusivity_scores
      barplot_data[top_cosmic_hubs] <- -barplot_data[top_cosmic_hubs]
      colors <- rep("black", length(barplot_data))
      names(colors) <- names(barplot_data)
      colors[top_cosmic_hubs] <- "red"
      png(paste0(params[[7]], aracne_name, "_top_", params[[3]], "_cosmic_hubs_exclusivity_scores.png"),
          width = 2000, height = 1400, res = 130)
      barplot(barplot_data, col = colors, border = colors,
              main = paste(toupper(aracne_name), "Cosmic Enriched Hubs on Exclusivity Scores", "\n",
                           params[[4]], " Permutation P-Value = ", pVal),
              xaxt = "n", ylim = c(min(barplot_data, na.rm=TRUE)*1.3, max(barplot_data, na.rm=TRUE)*1.3),
              xlab = "All the hubs in the Aracne network",
              ylab = "Regulon Exclusivity Scores")
      legend("topright", legend = c(paste0("Top ", params[[3]], " Cosmic-enriched Hubs"), "Others"), col=c("red", "black"), lty = 1)
      dev.off()
    }
    
    ### compare comsmic-enriched hubs between GTEx and TCGA
    common_cosmic_enriched_hubs <- data.frame(Common_Count=rep(NA, nrow(GTEx_TCGA_Map)), Common_Hubs=NA, Common_Regulon_Num=NA,
                                              stringsAsFactors = FALSE, check.names = FALSE)
    rownames(common_cosmic_enriched_hubs) <- paste0("GTEx_", GTEx_TCGA_Map[,1], "_TCGA_", toupper(GTEx_TCGA_Map[,2]))
    for(i in 1:nrow(GTEx_TCGA_Map)) {
      ### get regulon info of GTEx and TCGA
      gtex_regulon <- cosmic_enriched_regulons[[GTEx_TCGA_Map[i,"GTEx"]]]
      tcga_regulon <- cosmic_enriched_regulons[[paste0("tcga_", GTEx_TCGA_Map[i,"TCGA"])]]
      
      ### get common hubs
      common_hubs <- intersect(names(gtex_regulon), names(tcga_regulon))
      union_hubs <- union(names(gtex_regulon), names(tcga_regulon))
      common_cosmic_enriched_hubs[i,"Common_Count"] <- paste0(length(common_hubs), "/", length(union_hubs))
      common_cosmic_enriched_hubs[i,"Common_Hubs"] <- paste(common_hubs, collapse = ";")
      
      ### number of common target genes of the common hubs
      common_cosmic_enriched_hubs[i,"Common_Regulon_Num"] <- paste(sapply(common_hubs, function(x) {
        return(paste0(length(intersect(rownames(gtex_regulon[[x]]), rownames(tcga_regulon[[x]]))),
                      "/", length(union(rownames(gtex_regulon[[x]]), rownames(tcga_regulon[[x]])))))
      }), collapse = ";")
    }
    write.xlsx2(data.frame(Tissue=rownames(common_cosmic_enriched_hubs), common_cosmic_enriched_hubs,
                           stringsAsFactors = FALSE, check.names = FALSE),
                file = paste0(params[[7]], "top_", params[[3]], "_cosmic_hubs_common_between_GTex_vs_TCGA.xlsx"), row.names = FALSE)
    
    ### perform analysis with random Aracne networks
    set.seed(1234)
    for(aracne_name in tcga_aracne_names) {
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
        target_genes <- interactome_total_genes[sample(length(interactome_total_genes), nrow(aracne[[2]][[hub]]))]
        
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
      top_cosmic_hubs <- intersect(names(all_enrichment_pvs),
                                   names(reg_exclusivity_scores))[1:params[[3]]]
      
      ### permutation test for getting a p-value (1-tail)
      set.seed(1234)
      permu_result <- sapply(1:(params[[4]]-1), function(x) {
        total_hubs <- intersect(names(all_enrichment_pvs), names(reg_exclusivity_scores))
        random_hubs <- total_hubs[sample(length(all_enrichment_pvs), params[[3]])]
        return(sum(reg_exclusivity_scores[random_hubs]))
      })
      pVal <- (length(which(permu_result > sum(reg_exclusivity_scores[top_cosmic_hubs])))+1) / params[[4]]
      
      ### draw a plot
      barplot_data <- -reg_exclusivity_scores
      barplot_data[top_cosmic_hubs] <- -barplot_data[top_cosmic_hubs]
      colors <- rep("black", length(barplot_data))
      names(colors) <- names(barplot_data)
      colors[top_cosmic_hubs] <- "red"
      png(paste0(params[[7]], "random_", aracne_name, "_top_", params[[3]], "_cosmic_hubs_exclusivity_scores.png"),
          width = 2000, height = 1400, res = 130)
      barplot(barplot_data, col = colors, border = colors,
              main = paste("Random", toupper(aracne_name), "Cosmic Enriched Hubs on Exclusivity Scores", "\n",
                           params[[4]], " Permutation P-Value = ", pVal),
              xaxt = "n", ylim = c(min(barplot_data, na.rm=TRUE)*1.3, max(barplot_data, na.rm=TRUE)*1.3),
              xlab = "All the hubs in the Aracne network",
              ylab = "Regulon Exclusivity Scores")
      legend("topright", legend = c(paste0("Top ", params[[3]], " Cosmic-enriched Hubs"), "Others"), col=c("red", "black"), lty = 1)
      dev.off()
    }
    
  }
  
  # ******************* which = highly_mutated_hubs_enrichment_with_exclusivity_score *******************
  # This function examines exclusivity scores of the top highly mutated hubs.
  # The top highly mutated hubs are computed by mutation counts per gene normalized by gene length.
  # The results can tell whether highly mutated hubs also have large exclusivity scores or not.
  #
  # params[[1]]: The file path of the Aracne RDA file (All_62_ARACNE.rda)
  #              (a character vector of length 1)
  # params[[2]]: The file path of highly mutated genes list of TCGA tissues (TCGA_33_Top_Mutated_Genes.rda)
  #              (a character vector of length 1)
  # params[[3]]: The number of top highly mutated hubs that will be tested
  #              (a number)
  # params[[4]]: The number of permutation test for computing p-value
  #              (an integer)
  # params[[5]]: Character string, specifying which enrichment score computations option
  #		           to use when computing exclusivity scores. Should be either "FET_SUM" or "FET_COUNT".
  #              For more details, please refer to oneOffs(which="generate_exclusivity_scores").
  #              (a character vector of length 1)
  # params[[6]]: Directory path for the results
  #              (a character vector of length 1)
  #
  # e.g., params=list("//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/All_62_ARACNE.rda",
  #                   "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/TCGA_33_Top_Mutated_Genes.rda",
  #                   100, 10000,
  #                   "FET_SUM",
  #                   "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/exclusive_conservation/ECH/reg_exclusivity/highly_mutated_hubs_enrichment_with_exclusivity_score/")
  # e.g., params=list("./data/RDA_Files/All_62_ARACNE.rda",
  #                   "./data/RDA_Files/TCGA_33_Top_Mutated_Genes.rda",
  #                   100, 10000,
  #                   "FET_SUM",
  #                   "./results/exclusive_conservation/ECH/reg_exclusivity/highly_mutated_hubs_enrichment_with_exclusivity_score/")
  
  if (which == "highly_mutated_hubs_enrichment_with_exclusivity_score") {
    
    ### argument checking
    assertString(params[[1]])
    assertString(params[[2]])
    assertNumeric(params[[3]])
    assertIntegerish(params[[4]])
    assertChoice(params[[5]], c("FET_SUM", "FET_COUNT"))
    assertString(params[[6]])
    
    ### load data
    load(params[[1]], envir = globalenv())
    load(params[[2]])
    Sys.sleep(3)
    
    ### compute exclusivity scores of all the hubs
    if(params[[5]] == "FET_SUM") {
      oneOffs("generate_exclusivity_scores", params = list("FET_SUM"))
    } else {
      oneOffs("generate_exclusivity_scores", params = list("FET_COUNT", 1e-20))
    }
    
    ### get TCGA Aracne names
    tcga_aracne_names <- varNames[grep("tcga", varNames)]
    
    ### perform analysis for each TCGA tissue
    for(aracne_name in tcga_aracne_names) {
      ### get Aracne network
      aracne <- get(aracne_name)
      
      ### get highly mutated hubs for the tissue
      mutated_hubs <- geneSymbolToEntrezId(names(tcga_highly_mutated_genes[[toupper(aracne_name)]]))
      mutated_hubs <- intersect(intersect(as.character(mutated_hubs), rownames(aracne[[1]])),
                                names(reg_exclusivity_scores))
      top_mutated_hubs <- mutated_hubs[1:params[[3]]]
      
      ### permutation test for getting a p-value (1-tail)
      set.seed(1234)
      permu_result <- sapply(1:(params[[4]]-1), function(x) {
        random_hubs <- mutated_hubs[sample(length(mutated_hubs), params[[3]])]
        return(sum(reg_exclusivity_scores[random_hubs]))
      })
      pVal <- (length(which(permu_result > sum(reg_exclusivity_scores[top_mutated_hubs])))+1) / params[[4]]
      
      ### draw a plot
      barplot_data <- -reg_exclusivity_scores
      barplot_data[top_mutated_hubs] <- -barplot_data[top_mutated_hubs]
      colors <- rep("black", length(barplot_data))
      names(colors) <- names(barplot_data)
      colors[top_mutated_hubs] <- "red"
      png(paste0(params[[6]], aracne_name, "_top_", params[[3]], "_top_mutated_hubs_exclusivity_scores.png"),
          width = 2000, height = 1400, res = 130)
      barplot(barplot_data, col = colors, border = colors,
              main = paste(toupper(aracne_name), "Top Mutated Hubs on Exclusivity Scores", "\n",
                           params[[4]], " Permutation P-Value = ", pVal),
              xaxt = "n", ylim = c(min(barplot_data, na.rm=TRUE)*1.3, max(barplot_data, na.rm=TRUE)*1.3),
              xlab = "All the hubs in the Aracne network",
              ylab = "Regulon Exclusivity Scores")
      legend("topright", legend = c(paste0("Top ", params[[3]], " Top-mutated Hubs"), "Others"), col=c("red", "black"), lty = 1)
      dev.off()
    }
    
  }
  
  # ******************* which = highly_mutated_regulons_enrichment_with_exclusivity_score *******************
  # This function examines exclusivity scores of the top highly mutated hubs.
  # The top highly mutated hubs are computed by mutation counts per gene normalized by gene length.
  # The results can tell whether highly mutated hubs also have large exclusivity scores or not.
  #
  # params[[1]]: The file path of the Aracne RDA file (All_62_ARACNE.rda)
  #              (a character vector of length 1)
  # params[[2]]: The file path of highly mutated regulon list of TCGA tissues (TCGA_26_Top_Mutated_Regulons.rda)
  #              (a character vector of length 1)
  # params[[3]]: The number of top highly mutated regulons that will be tested
  #              (a number)
  # params[[4]]: The number of permutation test for computing p-value
  #              (an integer)
  # params[[5]]: Character string, specifying which enrichment score computations option
  #		           to use when computing exclusivity scores. Should be either "FET_SUM" or "FET_COUNT".
  #              For more details, please refer to oneOffs(which="generate_exclusivity_scores").
  #              (a character vector of length 1)
  # params[[6]]: Directory path for the results
  #              (a character vector of length 1)
  #
  # e.g., params=list("//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/All_62_ARACNE.rda",
  #                   "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/TCGA_26_Top_Mutated_Regulons.rda",
  #                   100, 10000,
  #                   "FET_SUM",
  #                   "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/exclusive_conservation/ECH/reg_exclusivity/highly_mutated_regulons_enrichment_with_exclusivity_score/")
  # e.g., params=list("./data/RDA_Files/All_62_ARACNE.rda",
  #                   "./data/RDA_Files/TCGA_26_Top_Mutated_Regulons.rda",
  #                   100, 10000,
  #                   "FET_SUM",
  #                   "./results/exclusive_conservation/ECH/reg_exclusivity/highly_mutated_regulons_enrichment_with_exclusivity_score/")
  
  if (which == "highly_mutated_regulons_enrichment_with_exclusivity_score") {
    
    ### argument checking
    assertString(params[[1]])
    assertString(params[[2]])
    assertNumeric(params[[3]])
    assertIntegerish(params[[4]])
    assertChoice(params[[5]], c("FET_SUM", "FET_COUNT"))
    assertString(params[[6]])
    
    ### load data
    load(params[[1]], envir = globalenv())
    load(params[[2]])
    Sys.sleep(3)
    
    ### compute exclusivity scores of all the hubs
    if(params[[5]] == "FET_SUM") {
      oneOffs("generate_exclusivity_scores", params = list("FET_SUM"))
    } else {
      oneOffs("generate_exclusivity_scores", params = list("FET_COUNT", 1e-20))
    }
    
    ### get TCGA Aracne names
    tcga_aracne_names <- varNames[grep("tcga", varNames)]
    
    ### perform analysis for each TCGA tissue
    for(aracne_name in tcga_aracne_names) {
      ### get highly mutated regulons for the tissue
      mutated_regulons <- intersect(rownames(tcga_highly_mutated_regulons[[aracne_name]]),
                                    names(reg_exclusivity_scores))
      top_mutated_regulons <- mutated_regulons[1:params[[3]]]
      
      ### permutation test for getting a p-value (1-tail) (top mutated hubs on exclusivity scores)
      set.seed(1234)
      permu_result <- sapply(1:(params[[4]]-1), function(x) {
        random_hubs <- mutated_regulons[sample(length(mutated_regulons), params[[3]])]
        return(sum(reg_exclusivity_scores[random_hubs]))
      })
      pVal <- (length(which(permu_result > sum(reg_exclusivity_scores[top_mutated_regulons])))+1) / params[[4]]
      
      ### draw a plot (top mutated hubs on exclusivity scores)
      barplot_data <- -reg_exclusivity_scores
      barplot_data[top_mutated_regulons] <- -barplot_data[top_mutated_regulons]
      colors <- rep("black", length(barplot_data))
      names(colors) <- names(barplot_data)
      colors[top_mutated_regulons] <- "red"
      png(paste0(params[[6]], aracne_name, "_top_", params[[3]], "_mutated_regulons_exclusivity_scores.png"),
          width = 2000, height = 1400, res = 130)
      barplot(barplot_data, col = colors, border = colors,
              main = paste(toupper(aracne_name), "Top Mutated Regulons on Exclusivity Scores", "\n",
                           params[[4]], " Permutation P-Value = ", pVal),
              xaxt = "n", ylim = c(min(barplot_data, na.rm=TRUE)*1.3, max(barplot_data, na.rm=TRUE)*1.3),
              xlab = "All the Hubs in the Aracne Network",
              ylab = "Regulon Exclusivity Scores")
      legend("topright", legend = c(paste0("Top ", params[[3]], " Top-mutated Regulons"), "Others"), col=c("red", "black"), lty = 1)
      dev.off()
      
      ### get top exclusively conserved hubs
      top_echs <- intersect(names(reg_exclusivity_scores),
                            rownames(tcga_highly_mutated_regulons[[aracne_name]]))[1:params[[3]]]
      
      ### permutation test for getting a p-value (1-tail) (top exclusively conserved hubs on regulon mutation statistics)
      set.seed(1234)
      permu_result <- sapply(1:(params[[4]]-1), function(x) {
        total_hubs <- intersect(names(reg_exclusivity_scores), rownames(tcga_highly_mutated_regulons[[aracne_name]]))
        random_hubs <- total_hubs[sample(length(total_hubs), params[[3]])]
        return(sum(tcga_highly_mutated_regulons[[aracne_name]][random_hubs,"S"]))
      })
      pVal <- (length(which(permu_result > sum(tcga_highly_mutated_regulons[[aracne_name]][top_echs,"S"])))+1) / params[[4]]
      
      ### draw a plot (top exclusively conserved hubs on regulon mutation statistics)
      barplot_data <- -tcga_highly_mutated_regulons[[aracne_name]][,"S"]
      names(barplot_data) <- rownames(tcga_highly_mutated_regulons[[aracne_name]])
      barplot_data[top_echs] <- -barplot_data[top_echs]
      colors <- rep("black", length(barplot_data))
      names(colors) <- names(barplot_data)
      colors[top_echs] <- "red"
      png(paste0(params[[6]], aracne_name, "_top_", params[[3]], "_ECHs_Mutation_Enrichment.png"),
          width = 2000, height = 1400, res = 130)
      barplot(barplot_data, col = colors, border = colors,
              main = paste(toupper(aracne_name), "Top ECHs on Mutation Enrichment", "\n",
                           params[[4]], " Permutation P-Value = ", pVal),
              xaxt = "n", ylim = c(min(barplot_data, na.rm=TRUE)*1.3, max(barplot_data, na.rm=TRUE)*1.3),
              xlab = "All the Hubs in the Aracne Network",
              ylab = "Mutation Enrichment Chi-square Statistics")
      legend("topright", legend = c(paste0("Top ", params[[3]], " Top-mutated Regulons"), "Others"), col=c("red", "black"), lty = 1)
      dev.off()
    }
    
  }
  
  # ******************* which = exclusivity_rank_comparison *******************
  # Hubs are given ranks based on exclusivity scores/Cosmic enrichment/mutation enrichment.
  # We would like to know where the top hubs from each are located in other's rank range.
  # And also want to know correlations among them.
  #
  #
  # params[[1]]: The file path of the Aracne RDA file (All_62_ARACNE.rda)
  #              (a character vector of length 1)
  # params[[2]]: The file path of Cosmic-enriched regulon list of All GTEx and TCGA tissues (All_62_Top_Cosmic_Regulons.rda)
  #              (a character vector of length 1)
  # params[[3]]: The file path of highly mutated regulon list of TCGA tissues (TCGA_26_Top_Mutated_Regulons.rda)
  #              (a character vector of length 1)
  # params[[4]]: The number of top hubs that will be tested
  #              (a number)
  # params[[5]]: The number of permutation test for computing p-value
  #              (an integer)
  # params[[6]]: Character string, specifying which enrichment score computations option
  #		           to use when computing exclusivity scores. Should be either "FET_SUM" or "FET_COUNT".
  #              For more details, please refer to oneOffs(which="generate_exclusivity_scores").
  #              (a character vector of length 1)
  # params[[7]]: Directory path for the results
  #              (a character vector of length 1)
  #
  # e.g., params=list("//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/All_62_ARACNE.rda",
  #                   "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/All_62_Top_Cosmic_Regulons.rda",
  #                   "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/TCGA_26_Top_Mutated_Regulons.rda",
  #                   100, 10000,
  #                   "FET_SUM",
  #                   "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/exclusive_conservation/ECH/reg_exclusivity/rank_comparison/")
  # e.g., params=list("./data/RDA_Files/All_62_ARACNE.rda",
  #                   "./data/RDA_Files/All_62_Top_Cosmic_Regulons.rda",
  #                   "./data/RDA_Files/TCGA_26_Top_Mutated_Regulons.rda",
  #                   100, 10000,
  #                   "FET_SUM",
  #                   "./results/exclusive_conservation/ECH/reg_exclusivity/rank_comparison/")
  
  if (which == "exclusivity_rank_comparison") {
    
    ### argument checking
    assertString(params[[1]])
    assertString(params[[2]])
    assertString(params[[3]])
    assertNumeric(params[[4]])
    assertIntegerish(params[[5]])
    assertChoice(params[[6]], c("FET_SUM", "FET_COUNT"))
    assertString(params[[7]])
    
    ### load data
    load(params[[1]], envir = globalenv())
    load(params[[2]])
    load(params[[3]])
    Sys.sleep(3)
    
    ### compute exclusivity scores of all the hubs
    if(params[[6]] == "FET_SUM") {
      oneOffs("generate_exclusivity_scores", params = list("FET_SUM"))
    } else {
      oneOffs("generate_exclusivity_scores", params = list("FET_COUNT", 1e-20))
    }
    
    ### rank comparison between exclusivity score and Cosmic enrichment
    for(tissue_name in names(cosmic_enriched_regulons)) {
      ### get common hubs between exclusivity score and Cosmic enrichment
      common_hubs <- intersect(names(reg_exclusivity_scores), rownames(cosmic_enriched_regulons[[tissue_name]]))
      
      ### prepare two vectors for comparison
      A <- reg_exclusivity_scores[common_hubs]
      B <- -log10(cosmic_enriched_regulons[[tissue_name]][common_hubs,"FDR"])
      names(B) <- rownames(cosmic_enriched_regulons[[tissue_name]][common_hubs,])
      
      ### print plot of the rank comparison between exclusivity and cosmic enrichment
      compare_two_different_ranks(A = A,
                                  B = B,
                                  A_name = "Exclusivity_Score",
                                  B_name = "-log10(Cosmic_Enrichment_FDR)",
                                  ordering = "decreasing",
                                  top = params[[4]],
                                  alternative = "greater",
                                  permutation = 10000,
                                  fileName = paste0("Exclusivity_VS_Cosmic_Enrichment_", tissue_name),
                                  printPath = params[[7]],
                                  width = 24,
                                  height = 10)
    }
    
    ### rank comparison between exclusivity score and mutation enrichment
    for(tissue_name in names(tcga_highly_mutated_regulons)) {
      ### get common hubs between exclusivity score and mutation enrichment
      common_hubs <- intersect(names(reg_exclusivity_scores), rownames(tcga_highly_mutated_regulons[[tissue_name]]))
      
      ### prepare two vectors for comparison
      A <- reg_exclusivity_scores[common_hubs]
      B <- tcga_highly_mutated_regulons[[tissue_name]][common_hubs,"S"]
      names(B) <- rownames(tcga_highly_mutated_regulons[[tissue_name]][common_hubs,])
      
      ### print plot of the rank comparison between exclusivity and cosmic enrichment
      compare_two_different_ranks(A = A,
                                  B = B,
                                  A_name = "Exclusivity_Score",
                                  B_name = "Mutation_Enrichment_Statistics",
                                  ordering = "decreasing",
                                  top = params[[4]],
                                  alternative = "greater",
                                  permutation = 10000,
                                  fileName = paste0("Exclusivity_VS_Mutation_Enrichment_", tissue_name),
                                  printPath = params[[7]],
                                  width = 24,
                                  height = 10)
    }
    
  }
  
  # ******************* which = differentially_conserved_hubs *******************
  # This function identifies hubs that are conserved well in both GTEx and TCGA
  # but their target genes are very different between the two.
  #
  # params[[1]]: The file path of the Aracne RDA file (All_62_ARACNE.rda)
  #              (a character vector of length 1)
  # params[[2]]: Character string, specifying which enrichment score computations option
  #		           to use when computing exclusivity scores. Should be either "FET_SUM" or "FET_COUNT".
  #              For more details, please refer to oneOffs(which="generate_exclusivity_scores").
  #              (a character vector of length 1)
  # params[[3]]: The file path of the GTEx-TCGA tissue mapping info (GTEx_TCGA_Map.rda)
  #              (a character vector of length 1)
  # params[[4]]: Directory path for the results
  #              (a character vector of length 1)
  #
  # e.g., params=list("//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/All_62_ARACNE.rda",
  #                   "FET_SUM",
  #                   "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/GTEx_TCGA_Map.rda",
  #                   "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/exclusive_conservation/ECH/differentially_conserved_hubs/")
  # e.g., params=list("./data/RDA_Files/All_62_ARACNE.rda",
  #                   "FET_SUM",
  #                   "./data/RDA_Files/GTEx_TCGA_Map.rda",
  #                   "./results/exclusive_conservation/ECH/differentially_conserved_hubs/")
  
  if (which == "differentially_conserved_hubs") {
    
    ### argument checking
    assertString(params[[1]])
    assertChoice(params[[2]], c("FET_SUM", "FET_COUNT"))
    assertString(params[[3]])
    assertString(params[[4]])
    
    ### load library
    if(!require(xlsx, quietly = TRUE)) {
      install.packages("xlsx")
      require(xlsx, quietly = TRUE)
    }
    
    ### load data
    load(params[[1]], envir = globalenv())
    load(params[[3]])
    Sys.sleep(3)
    
    ### compute exclusivity scores of all the hubs
    if(params[[2]] == "FET_SUM") {
      oneOffs("generate_exclusivity_scores", params = list("FET_SUM"))
    } else {
      oneOffs("generate_exclusivity_scores", params = list("FET_COUNT", 1e-20))
    }
    
    ### prepare a result table
    result_table <- data.frame(reg_exclusivity_scores_collection,
                               stringsAsFactors = FALSE, check.names = FALSE)
    colnames(result_table)[3] <- "BOTH_All"
    
    ### prepare colnames that are from same tissue but different collections
    target_colnames <- apply(GTEx_TCGA_Map, 1, function(x) {
      return(paste0(x["GTEx"], "%tcga_", x["TCGA"]))
    })
    
    ### scale factors
    homog_num = choose(sum(grepl("tcga", varNames)),2)+choose(length(varNames)-sum(grepl("tcga", varNames)), 2)
    heter_num = choose(length(varNames), 2) - homog_num
    same_scale = max(homog_num, heter_num) / length(target_colnames)
    
    ### Replace -Inf with the smallest recorded non-infinite log(p-value).
    minp = min(reg_fet_cons_mat[reg_fet_cons_mat != -Inf])
    
    ### calculate the sum(FET pVal) of interactome pairs that are from the same tissue of each collection
    result_table$BOTH_Same_Tissue <- NA
    for(hub in rownames(result_table)) {
      fets = reg_fet_cons_mat[hub, ]
      fets = replace(fets, fets == -Inf, minp)
      if(params[[2]] == "FET_SUM") {
        result_table[hub, "BOTH_Same_Tissue"] <- -sum(fets[target_colnames])*same_scale
      } else {
        result_table[hub, "BOTH_Same_Tissue"] <- sum(fets[target_colnames] <= log10(thresh))*same_scale
      }
    }
    
    ### calculate exclusivity scores (the original exclusivity scores)
    ### GTEX + TCGA - BOTH_All
    result_table$Exclusivity_Score <- round(result_table$GTEX + result_table$TCGA - result_table$BOTH_All)
    
    ### calculate the differential exclusivity scores
    ### conserved in the same collections but the hubs have very different targets
    ### in the same tissues between GTEx and TCGA
    ### GTEX + TCGA - BOTH_Same_Tissue
    result_table$Differential_Exclusivity_Score <- round(result_table$GTEX + result_table$TCGA - result_table$BOTH_Same_Tissue)
    
    ### there are some hubs only exist in GTEx Aracne networks and
    ### in only some of TCGA networks. Get hub names that exist in
    ### all the GTEx and TCGA tissues, and only use them.
    all_intersect_hubs = as.character(getInteractomeGenes(nets=varNames, count=FALSE, hubs_only=TRUE, common=TRUE))
    result_table <- result_table[all_intersect_hubs,]
    
    ### order the hubs based on the differential exclusivity scores
    result_table <- result_table[order(-result_table$Differential_Exclusivity_Score),]
    
    ### add gene symbol to the result_table
    result_table <- data.frame(Entrez_ID=rownames(result_table),
                               Gene_Symbol=entrezIDtoSymbol(rownames(result_table)),
                               result_table,
                               stringsAsFactors = FALSE, check.names = FALSE)
    
    write.xlsx2(result_table, file = paste0(params[[4]], "Differentially_Conserved_Hubs_Statics_Table.xlsx"), row.names = FALSE)
    
    ### with a given hub list and threshold (0.01), make a pathway table that shows
    ### which pathways are most enriched with the regulons of all the given hubs
    ### implement it in aracne_go.R
    
    
  }
  
  # ******************* which = network_similarity *******************
  # This function calculates similarities among all the existing iteractomes
  # pairwisely with using graph kernels.
  #
  # params[[1]]: The file path of the Aracne RDA file (All_62_ARACNE.rda)
  #              (a character vector of length 1)
  # params[[2]]: Graph kernel to be used ("GeometricRandomWalk", "Graphlet", "ShortestPath")
  #              (a character vector of length 1)
  # params[[3]]: The file path of the interactome similarity result
  #              (a character vector of length 1)
  #
  # e.g., params=list("//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/All_62_ARACNE.rda",
  #                   "GeometricRandomWalk",
  #                   "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/All_62_Similarity.rda")
  # e.g., params=list("./data/RDA_Files/All_62_ARACNE.rda",
  #                   "GeometricRandomWalk",
  #                   "./data/RDA_Files/All_62_Similarity.rda")
  
  if (which == "network_similarity") {
    
    ### argument checking
    assertString(params[[1]])
    assertChoice(params[[2]], c("GeometricRandomWalk", "Graphlet", "ShortestPath"))
    assertString(params[[3]])
    
    ### load library
    if(!require(graphkernels, quietly = TRUE)) {
      install.packages("graphkernels")
      require(graphkernels, quietly = TRUE)
    }
    if(!require(igraph, quietly = TRUE)) {
      install.packages("igraph")
      require(igraph, quietly = TRUE)
    }
    
    ### load data
    load(params[[1]], envir = globalenv())
    Sys.sleep(3)
    
    ### transform the interactome networks to igraphs
    igs <- vector("list", length = length(varNames))
    names(igs) <- varNames
    for(net_name in varNames) {
      ### get the network
      net <- get(net_name)
      
      ### gather all interactions for the given network
      edge_list <- data.frame(Hub=names(net[[2]])[1],
                              Target=as.character(abs(net[[2]][[1]][,1])),
                              weight=net[[2]][[1]][,2],
                              stringsAsFactors = FALSE)
      for(i in 2:length(net[[2]])) {
        edge_list <- rbind(edge_list,
                           data.frame(Hub=names(net[[2]])[i],
                                      Target=as.character(abs(net[[2]][[i]][,1])),
                                      weight=net[[2]][[i]][,2],
                                      stringsAsFactors = FALSE))
      }
      
      ### transform into an igraph
      igs[[net_name]] <- graph.data.frame(edge_list, directed = FALSE)
      
      ### simplify the igraph (remove duplicate edges)
      igs[[net_name]] <- simplify(igs[[net_name]],
                                  remove.multiple = TRUE,
                                  remove.loops = TRUE,
                                  edge.attr.comb = "first")
      
      ### garbage collection
      rm(edge_list)
      rm(net_name, envir = globalenv())
      gc()
    }
    
    ### remove original Aracne networks from memory
    rm(list = c(varNames), envir = globalenv())
    gc()
    
    ### pairwisely calculate the similarities
    ### the run using the graphkernels can be vectorized,
    ### but it would take a long time, and need large memory
    ### so we run one at a time and save each step to be safe
    similarity_mat <- matrix(0, length(varNames), length(varNames))
    rownames(similarity_mat) <- varNames
    colnames(similarity_mat) <- varNames
    if(params[[2]] == "GeometricRandomWalk") {
      for(net1 in varNames[-length(varNames)]) {
        for(net2 in varNames[seq(which(varNames == net1)+1, length(varNames))]) {
          ### try the graph kernel - one pair at a time
          kernel_mat <- CalculateGeometricRandomWalkKernel(G = list(igs[[net1]], igs[[net2]]), par = 0.1)
          
          ### save the calculated similarities to the result matrix
          similarity_mat[net1, net1] <- kernel_mat[1, 1]
          similarity_mat[net1, net2] <- kernel_mat[1, 2]
          similarity_mat[net2, net1] <- kernel_mat[2, 1]
          similarity_mat[net2, net2] <- kernel_mat[2, 2]
        }
      }
    } else if(params[[2]] == "Graphlet") {
      for(net1 in varNames[-length(varNames)]) {
        for(net2 in varNames[seq(which(varNames == net1)+1, length(varNames))]) {
          ### try the graph kernel - one pair at a time
          kernel_mat <- CalculateConnectedGraphletKernel(G = list(igs[[net1]], igs[[net2]]), par = 4)
          
          ### save the calculated similarities to the result matrix
          similarity_mat[net1, net1] <- kernel_mat[1, 1]
          similarity_mat[net1, net2] <- kernel_mat[1, 2]
          similarity_mat[net2, net1] <- kernel_mat[2, 1]
          similarity_mat[net2, net2] <- kernel_mat[2, 2]
        }
      }
    } else if(params[[2]] == "ShortestPath") {
      for(net1 in varNames[-length(varNames)]) {
        for(net2 in varNames[seq(which(varNames == net1)+1, length(varNames))]) {
          ### try the graph kernel - one pair at a time
          kernel_mat <- CalculateShortestPathKernel(G = list(igs[[net1]], igs[[net2]]))
          
          ### save the calculated similarities to the result matrix
          similarity_mat[net1, net1] <- kernel_mat[1, 1]
          similarity_mat[net1, net2] <- kernel_mat[1, 2]
          similarity_mat[net2, net1] <- kernel_mat[2, 1]
          similarity_mat[net2, net2] <- kernel_mat[2, 2]
        }
      }
    } else {
      stop("ERROR: params[[2]] should be either \"GeometricRandomWalk\", \"Graphlet\", or \"ShortestPath\".")
    }
    
    ### set README function
    README <- function() {
      writeLines(paste(rep("#", 100), collapse = ""))
      writeLines("The \"similarity_mat\" is a 62 x 62 numeric matrix.")
      writeLines("It is a symmetric matrix that rows and columns both represent")
      writeLines("62 GTEx + TCGA tissues. A value in a cell indicates similarity")
      writeLines("between the two tissue Aracne networks.")
      writeLines("The matrix was generated by which = \"network_similarity\"")
      writeLines("of oneOffs() in aracne.R.")
      writeLines(paste(rep("#", 100), collapse = ""))
    }
    
    ### save as RDA
    save(list = c("similarity_mat", "README"), file = params[[3]])
    
  }
  
}


#################################################################################
# Generate a matrix with one row per hub gene and three columns titled "GTEX",
# "TCGA", and "BOTH", listing the number of interactome pairs where the hub gene 
# regulon is conserved at a FET-pvalue above a given p-value threshold. The "GTEX"
# column counts the number of pairs involving only GTEx interactomes; the column
# "TCGA" the number of pairs involving only TCGA interactomes; and the column
# "BOTH" the number of pairs involving one GTEx and one TCGA interactome.
# 
# Uses the data in tfPairEnrich[[2]]. Specifically, it first defines a p-value
# threshold as follows:
# 1. Aggregate all FET p-values from all matrices in tfPairEnrich[[2]], and order
#		them in increasing order.
# 2. Use the function argument "cutoff_ratio" to define a singnificance threshold.
#		E.g., if there are 1,000,000 p-values in the ordered list from step 1 and
#		cutoff_ratio = 0.001, the select the 1000-th p-value from the order list
#		as the p-value significance threshold T.
# 
# The threshold T computed above is used for filtering the pairwise interactome
# conservation: specifically, the pairwise counts reported in the methods results
# matrix are computed by considering only piarwise comparisons where the FET p-value
# is less than T.
#
# NOTE: Before using this method the following function call must be made, at it
# generates variables that this method utilizes:
#		oneOffs((which == "generate_exclusivity_scores")
#################################################################################
regulonConservationMode <- function(cutoff_ratio){
	pvals = sort(unlist(sapply(tfPairEnrich[[2]], function(m){return(m[, "log10_pval"])})))
	thresh = pvals[ceiling(length(pvals)*cutoff_ratio)]
	res = apply(reg_fet_cons_mat, 1, function(h_row){
				h_row = h_row[h_row <= thresh]
				if(length(h_row) == 0)
					return(NA)
				res = sapply(c("GTEX", "TCGA", "BOTH"), function(type){
							return(sum(interactome_pair_map[names(h_row)] == type))
						})
				return(res)
			})
	res = res[sapply(res, function(v){return(!is.na(v)[1])})]
	n = names(res)
	res = Reduce(rbind, res)
	rownames(res) = entrezIDtoSymbol(n)
	return(res)
} 


uniqueGenomicAlterations <- function(rootDir = "//isilon.c2b2.columbia.edu/ifs/scratch/c2b2/ac_lab/fmg2117/projects/pancancer/cptac/citrusFiles/"){
	fileNames = c("genomicEvents-blca.dat", "genomicEvents-brca.dat", "genomicEvents-coad.dat", 
			"genomicEvents-gbm.dat", "genomicEvents-hnsc.dat", "genomicEvents-kirc.dat", 
			"genomicEvents-kirp.dat", "genomicEvents-laml.dat", "genomicEvents-lgg.dat", 
			"genomicEvents-lihc.dat", "genomicEvents-luad.dat", "genomicEvents-lusc.dat", 
			"genomicEvents-ov.dat", "genomicEvents-prad.dat", "genomicEvents-read.dat", 
			"genomicEvents-sarc.dat", "genomicEvents-skcm.dat", "genomicEvents-stad.dat", 
			"genomicEvents-thca.dat", "genomicEvents-ucec.dat")
	
	res = vector(mode = "integer", length = length(fileNames))
	names(res) = sapply(fileNames, function(s){return(gsub(".dat", "", 
								gsub("genomicEvents-", "", s)))})
	
	for (i in 1:length(fileNames)){
		gFile = fileNames[i]
		logLines(paste("\nProcessing file -> ", gFile))
		data = as.matrix(read.table(paste(rootDir, gFile, sep=""), header = TRUE))
		res[i] = length(unique(data[,1]))
	}
	return(res)
}


filterDGs <- function(altsPerTumor, threshold = 0.01){
	for (index in 1:length(varNamesDG)){
		dg_data = get(varNamesDG[index])
		numRegs = length(dg_data[[1]][,1])
		select = rep(TRUE, length=numRegs)
		correctFactor = numRegs * altsPerTumor[gsub("DG", "", varNamesDG[index])]
		for (i in 1:numRegs){
			mods = dg_data[[2]][[dg_data[[1]][i,2]]]
			mods = mods[mods[,6] < (threshold/correctFactor),,drop=FALSE]
			if (length(mods) == 0)
				select[i] = FALSE
			else
				dg_data[[2]][[dg_data[[1]][i,2]]] = mods
		}
		l1 = dg_data[[1]][select,]
		l2 = list()
		# writeLines(paste("Length = ", length(l1)))
		if (length(l1) > 0){
			for (i in 1:length(l1[,1])){
				# writeLines(as.character(i))
				l2[[i]] = dg_data[[2]][[l1[i,2]]]
				l1[i,2] = i
				l1[i,3] = length(l2[[i]][,1]) 
			}
			names(l2) = l1[,1]
			l1 = l1[order(l1[,3], decreasing=TRUE),]
			rownames(l1) = l1[,1]
		}
		
		res = list()
		res[[1]] = l1
		res[[2]] = l2
		assign(paste(varNamesDG[index],"F", sep=""), res, envir = globalenv())
	}
}


###### FIXME
# A quick and dirty way to extract a table with 20 rows (one for each tunmor) 
# and columns representing TFs that are highly modulated across mulitple tumors
# where entry [i, j] contains the compbinedP p-value associated with modulator 
# modID and the j-th TF on the i-th tumor. Values are temporariry stored as 
# string to facilitate export to an Excel file.
getModInfo <- function(modID, top20CommontTFs){
	results = matrix("", length(varNamesDG), length(top20CommontTFs))
	colnames(results) = entrezIDtoSymbol(top20CommontTFs)
	rownames(results) = varNamesDG
	for (i in 1:length(top20CommontTFs)){
		l = tfDGDetails(top20CommontTFs[i], modGeneID=modID)
		for (j in 1:length(varNamesDG)){
			if (is.null(l[[j]]))
				results[i, j] = ""
			else
				results[i, j] = toString(l[[j]][1,6])
		}
	}
	return(results)
}


# *****************************************************************************
# Convert the VIPER files sent by Federico on 7/10/2016 to the "standard" 3
# column format of the origina files founf here:
#		/ifs/scratch/c2b2/ac_lab/fmg2117/projects/pancancer/cptac/citrusFiles/
# *****************************************************************************	
convertViper <- function(workingDir = "//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/CPTAC/integrated_refined_viper"){
	setwd(workingDir)
	fNames = dir()
	# Directory where to copy Federico's original files
	newDir = "Originals/"
	dir.create(newDir)
	
	for (fileName in fNames){
		writeLines(paste("--> Now processing file:", fileName))
		vipMat = as.matrix(read.table(fileName, header = TRUE, check.names = FALSE))
		file.rename(fileName, paste(newDir, fileName, sep=""))
		# Rearrange matrix, to the format we want
		res = matrix(nrow = length(vipMat), ncol = 3)
		colnames(res) = c("Gene", "Sample", "ViperValue")
		k = 1
		rowNames = rownames(vipMat)
		colNames = colnames(vipMat)
		for (i in 1:length(vipMat[,1]))
			for (j in 1:length(vipMat[1,])){
				res[k,1] = rowNames[i]
				res[k,2] = colNames[j]
				res[k,3] = vipMat[i,j]
				k = k+1
			}
		write.table(res, file=fileName, row.names = FALSE, col.names = TRUE, quote=FALSE, sep="\t")
	}
}


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


# *****************************************************************************
# Plot regulon sizes in decreasing order, with vertical lines corresponsing to
# hub type (TF = red, coTF = blue, signaling = green). To make the plots easier
# to read we are plotting the value "log(regulon_size)" inseast of "regulon_size"
#
# net:	Either a string (e.g., "brca") or an actual variable name corresponding
#		to an ARACNe network.
# *****************************************************************************	
plotReg <- function(net){
	if (is.character(net))
		net = get(net)
	colTF="red"
	colCoTF = "blue"
	colSign = "green"
	plot(log(net[[1]][,3]), t="l", main="", xlab = "Ordered hub gene index", ylab = "log(Regulon size)")
	for (i in 1:nrow(net[[1]])){
		if (is.tf(net[[1]][i, 1]))
			abline(v = i, col = colTF)
		if (is.cotf(net[[1]][i, 1]))
			abline(v = i, col = colCoTF)
		if (is.sign(net[[1]][i, 1]))
			abline(v = i, col = colSign)
	}
	lines(log(net[[1]][,3]), type = "p")
}



# *****************************************************************************
# Create summary table listing key statistics of shared eriched GO terms and 
# interactionsfor hub genes from 2 interactomes.
#
# ARGUMENTS
# * net1GO:	A named list with one entry per hub gene in the first interactome.
#			The entry for gene G is the data.frame returned by the method
#					Utils::goEnrichment(regulon(G))[[1]]
#			Only genes G for which this data.frame has at least one row should
#			be incuded.
# * net2GO:	Same as above, for the secong interactome.
#
# RETURN VALUE
# A data frame with one row per hub gene G present in both 'net1GO' and 
# 'net2GO'. Rows are named after the Entrez ID of the genes G. Columns contain:
# * The gene symbol for G.
# * The number of GO terms in net1GO enriched for the regulon of G.
# * The number of GO terms in net2GO enriched for the regulon of G .
# * The number of GO terms enriched for the regulon of G in both net1GO and net2GO.
# * The size of the regulon of G in net1.
# * The size of the regulon of G in net2.
# * The number of common targets shared by the regulos of G in net1 and net2.
# * The description of gene G.
# *****************************************************************************	
summarizeGO <- function(net1GO, net2GO){
	commonHubs = intersect(names(net1GO), names(net2GO))
	sumGO = matrix(nrow = length(commonHubs), ncol = 8)
	colnames(sumGO) = c("Gene", "net1GO", "net2GO", "SameGO", "BrcaReg", "NbrcaReg", "Common", "Description")
	for (i in 1:length(commonHubs)){
		eid = commonHubs[i]
		sumGO[i, 1] = entrezIDtoSymbol(eid)
		sumGO[i, 2] = nrow(net1GO[[eid]])
		sumGO[i, 3] = nrow(net2GO[[eid]])
		sumGO[i, 4] = length(intersect(net1GO[[eid]]$GO.ID, net2GO[[eid]]$GO.ID))
		sumGO[i, 5] = nrow(brca[[2]][[eid]])
		sumGO[i, 6] = nrow(nbrca[[2]][[eid]])
		sumGO[i, 7] = length(intersect(abs(brca[[2]][[eid]][,1]), abs(nbrca[[2]][[eid]][,1])))
		sumGO[i, 8] = entrezIDtoDescription(eid)
	}
	rownames(sumGO) = commonHubs
	sumGO = data.frame(sumGO)
	return(sumGO)
}


compareGOterms <- function(net1GO, net2GO, top=20){
	# Count how many times each GO term is found enriched. In other words, for each enriched
	# GO terms, find how many hub genes have regulons enriched for genes annotated to that 
	# term
	n1Top = sort(table(unlist(sapply(net1GO, function(x){return(x[,1])}))), decreasing = T)
	n2Top = sort(table(unlist(sapply(net2GO, function(x){return(x[,1])}))), decreasing = T)
	
	common = intersect(names(n1Top[1:top]), names(n2Top[1:top]))
	n1Only = setdiff(names(n1Top[1:top]), names(n2Top[1:top]))
	n2Only = setdiff(names(n2Top[1:top]), names(n1Top[1:top]))
	res = lapply(list(common, n1Only, n2Only), function(termList){
				res = matrix(nrow = length(termList), ncol = 3)
				colnames(res) = c("n1Count", "n2Count", "Description")
				rownames(res) = termList
				res[, 1] = n1Top[termList]
				res[, 2] = n2Top[termList]
				res[, 3] = goTermMap[termList]
				res[is.na(res)] = 0
				return(data.frame(res))
			})
	return(res)
}


# Get the list of regulators t that are modulated by gene 6794 in luad. Then retrieve
# in t1 the VIPER profiles of these regulators and draw a heatmap corresponding to 
# the Pearson's correlation of these profiles. The heatmap identifies 2 groups of 
# regulators whose activity is highly correlated. Slice out the two groups and write
# each in a CSV file, where we can copy the genes from and paste them, e.g., into
# enrichr.
# t = getModulatedRegulators(6794, luadDG)
# t1 = getViperProfiles(names(t), luadVP)
# corMat = cor(t(t1))
# L = heatmap(corMat)
# grp1 = rownames(t1)[L[[1]][1:29]]
# grp2 = rownames(t1)[L[[1]][30:61]]
# write.csv(entrezIDtoSymbol(grp1), file="tmp1.csv")
# write.csv(entrezIDtoSymbol(grp2), file="tmp2.csv")


#fun <- function(x1, x2){
#	x =  sapply(seq(1,20), function (e, A = x1, B = x2){return(setdiff(B[[e]],A[[e]]))})
#	names(x) = names(x1)
#	return(x)
#}
#x100 = getTopRegsAndMods(); x200 = getTopRegsAndMods(200);  x300 = getTopRegsAndMods(300);  x400 = getTopRegsAndMods(400)
#x100 = x100[names(eventsPerTumor)]
#x200 = x200[names(eventsPerTumor)]
#x300 = x300[names(eventsPerTumor)]
#x400 = x400[names(eventsPerTumor)]
#x2_1 = fun(x100, x200)
#x3_2 = fun(x200, x300)
#x4_3 = fun(x300, x400)
#sx100 = sapply(x100, function(e){return(paste(e, collapse = ", "))})
#sx2_1 = sapply(x2_1, function(e){return(paste(e, collapse = ", "))})
#sx3_2 = sapply(x3_2, function(e){return(paste(e, collapse = ", "))})
#sx4_3 = sapply(x4_3, function(e){return(paste(e, collapse = ", "))})
#z = rbind(sx100, sx2_1, sx3_2, sx4_3)
#rownames(z) = c("N = 100", "N = 200", "N = 300", "N = 400")
#write.xlsx(t(z), file = OUTPUT_FILE, row.names = TRUE, col.names = TRUE)


#genes = c(2099, 3169, 2625, 2305)
#xx = sapply(varNamesDG, function(e){return(which(get(e)[[1]][,1] %in% genes))})
#t <- function(e){
#	res = vector(mode="integer", length=4)
#	e = get(e)[[1]][,1]
#	j = 1
#	for (i in genes){
#		x = which(e == i)
#		if (length(x) == 0)
#			res[j] = 0
#		else
#			res[j] = x
#		j = j+1
#	}
#	return(res)
#}

# The code below identifies genes with high modulatory activity across all tunors.
#x = sapply(seq(1,20), function (e, A = x100, B = x200){return(setdiff(B[[e]],A[[e]]))})
# Return a 3 x 20 matrix listing for each tumor (column) the count of genomic events for each
# of the 3 possible types: AMP, DEL, SNV
#sapply(varNamesDG, function(x){return(rowSums(sapply(get(x)[[2]], function(e){return(c(sum(e[,2] == 1), sum(e[,2] == 2), sum(e[,2] == 3)))})))})

#tmp = as.matrix(read.table("viperValues-gbm.dat", header = TRUE))
#x = cbind(as.integer(tmp[,1]), as.numeric(tmp[,3]))
#x = x[order(abs(x[,2]), decreasing = TRUE),]   ----? Order in decreasing order of absolute activity
#plot(density(x[x[,1] == tmp[[1]][60,1] ,2]))   ---> plot density of activity measurements for the 60-th most mdulated regulator
#y = x[abs(x[,2]) > 3,] ---> Choose regulator/sample pairs with strong regulator activity.
#z = sort(table(y[,1]), decreasing = TRUE)  ---> Identify regulators with strong activity in many samples.
#intersect(strtoi(names(z[1:100])), gbmDG[[1]][1:100,1])  ---> See which of the most widely active regulators are also highly modulated

#x = sapply(ovDG[[1]][1:10, 1], function(e){v = coModulatedRegulators(e, ovDG, TRUE); v = v[v < 0.01]; return(names(v))})

#k = 1
#min = brcaDG[[2]][[1]][1,6]
#for (i in 2:length(brcaDG[[2]])){
#	if (brcaDG[[2]][[i]][1,6] < min){
#		min = brcaDG[[2]][[i]][1,6]
#		k = i
#	}
#}
#entrezIDtoSymbol(names(brcaDG[[2]])[k])

#TOP_NUM_TFS = 50
#topTFs = matrix(nrow=TOP_NUM_TFS, ncol=length(varNamesDG))
#colnames(topTFs) = varNamesDG
#for(i in 1:length(varNamesDG)) topTFs[, i] = get(varNamesDG[i])[[1]][1:TOP_NUM_TFS, 1]
#x = sort(table(as.vector(topTFs)), decreasing = TRUE)

#top20CommontTFs = sort(table(as.vector(topTFs)), decreasing = TRUE)[1:20]
#top20GlobalMods = as.integer(sapply(names(top20CommontTFs), tmp <- function(x){names(geneTableOfModulators(x))[1:20]}))
#top10Mods = as.integer(names(sort(table(top20GlobalMods), decreasing = TRUE)[1:10]))


#destination = "/data1/federico/local_tcga_data/GeneExpression/"
#source = "/ifs/scratch/c2b2/ac_lab/fmg2117/projects/pancancer/tumAcros005/"
#dirs = c("blca",  "brca",  "ccle", "coad",  "gbm",  "hnsc",  "kirc",  "kirp",  
# 	"laml",  "lgg",  "lihc",  "luad",  "lusc",  "metabric",  "ov", 	"paad",  
#	"prad",  "read",  "sarc",  "skcm",  "stad",  "taylor",  "thca",  "ucec")
#for (i in 1:length(dirs)){
#	sFile = paste(source, dirs[i], "/*expmat.rda", sep="")
#	system(paste("cp", sFile, destination))
#}

#tmp = sort(table(unlist(sapply(tfPairEnrich[[2]], function(e, TOP = 10){
#			return(e[1:TOP, 1])
#		}))), decreasing = TRUE)

# Find the regulators that regulate gene 942. This code is a sanity test,
# to check that the method getRegulators() works correctly.
#k = 0
#res = vector(mode="integer", length = 10000)
#for (i in 1:length(panCancer[[2]])){
#	ind = which(panCancer[[2]][[i]][,1] == 942)
#	if (length(ind) >0){
#		k = k+1
#		res[k] = panCancer[[1]][i,1]
#	}
#}
#res = res[1:k]

#t1 = t(sapply(tfPairEnrich[[2]], function(e){
#				ind = which(e[,1] ==  geneSymbolToEntrezId("CD86"))
#				if (length(ind) > 0)
#					return(e[ind,3:5])
#				else
#					return(NULL)
#			}))

temp <- function(gid){
	if (is.character(gid))
		gid = geneSymbolToEntrezId(gid)
	
	par(mfrow=c(4,5))
	sapply(varNames, function(var){
				vipers = get(paste(var, "VP", sep=""))[[1]]
				vipers = vipers[vipers[, 1] == gid,2]
				plot(density(vipers), main=var)
			})
	
	x = panCancer[[2]][[which(panCancer[[1]][,1] == geneSymbolToEntrezId("MYSM1"))]]
	y = sapply(unique(x[, 1]), function(gid){return(length(which(x[,1] == gid)))})
	names(y) = unique(x[, 1])
	y = sort(y, decreasing=TRUE)
	
	sapply(t, function(geneId){
				return(sapply(varNames, function(net){
									return(which(get(net)[[1]][,1] == geneId))
								}))})
	
	t = sapply(normal[[2]], function(e){y=e[,1]; return(length(y[y>0]))})
	
}

#This function clusters the GTEx Aracne networks by the top 100 TRs (i.e. those with the highest degree)  and as a null hypothesis, random 100 TRs. 
clusterGTexbytop100MRs <- function(){
	source("~/jb3401/scripts/Utils/rFunctions.R")
	load("/ifs/home/c2b2/ac_lab/jb3401/archive-af_lab/GTEx/Aracne/GTEx_Aracne.rda")
	#get the names of all the tissues
	tissue_names=varNames
	df=data.frame(matrix(NA, nrow = length(tissue_names), ncol = length(tissue_names)))
	rownames(df)=tissue_names
	colnames(df)=tissue_names
	#get all top 100 master regulators for each 36 tissues and count how many overlap in both
	for (i in (1:length(tissue_names))) {
		tissue_1=eval(as.name(tissue_names[i]))
		top_100_tissue_1=sort(tissue_1[[1]][,3],decreasing = TRUE)[0:100]
		for (j in (1:length(tissue_names))) {
			tissue_2=eval(as.name(tissue_names[j]))
			top_100_tissue_2=sort(tissue_2[[1]][,3],decreasing = TRUE)[0:100]
			#get size of intersection between 100 TRs in tissue_1 and 100 TRs in tissue_2 
			size_of_intersection=length(intersect(names(top_100_tissue_1),names(top_100_tissue_2)))
			df[i,j]=size_of_intersection
			#print(size_of_intersection)
			
		}
	}
	#cluster and create dendrogram
	hc=hclust(dist(df));dend=as.dendrogram(hc);y=labels(dend);y[]='black';y=cutree(hc,1);plot(hc)
	#put pdf in plot directory
	myPlot(hc,'dendrogram_top_100_MRs.pdf')
	dev.off()
	write.table(df,file = "Overlap_100_MR.tsv",sep = "\t",quote = FALSE)
	
	df_orig=df
	#now do randomly selected candidate transcriptional regulators
	df=data.frame(matrix(NA, nrow = length(tissue_names), ncol = length(tissue_names)))
	rownames(df)=tissue_names
	colnames(df)=tissue_names
	for (i in (1:length(tissue_names))) {
		#sample take 100 random ones
		tissue_1=eval(as.name(tissue_names[i]))
		top_100_tissue_1=sample(sort(tissue_1[[1]][,3],decreasing = TRUE))[0:100]
		for (j in (1:length(tissue_names))) {
			tissue_2=eval(as.name(tissue_names[j]))
			top_100_tissue_2=sample(sort(tissue_2[[1]][,3],decreasing = TRUE))[0:100]
			size_of_intersection=length(intersect(names(top_100_tissue_1),names(top_100_tissue_2)))
			df[i,j]=size_of_intersection
			#print(size_of_intersection)
			
		}
	}
	df_null=df
	#cluster and create dendrogram
	hc_null=hclust(dist(df_null));dend=as.dendrogram(hc_null);y=labels(dend);y[]='black';y=cutree(hc_null,1);plot(hc_null)
	#put pdf in plot directory
	myPlot(hc_null,'dendrogram_random_100_MRs.pdf')
	write.table(df_null,file = "Overlap_100_MR_random.tsv",sep = "\t",quote = FALSE)
	dev.off()
	
}

GTEx_network_centrality_analysis <- function (network_name, somatic_mutation_file_name_1,somatic_mutation_file_name_2, gene_lengths) {
	#This function calculations a number of network centrality properties for the hub proteins, and, based on the somatic mutation data (somatic_mutation_file_name_1), the mutation frequency of the regulon. It then performs a regression analysis (with and 80%/20% train test split). This script also performs the same analysis for the second somatic mutation file (somatic_mutation_file_name_2). Finally, It creates a null hypothesis for somatic_mutation_file_name_1 by permuting the network 100 times and calcuting the regression correlation.
	#This function also reads in and calculates the total gene lengths of the regulons (using amino acid length) (file: gene_lengths). In this tab delimited file (without a header) the first element is the entrezgene name, and the second is the length of that gene.
	#In the end, I did not use this data for the function (for now at least)
	
	#The network read in is in the 6cols format that includes pvalues. The somatic mutation data is from the TCGA
	
	
	#needed for centralities
	library(igraph)
	#train/test split
	library(caret)
	#regression
	library(randomForest)
	#other utilities
	source("~/jb3401/scripts/Utils/rFunctions.R")
	#load("/Volumes/jb3401/archive-af_lab/GTEx/Aracne/Lung_vst/Lung_vst_6cols.rda",verbose=T)
	#network <- read.delim("~/jb3401/archive-af_lab/GTEx/Aracne/Lung_vst/Lung_vst_6cols.txt")
	network <- read.delim(network_name)
	#for pilot analysis, just use edges with pvalue=0
	network_old=network
	network=network[network$pvalue==0,]
	#load("~/jb3401/luad-somut.rda",verbose =T)
	load(somatic_mutation_file_name_1,verbose =T)
	network$neg_log_pvalue=1
	gene_lengths=read.delim("~/jb3401/scratch/GTEx_project/entrezgene-proteome872015_2_length.txt",header=F)
	rownames(gene_lengths)=gene_lengths$V1
	#create the graph, weighted by pvalue
	#rel=network[,c("Hub","Target","neg_log_pvalue")]
	rel=network[,c("Hub","Target")]
	rel_df=as.data.frame(rel)
	#colnames(rel_df)=c("Hub","Target","weight")
	colnames(rel_df)=c("Hub","Target")
	g= graph.empty(n=0, directed=FALSE)
	g=graph_from_data_frame(d=rel_df,directed = FALSE)
	#remove duplicate edges and loops etc.
	g=simplify(g,edge.attr.comb = "min")
	
	
	
	#remove genes that have no annotated corresponding length
	no_length=setdiff(V(g)$name,rownames(gene_lengths))
	g=delete_vertices(g,no_length)
	Hubs=unique(sorted(rel[,1]))
	Hubs=intersection(Hubs,V(g)$name)
	
	#Now, for each Hub, we will gather its regulons and count how many mutations are in it or its regulon based on a gene expression profile of somatic mutations
	#gather the regulons
	Hubs_df=data.frame(matrix(ncol=1,nrow=len(Hubs)))
	rownames(Hubs_df)=Hubs
	colnames(Hubs_df)=c("Regulon_num_muts")
	Hubs_df[,2]=NULL
	#colnames(Hubs_df)[2]="regulon_mutations"
	Hubs_df$regulon_gene_length=0
	Hubs_df$regulon_num_muts_normalized=0
	Hubs_df$regulon_size=0
	#Now, for each protein, get the number of its mutations normalized for the size of the regulon
	for (i in rownames(Hubs_df)){
		#get the regulon, including itself
		x=adjacent_vertices(g,i)
		x=names(x[[1]])
		x=append(x,i)
		x=unique(sorted(x))
		#get number of mutations
		num_mutations=sum(somut[intersection(rownames(somut),x),])/dim(somut)[2]
		#get length(in amino acids) of regulon
		genes_length=sum(gene_lengths[x,2])
		Hubs_df[i,1]=num_mutations
		Hubs_df[i,'regulon_gene_length']=genes_length
		Hubs_df[i,'regulon_size']=length(x)
		#Hubs_df[i,'regulon_num_muts_normalized']=num_mutations/genes_length
		Hubs_df[i,'regulon_num_muts_normalized']=num_mutations/length(x)
		print(i)
	}
	
	#Start calculating centralities
	#vertices_named_list=as.vector(flattenList(V(g)))
	#names(vertices_named_list)=rownames(Hubs_df)
	#all_degree=degree(g,V(g)$name)
	#all_alpha=alpha_centrality(g,V(g))
	all_articulation=articulation_points(g)
	all_authority=authority_score(g)
	all_betweenness=estimate_betweenness(g,V(g),cutoff=3)
	#betweenness=betweenness(g,V(g))
	all_closeness=estimate_closeness(g,V(g),cutoff=3)
	#closeness=closeness(g,V(g),cutoff=3)
	all_constraint=constraint(g,V(g))
	all_coreness=coreness(g)
	#all_diversity=diversity(g,V(g))
	all_eccentricity=eccentricity(g,V(g))
	all_eigen_centrality=eigen_centrality(g,V(g))
	all_hub_score=hub_score(g)
	all_knn=knn(g,V(g))
	all_local_scan_1=local_scan(g,k=1)
	#all_local_scan_2=local_scan(g,k=2)
	#all_min_seperators=min_separators(g)
	all_page_rank=page_rank(g)
	#all_power_centrality=power_centrality(g,V(g))
	#all_subgraph_centrality=subgraph_centrality(g)
	all_transitivity=transitivity(g,vids =V(g),type = "localundirected")
	
	the_Hubs=rownames(Hubs_df)
	#now collect these into a dataframe,keep in mind some failed
	#articulation points, articulation is binary, 1 for yes, 0 for no
	Hubs_df$articulation=0
	Hubs_df[intersection(all_articulation,rownames(Hubs_df)),'articulation']=1
	#diversity failed, alpha centrality failed min_seperators found none, power centrality failed
	Hubs_df$authority=all_authority[[1]][rownames(Hubs_df)]
	Hubs_df$betweenness=all_betweenness[the_Hubs]
	Hubs_df$closeness=all_closeness[the_Hubs]
	Hubs_df$constraint=all_constraint[the_Hubs]
	Hubs_df$coreness=all_coreness[the_Hubs]
	#Hubs_df$degree=all_degree[the_Hubs]
	Hubs_df$eccentricity=all_eccentricity[the_Hubs]
	Hubs_df$eigen=all_eigen_centrality[[1]][the_Hubs]
	Hubs_df$hub_score=all_hub_score[[1]][the_Hubs]
	Hubs_df$knn=all_knn[[1]][the_Hubs]
	Hubs_df$local_scan_1=all_local_scan_1[the_Hubs]
	#Hubs_df$local_scan_2=all_local_scan_2[the_Hubs]
	Hubs_df$page_rank=all_page_rank[[1]][the_Hubs]
	#Hubs_df$subgraph_centrality=all_subgraph_centrality[the_Hubs]
	names(all_transitivity)=V(g)$name
	Hubs_df$transitivity=all_transitivity[the_Hubs]
	
	#only get complete cases, so no NA's
	Hubs_df=Hubs_df[complete.cases(Hubs_df),]
	train_set=createDataPartition(Hubs_df$regulon_num_muts_normalized,p=.8,list = FALSE,times = 1)
	independant_variables=colnames(Hubs_df[5:dim(Hubs_df)[2]])
	response_variable="regulon_num_muts_normalized"
	#train random forest classifier
	result=randomForest(Hubs_df[train_set,independant_variables],Hubs_df[train_set,response_variable],ntree=100)
	predictions=predict(result,Hubs_df[-train_set,independant_variables])
	
	cor.test(Hubs_df[-train_set,response_variable],predictions)
	
	#generate a plot object with a linear fit
	df_to_plot=twoVectodataframe(Hubs_df[-train_set,response_variable],predictions,'actual','predictions')
	z=ggplot(df_to_plot,aes(x=actual, y=predictions)) + geom_point(size=2,shape=1) +geom_smooth(method=lm)
	myPlot(z,"Random_Forest_results_LUAD_LUNG.pdf")
	
	#train random forest classifier on train on randomized response variable
	result=randomForest(Hubs_df[train_set,independant_variables],Hubs_df[sample(train_set),response_variable],ntree=100)
	predictions=predict(result,Hubs_df[-train_set,independant_variables])
	
	cor.test(Hubs_df[-train_set,response_variable],predictions)
	df_to_plot=twoVectodataframe(Hubs_df[-train_set,response_variable],predictions,'actual','predictions')
	z=ggplot(df_to_plot,aes(x=actual, y=predictions)) + geom_point(size=2,shape=1)
	myPlot(z,"Random_Forest_results_LUAD_LUNG_dummy.pdf")
	
	#Save the Hubs feature data frame for lung:
	Hubs_df_real_lung=Hubs_df
	train_set_lung=train_set
	true_cor=cor.test(Hubs_df[-train_set,response_variable],predictions)$estimate
	
	
	#now let us try it with another cancer (skcm)
	
	load(somatic_mutation_file_name_2,verbose = T)
	#load("~/jb3401/gbm-somut.rda",verbose = T)
	#load("~/jb3401/coad-somut.rda",verbose = T)
	Hubs_df_dummy=Hubs_df
	Hubs_df_dummy$Regulon_num_muts=0
	Hubs_df_dummy$regulon_gene_length=0
	Hubs_df_dummy$regulon_num_muts_normalized=0
	#Only get hubs reported in the other tissue
	Hubs_df_dummy=Hubs_df_dummy[intersection(rownames(Hubs_df_dummy),rownames(somut)),]
	for (i in rownames(Hubs_df_dummy)){
		#get the regulon, including itself
		x=adjacent_vertices(g,i)
		x=names(x[[1]])
		x=append(x,i)
		x=unique(sorted(x))
		#get number of mutations
		num_mutations=sum(somut[intersection(rownames(somut),x),])/dim(somut)[2]
		#get length(in amino acids) of regulon
		genes_length=sum(gene_lengths[x,2])
		Hubs_df_dummy[i,1]=num_mutations
		Hubs_df_dummy[i,'regulon_gene_length']=genes_length
		#Hubs_df_dummy[i,'regulon_num_muts_normalized']=num_mutations/genes_length
		Hubs_df_dummy[i,'regulon_num_muts_normalized']=num_mutations/length(x)
		print(i)
		
		
	}
	Hubs_df_dummy=Hubs_df[complete.cases(Hubs_df_dummy),]
	#train_set=createDataPartition(Hubs_df_dummy$regulon_num_muts_normalized,p=.8,list = FALSE,times = 1)
	result=randomForest(Hubs_df_dummy[train_set,independant_variables],Hubs_df_dummy[train_set,response_variable],ntree=100)
	predictions=predict(result,Hubs_df_dummy[-train_set,independant_variables])
	
	cor.test(Hubs_df_dummy[-train_set,response_variable],predictions)
	#generate a plot object
	df_to_plot=twoVectodataframe(Hubs_df_dummy[-train_set,response_variable],predictions,'actual','predictions')
	z=ggplot(df_to_plot,aes(x=actual, y=predictions)) + geom_point(size=2,shape=1)+geom_smooth(method=lm)
	myPlot(z,"Random_Forest_results_skcm_LUNG.pdf")
	
	#Now we are going to do the same analysis but permuting the network 100 times
	train_set=train_set_lung
	null_hypothesis_cors=list()
	load("~/jb3401/luad-somut.rda",verbose =T)
	rel_df=as.data.frame(rel)
	colnames(rel_df)=c("Hub","Target")
	for (null_iteration in 1:100){
		print("Begin Iteration")
		#this scrambles the network
		rel_df$Target=sample(rel_df$Target)
		g= graph.empty(n=0, directed=FALSE)
		g=graph_from_data_frame(d=rel_df,directed = FALSE)
		g=simplify(g,edge.attr.comb = "min")
		#remove genes that have no annotated corresponding length
		no_length=setdiff(V(g)$name,rownames(gene_lengths))
		g=delete_vertices(g,no_length)
		Hubs=unique(sorted(rel[,1]))
		Hubs=intersection(Hubs,V(g)$name)
		
		Hubs_df=data.frame(matrix(ncol=1,nrow=len(Hubs)))
		rownames(Hubs_df)=Hubs
		colnames(Hubs_df)=c("Regulon_num_muts")
		Hubs_df[,2]=NULL
		#colnames(Hubs_df)[2]="regulon_mutations"
		Hubs_df$regulon_gene_length=0
		Hubs_df$regulon_num_muts_normalized=0
		Hubs_df$regulon_size=0
		#Now, for each protein, get the number of its mutations
		#we are going to correct for number of mutations for each regulon. Currently, the number of mutations will correlate with size of regulon (i.e. bigger regulon correlets with more random mutations). We will correct for this by looking up how long each member of the regulon is in amino acids, counting it up and dividing mutations by that number. E.g. a regulon that has a total of 1 mutation per 1000 amino acid length is the same as 10 mutations per 10,000
		for (i in rownames(Hubs_df)){
			#get the regulon, including itself
			x=adjacent_vertices(g,i)
			x=names(x[[1]])
			x=append(x,i)
			x=unique(sorted(x))
			#get number of mutations
			num_mutations=sum(somut[intersection(rownames(somut),x),])/dim(somut)[2]
			#get length(in amino acids) of regulon
			genes_length=sum(gene_lengths[x,2])
			Hubs_df[i,1]=num_mutations
			Hubs_df[i,'regulon_gene_length']=genes_length
			Hubs_df[i,'regulon_size']=length(x)
			#Hubs_df[i,'regulon_num_muts_normalized']=num_mutations/genes_length
			Hubs_df[i,'regulon_num_muts_normalized']=num_mutations/length(x)
			#print(i)
		}
		
		
		all_articulation=articulation_points(g)
		all_authority=authority_score(g)
		all_betweenness=estimate_betweenness(g,V(g),cutoff=3)
		all_closeness=estimate_closeness(g,V(g),cutoff=3)
		all_constraint=constraint(g,V(g))
		all_coreness=coreness(g)
		all_eccentricity=eccentricity(g,V(g))
		all_eigen_centrality=eigen_centrality(g,V(g))
		all_hub_score=hub_score(g)
		all_knn=knn(g,V(g))
		all_local_scan_1=local_scan(g,k=1)
		all_page_rank=page_rank(g)
		all_transitivity=transitivity(g,vids =V(g),type = "localundirected")
		the_Hubs=rownames(Hubs_df)
		Hubs_df$articulation=0
		Hubs_df[intersection(all_articulation,rownames(Hubs_df)),'articulation']=1
		Hubs_df$authority=all_authority[[1]][rownames(Hubs_df)]
		Hubs_df$betweenness=all_betweenness[the_Hubs]
		Hubs_df$closeness=all_closeness[the_Hubs]
		Hubs_df$constraint=all_constraint[the_Hubs]
		Hubs_df$coreness=all_coreness[the_Hubs]
		Hubs_df$eccentricity=all_eccentricity[the_Hubs]
		Hubs_df$eigen=all_eigen_centrality[[1]][the_Hubs]
		Hubs_df$hub_score=all_hub_score[[1]][the_Hubs]
		Hubs_df$knn=all_knn[[1]][the_Hubs]
		Hubs_df$local_scan_1=all_local_scan_1[the_Hubs]
		Hubs_df$page_rank=all_page_rank[[1]][the_Hubs]
		names(all_transitivity)=V(g)$name
		Hubs_df$transitivity=all_transitivity[the_Hubs]
		
		Hubs_df=Hubs_df[complete.cases(Hubs_df),]
		#train_set=createDataPartition(Hubs_df$regulon_num_muts_normalized,p=.8,list = FALSE,times = 1)
		independant_variables=colnames(Hubs_df[5:dim(Hubs_df)[2]])
		response_variable="regulon_num_muts_normalized"
		#train random forest classifier
		result=randomForest(Hubs_df[train_set,independant_variables],Hubs_df[train_set,response_variable],ntree=100)
		predictions=predict(result,Hubs_df[-train_set,independant_variables])
		
		the_cor=cor.test(Hubs_df[-train_set,response_variable],predictions)$estimate
		null_hypothesis_cors=append(null_hypothesis_cors,the_cor)
		print(null_iteration)
		print(the_cor)
		
	}
	true_cor_df=as.data.frame(true_cor)
	true_cor_df$y_val=.1
	null_hypothesis_cors_2=flattenList(null_hypothesis_cors)
	#null_hypothesis_cors_2=null_hypothesis_cors_2**2
	null_hypothesis_cors_2_df=as.data.frame(null_hypothesis_cors_2)
	z=ggplot(data=null_hypothesis_cors_2_df,aes(x=null_hypothesis_cors_2))+geom_histogram()+geom_point(data=true_cor_df,aes(x=true_cor,y=y_val),color='blue')
	z=z+xlim(-.2,.5)
	z=z+ggtitle("Histogram of correlations")
	myPlot(z,"correlation_GTEX_network_null_hypotheses.pdf")
	the_mean=mean(null_hypothesis_cors_2_df$null_hypothesis_cors_2)
	the_sd=sd(null_hypothesis_cors_2_df$null_hypothesis_cors_2)
	the_pvalue=pnorm(true_cor,sd=the_sd,mean=the_mean,lower.tail = FALSE)
}


#This function analyses aracne tf/cotf networks to  see how well they overlap with a ppi network. 
# We pass to it a file name with the arance network for co-transcription factors (in the 4 column 
# format , protein 1, protein 2, MI and pvalue (tab delimted)), the aracne network for the 
# transcription factors (also in 4 column format) and a PPI netwokr (protein 1, protein 2, strength 
# of interaction).  Note that for this analysis to work, in the protein/protein interaction 
# network, the entries must be duplicate (i.e. A,B,strength and B,A,strength must BOTH be present)
# This function calculates all pairwise TF/coTF fisher exact test scores for overlap of the regulons 
# (pvalue=0) and then plots a histogram of the distribution, it overlaps this histogram with the 
# those pairs of TF/coTF just found in the PPI network as well.
# RETURNS:
# This function returns 2 objects in a list: 1. a dataframe where the rows are the names of the 
# aracne_tf_network (i.e. all tfs), columns are all cotfs, and the values are the FET scores 2. 
# a dataframe that shows which pairs are in the PPI network and their mutual information
AracnePrePPI_TF_coTF <-function (aracne_cotf_network,aracne_tf_network,ppi_network){
	library(igraph)
	library(caret)
	library(randomForest)
	library(pryr)
	source("~/jb3401/scripts/Utils/rFunctions.R")
	setwd("/Users/joshuabroyde/jb3401/scripts/R_sessions")
	#load("/Volumes/jb3401/archive-af_lab/GTEx/Aracne/Lung_vst/Lung_vst_6cols.rda",verbose=T)
	aracne_cotf<- read.delim(aracne_cotf_network,header=F)
	aracne_tf <- read.delim(aracne_tf_network,header=F)
	#now read in the PrePPI network
	preppi_df=read.delim(ppi_network,header=F)
	
	#only keep pvalues at zero
	aracne_cotf=aracne_cotf[aracne_cotf$V4==0,]
	aracne_tf=aracne_tf[aracne_tf$V4==0,]
	aracne_tf$V1=as.character(aracne_tf$V1)
	aracne_tf$V2=as.character(aracne_tf$V2)
	aracne_cotf$V2=as.character(aracne_cotf$V2)
	aracne_cotf$V1=as.character(aracne_cotf$V1)
	#do pairwise fets of shared regulons
	aracne_network=network[network$pvalue==0,]
	count=0
	all_tfs=unique(sort(aracne_tf$V1))
	all_tfs=as.character(all_tfs)
	all_cotfs=unique(sort(aracne_cotf$V1))
	all_cotfs=as.character(all_cotfs)
	all_gene_counts=len(unique(sorted(flattenList(c(aracne_tf$V2,aracne_cotf$V2)))))
	df_1=namedrowscolsDataframe(all_tfs,all_cotfs)
	#do all pairwise fisher exact tests, and store in a datafreame call df_1
	for (i in all_tfs){
		the_tf_regulon=aracne_tf[aracne_tf$V1==i,'V2']
		for (j in all_cotfs){
			the_co_tf_regulon=aracne_cotf[aracne_cotf$V1==j,'V2']
			the_intersection=intersection(the_tf_regulon,the_co_tf_regulon)
			the_pvalue=fisher.test(matrix(c(len(the_intersection),len(the_tf_regulon),len(the_co_tf_regulon),all_gene_counts),ncol = 2))$p.value
			df_1[i,j]=the_pvalue
			print(count) #print count to see progress
			count=count+1
		}
	}
	all_tfs_cotfs=unique(sorted(flattenList(c(all_tfs,all_cotfs))))
	
	preppi_df$V1=as.character(preppi_df$V1)
	preppi_df$V2=as.character(preppi_df$V2)
	preppi_df=preppi_df[preppi_df$V1 %in% all_tfs,]
	preppi_df=preppi_df[preppi_df$V2 %in% all_cotfs,]
	preppi_df$V4=1
	for (i in 1:nrow(preppi_df)){
		p1=preppi_df[i,1]
		p2=preppi_df[i,2]
		the_fet_score=df_1[p1,p2]
		preppi_df[i,4]=the_fet_score
		print(p1)
	} 
	pdf("TF-CoTF_HIST.pdf")
	hist(unlist(df_1),freq = FALSE,breaks = c(0,.05,.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.75,.8,.85,.9,1))
	hist(preppi_df$V4,freq = FALSE,add=T,col=rgb(0,0,1,0.5),breaks = c(0,.05,.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.75,.8,.85,.9,1))
	dev.off()
	all_values=unlist(df_1)
	#This is a fisher exact test for the .05 threshold
	preppi_below_.05=nrow(preppi_df[preppi_df$V4<=.05,])
	all_preppi=nrow(preppi_df)
	all_below_.05=len(all_values[all_values<=.05])
	all_all=len(all_values[all_values<=1.1])
	fisher_results=fisher.test(matrix(c(preppi_below_.05,all_preppi,all_below_.05,all_all),ncol=2))
	return(c(df_1,preppi_df))
	#Save the file if you want to go back to it.
	#save.image(file=".PrePPI_TF_coTF_Aracne_analysis.RData")
}

#This function creates two histograms (overlaid on each other) to compare Fisher Exact (FET)
#scores for proteins predicted in a ppi network, compared to all pairwise FET scores
#ppi_df: is a 3-column dataframe (the names of the columnsdo not matter), where the first column is
#the name of the first protein, the second column is the name of the third protein, and the third 
#column is  the strength of the interaction (this score is not used, but must still be provided)
#FET_df is a m x n matrix, where rows and columns are, respectivley transcription factors and 
#cofactor names, and the cells values represent the FET scores for the overalap of the regulon
#all_tfs and all_cotfs are lists of transcription factors and cofactors respectivley. 

#This function will first find all the fet scores for protein pairs in ppi_df (that are also tfs and cotfs)
#and plot a histogram of the FET scores, and create a similair plot for all FET scores in FET_df
#The blue bars are the histogram for ppi_df, clear is histogram for FET_df
#This function also returns the ppi_df with a new column, corresponding the the values in the FET matrix for those proteins.


network_arance_histogram <- function(ppi_df,FET_df,all_tfs,all_cotfs,outfile)
{
	preppi_df=ppi_df
	df_1=FET_df
	colnames(preppi_df)=c("V1","V2","V3")
	rownames(preppi_df)=NULL
	preppi_df$V1=as.character(preppi_df$V1)
	preppi_df$V2=as.character(preppi_df$V2)
	preppi_df=preppi_df[preppi_df$V1 %in% all_tfs,]
	preppi_df=preppi_df[preppi_df$V2 %in% all_cotfs,]
	preppi_df[,4]=1
	#return(preppi_df)
	for (i in 1:nrow(preppi_df)){
		p1=preppi_df[i,1]
		p2=preppi_df[i,2]
		the_fet_score=df_1[p1,p2]
		preppi_df[i,4]=the_fet_score
		print(p1)
	} 
	pdf(outfile)
	hist(unlist(df_1),freq = FALSE,breaks = c(0,.05,.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.75,.8,.85,.9,.95,1))
	hist(preppi_df$V4,freq = FALSE,add=T,col=rgb(0,0,1,0.5),breaks = c(0,.05,.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.75,.8,.85,.9,.95,1))
	dev.off()
	return(preppi_df)
	
}


#This function converts the regulon object to a dataframe
#e.g. Skin_not_sun_df=RegulonToDataframe(SkinNotSun)
RegulonToDataframe <- function(regulon_object){
	temp=regulon_object[[2]][1][[1]]
	for (i in 2:len(regulon_object[[2]])){
		temp=rbind(temp,regulon_object[[2]][i][[1]])
		print(i)
	}
	rownames(temp)=NULL
	temp_df=as.data.frame(temp)
	return(temp_df)
}


#This function takes in a dataframe that representents an aracne network and creates a pairwise m x m matrix 
#that is the FET score between each hub.
#Note that this function will take several hours to run for an interactome including 6000 hubs.
#Input: a dataframe of the form:
#Hub	Target
#19	360
#19	421
#19	604
#19	627
#19	858
#19	904

#output: a matrix of the form:
#	row.names	10001	10002	10009
#1	10001	9.16115e-85	1.000000e+00	1.00000e+00
#2	10002	1.00000e+00	1.604302e-58	1.00000e+00
#3	10009	1.00000e+00	1.000000e+00	2.91621e-222

pairwiseFET <- function (aracne_network_dataframe){
	
	library(igraph)
	library(reshape2)
	library(data.table)
	source("~/jb3401/scripts/Utils/rFunctions.R")
	#df_1=read.delim("~/jb3401/archive-af_lab/GTEx/Aracne2/Colon-Sigmoid_clean_vst/Colon-Sigmoid_clean_vst_signal.tsv",header=F)
	#df_1=read.delim(input_file,header=F)
	#only keep pvalue of less than .001 for the edges
	#aracne_network=df_1[df_1$V4<=.001,]
	
	the_hubs=as.character(unique(sorted(aracne_network$V1)))
	total_potential_regulons=as.character(unique(sorted(aracne_network$V2)))
	total_potential_regulon_size=len(total_potential_regulons)
	
	#convert the aracne network to a graph object in igraph
	rel_df=as.data.frame(aracne_network[,c("V1","V2")])
	colnames(rel_df)=c("Hub","Target")
	rel_df$probability=NULL
	g= graph.empty(n=0, directed=FALSE)
	g=graph_from_data_frame(d=rel_df,directed = FALSE)
	#remove duplicate edges if there are any.
	g=simplify(g)
	
	
	the_final_matrix=namedrowscolsDataframe(the_hubs,the_hubs)
	#the final matrix pairwise will contain the FET overlap and Jaccard index for all pairs of hub proteins
	the_final_matrix_pairwise=melt(as.matrix(the_final_matrix))
	colnames(the_final_matrix_pairwise)=c("hub_1","hub_2","hub_1_size")
	the_final_matrix_pairwise$hub_2_size=NA
	the_final_matrix_pairwise$the_intersection=NA
	the_final_matrix_pairwise$FET_pvalue=NA
	the_final_matrix_pairwise$the_total=total_potential_regulon_size
	the_final_matrix_pairwise$hub_1=as.character(the_final_matrix_pairwise$hub_1)
	the_final_matrix_pairwise$hub_2=as.character(the_final_matrix_pairwise$hub_2)
	
	#get the numbers needed for FET test
	the_final_matrix_pairwise$hub_1_size=degree(g,the_final_matrix_pairwise$hub_1)
	the_final_matrix_pairwise$hub_2_size=degree(g,the_final_matrix_pairwise$hub_2)
	
	to_extract=cbind(the_final_matrix_pairwise$hub_1,the_final_matrix_pairwise$hub_2)
	#the cocitation matrix generated by "cocitation" is the the same as the overlap of all the hubs.
	#It is an m x m matrix of the hubs and the overlaps.
	cocitation_matrix=cocitation(g)
	
	the_values=cocitation_matrix[to_extract]
	the_final_matrix_pairwise$the_intersection=the_values
	cocitation_matrix_hubs_df=melt(cocitation_matrix_hubs)
	
	#correct cocitation so that same ones have intersection of the degree
	#i.e. (cocitation of A and A is the size of A regulon, )
	blub=the_final_matrix_pairwise$hub_1==the_final_matrix_pairwise$hub_2
	the_final_matrix_pairwise$the_intersection[blub]=the_final_matrix_pairwise$hub_1_size[blub]
	
	the_final_matrix_pairwise_next=the_final_matrix_pairwise[,c('hub_1_size','hub_2_size','the_intersection','the_total')]
	
	the_final_matrix_pairwise_next$hub_1_size=the_final_matrix_pairwise_next$hub_1_size-the_final_matrix_pairwise_next$the_intersection
	the_final_matrix_pairwise_next$hub_2_size=the_final_matrix_pairwise_next$hub_2_size-the_final_matrix_pairwise_next$the_intersection
	the_final_matrix_pairwise_next$the_total=the_final_matrix_pairwise_next$the_total-the_final_matrix_pairwise_next$hub_1_size-the_final_matrix_pairwise_next$hub_2_size-the_final_matrix_pairwise_next$the_intersection
	
	#Get all FET Scores, this is the part of the function that takes a long time
	all_pvalues=apply(the_final_matrix_pairwise_next,1, 
			function(i)
			{
				x=unlist(i)
				fisher.test(matrix(c(x[3],x[1],x[2],x[4]),ncol=2))$p.value
			}
	)
	
	the_final_matrix_pairwise$FET_pvalue=all_pvalues
	
	
	return(the_final_matrix_pairwise)
}


#given an m x m adjaceny matrix that is the interactome, random walk with restart is performed, and neighborhoods
#are found for each member of the interactome.
#The strength of the interaction is assumed to be p-values, and thus the negative log is taken before doing random walk.
#Note that duplicate neighborhoods may be found, and, that neighborhoods may not be symetric (ie. A may be a neighbor of B)
#but B is not a neighbor of A).
#For actually performing random walk with restart, this function uses the personalized pagerank algorithm, which is equivalent.
#The random walk with restart result for each gene is z scaled before returning.
#If you want to exclude low confidence interations, then replace those with 0 before passing to the function
#usage example: e.g the_neighborhoods=randomWalkerRegulon(TF_CoTF_all_pairwise_FET)
#This function takes about 2 hours for a 2200 x 22000 matrix.
randomWalkerRegulon <- function (regulon_mat){
	source("~/jb3401/scripts/Utils/rFunctions.R")
	#diag(TF_CoTF_all_pairwise_FET)=1
	library(igraph)
	library(reshape2)
	log_transformed_matrix=-log10(regulon_mat)
	#replace diagonal with 0
	diag(log_transformed_matrix)=0
	#replace Inf with second highest value:
	log_transformed_matrix=replaceInfwithSecondhighest(log_transformed_matrix)
	
	#remove interactions that are zero
	log_transformed_matrix_df=melt(as.matrix(log_transformed_matrix))
	log_transformed_matrix_df_over_zero=log_transformed_matrix_df[log_transformed_matrix_df$value>0,]
	#create a graph out of the interaction matrix, where the weight is negative log pvalue
	colnames(log_transformed_matrix_df_over_zero)=c("protein1","protein2","weight")
	rel_df=as.data.frame(log_transformed_matrix_df_over_zero)
	colnames(rel_df)=c("Hub","Target","weight")
	g= graph.empty(n=0, directed=FALSE)
	g=graph_from_data_frame(d=rel_df,directed = FALSE)
	#do the random walk for each one
	pairwise_random_walk=personalizedPagerank_pairwise(g,V(g)$name,V(g)$name)
	pairwise_random_walk_df=as.data.frame(pairwise_random_walk)
	#scale the results by converting to z scores:
	scaled_pairwise_random_walk=t(scale(t(pairwise_random_walk)))
	scaled_pairwise_random_walk_df=as.data.frame(scaled_pairwise_random_walk)
	return(scaled_pairwise_random_walk_df)
	
}


#this function performs personalized pagerank/random walk with restart for each member of list_1 and returns a matrix of m x n, 
#where m is the the length of list 1 and m is the length of list 2. For each row of the returned matrix, the scores
#are the normalized random walk probability (they all sum to 1). For example, if the rowname is "geneA",
#then the scores of row gene A are the probably if getting to the other genes when doing pagerank/rwr when starting at "geneA".
#"g" is a graph object of the underlying network generated by igraph  
personalizedPagerank_pairwise <- function (g,list_1,list_2){
	pp_matrix=namedrowscolsDataframe(list_1,list_2)
	input_vect=as.vector(rep(0,len(V(g)$name)))
	names(input_vect)=V(g)$name
	input_vect_old=input_vect
	count=1
	for (i in list_1){
		input_vect=input_vect_old
		input_vect[i]=1
		z=page_rank(g,personalized = input_vect)[1]
		z=z[[1]]
		results=z[list_2]
		pp_matrix[i,names(results)]=results
		print(count)
		count=count+1
	}
	return(pp_matrix)
}


# function to replace Inf with the second highest value
replaceInfwithSecondhighest <- function(the_dataframe) {
	new_value=sort(unique(unlist(the_dataframe)),decreasing = T)[2]
	the_dataframe2=do.call(data.frame,lapply(the_dataframe, function(x) replace(x, is.infinite(x),new_value)))
	the_dataframe2=as.data.frame(the_dataframe2)
	rownames(the_dataframe2)=rownames(the_dataframe)
	colnames(the_dataframe2)=colnames(the_dataframe)
	return(the_dataframe2)
}

aracneConservationStatistics <- function(file_null,file_real,output_plot,output_rda){
	#This function reads in the results of conservation between 10 null/scrambled aracne networks and
	#the actual aracne network conservation, and generates plots about these results. As well as n RDA file
	#In the results, each row represents  how frequently an interaction if the aracne networks (with DPI) is seen across
	#ALL the aracne networks (without DPI)
	library(data.table)
	#df_null=read.delim("/ifs/scratch/c2b2/ac_lab/jb3401/temp_work_dir/scrambled_networks/null_networks_corrected.txt",header=T)
	#df_null=fread("/ifs/home/c2b2/ac_lab/jb3401/jb3401/scratch/GTEx_project/compare_to_noDPI/TCGA/all_networks_null.txt")
	df_null=fread(file_null)
	#df_real=read.delim("/ifs/scratch/c2b2/ac_lab/jb3401/temp_work_dir/num_interactions_real.txt",header=F)
	#df_real=fread("/ifs/home/c2b2/ac_lab/jb3401/jb3401/scratch/GTEx_project/compare_to_noDPI/TCGA/true_results.txt")
	df_real=fread(file_real)
	df_real$V1=NULL 
	df_null$ID=NULL
	normalized_null=rowSums(df_null)/sum(df_null)
	normalized_real=rowSums(df_real)/sum(df_real)
	
	normalized_null_2=rowSums(df_null)
	normalized_real_2=rowSums(df_real)
	pdf(output_plot)
	plot(2:65,log10(normalized_real_2[2:65]),col='pink',ylim=c(0,9))
	points(2:65,log10(normalized_null_2[2:65]))
	dev.off()
	probability_ratio=normalized_real/normalized_null
	probability_ratio_2=na.omit(probability_ratio)
	probability_ratio_2=replaceInfwithSecondhighestVector(probability_ratio_2)
	plot(1:len(probability_ratio_2),log10(probability_ratio_2),pch=21,bg="lightgreen")
	dev.off()
	
	probabilty_ratio_data_tissue_conservation=as.data.frame(cbind(normalized_null,normalized_real,normalized_real/normalized_null))
	colnames(probabilty_ratio_data_tissue_conservation)=c("random_probability","real_probability","probability_ratio")
	save(probabilty_ratio_data_tissue_conservation,file=output_rda)
}

aracneGTExTCGAregressionAnalysis <-function(){
	#This function performs an extrapolation analysis and outputs r objects (and plots) of how conservation level corresponds to probability of an interaction
	#being in a true aracne network vs a scrambled network. 
	#Note that it calculates also prior and posterior probabilities, but this is not used actually.
	prior=34417467/125423879 #total number of detected interactions divided by total number of theoretically possible interactions (6109*20531)
	pdf("GTEx_TCGA_regression_analysis.pdf")
	load("~/jb3401/scratch/RDA_files/probabilty_ratio_data_tissue_conservation_GTEx_only.rda",verbose = T)
	
	posterior=prior*(probabilty_ratio_data_tissue_conservation$real_probability)/(prior*probabilty_ratio_data_tissue_conservation$real_probability+(1-prior)*probabilty_ratio_data_tissue_conservation$random_probability)
	
	probabilty_ratio_data_tissue_conservation$posterior=posterior
	probabilty_ratio_data_tissue_conservation$conditional_probability=probabilty_ratio_data_tissue_conservation$real_probability/(probabilty_ratio_data_tissue_conservation$real_probability+probabilty_ratio_data_tissue_conservation$random_probability)
	probabilty_ratio_data_tissue_conservation$num_tissues=1:(nrow(probabilty_ratio_data_tissue_conservation))
	probabilty_ratio_data_tissue_conservation=probabilty_ratio_data_tissue_conservation[complete.cases(probabilty_ratio_data_tissue_conservation),]
	x=probabilty_ratio_data_tissue_conservation$num_tissues
	y=probabilty_ratio_data_tissue_conservation$posterior
	#y=probabilty_ratio_data_tissue_conservation$probability_ratio
	
	#fit to logistic function with asymptote at 1
	fit <- nls(y~1/(1+ exp(-b * (x-c))), start=list(b=.5,c=1))
	new_y=predict(fit,newdata=x)
	#for x's where we could not empirically estimate y, impute using this regression model.
	new_values=y[y<1]
	new_values=append(new_values,new_y[!y<1])
	probabilty_ratio_data_tissue_conservation$posterior_regression=new_values
	
	#do the same with conditional probability
	x=probabilty_ratio_data_tissue_conservation$num_tissues
	y=probabilty_ratio_data_tissue_conservation$conditional_probability
	
	#y=probabilty_ratio_data_tissue_conservation$probability_ratio
	
	#fit to logistic function with asymptote at 1
	fit <- nls(y~1/(1+ exp(-b * (x-c))), start=list(b=.5,c=1))
	new_y=predict(fit,newdata=x)
	#for x's where we could not empirically estimate y, impute using this regression model.
	new_values=y[y<1]
	new_values=append(new_values,new_y[!y<1])
	probabilty_ratio_data_tissue_conservation$conditional_probability=new_values
	probabilty_ratio_data_tissue_conservation$prior=prior
	probabilty_ratio_data_tissue_conservation$random_probability_regression=probabilty_ratio_data_tissue_conservation$real_probability/probabilty_ratio_data_tissue_conservation$conditional_probability-probabilty_ratio_data_tissue_conservation$real_probability
	probabilty_ratio_data_tissue_conservation$random_probability_regression=probabilty_ratio_data_tissue_conservation$random_probability_regression/sum(probabilty_ratio_data_tissue_conservation$random_probability_regression)
	plot(probabilty_ratio_data_tissue_conservation$num_tissues,probabilty_ratio_data_tissue_conservation$random_probability)
	lines(probabilty_ratio_data_tissue_conservation$num_tissues,probabilty_ratio_data_tissue_conservation$random_probability_regression,col='blue')
	plot(probabilty_ratio_data_tissue_conservation$num_tissues,probabilty_ratio_data_tissue_conservation$random_probability,ylim=c(0,.0001))
	
	lines(probabilty_ratio_data_tissue_conservation$num_tissues,probabilty_ratio_data_tissue_conservation$random_probability_regression,ylim=c(0,.0001),col='blue')
	plot(probabilty_ratio_data_tissue_conservation$num_tissues,probabilty_ratio_data_tissue_conservation$conditional_probability,pch=20,col='green')
	
	probabilty_ratio_data_tissue_conservation=probabilty_ratio_data_tissue_conservation[,c('num_tissues','random_probability_regression','real_probability','prior','posterior','conditional_probability')]
	save(probabilty_ratio_data_tissue_conservation,file="~/jb3401/scratch/RDA_files/probabilty_ratio_data_tissue_conservation_GTEx_only_2.rda")
	
	load("~/jb3401/scratch/RDA_files/probabilty_ratio_data_tissue_conservation_TCGA_only.rda",verbose = T)
	prior=20163141/125423879
	posterior=prior*(probabilty_ratio_data_tissue_conservation$real_probability)/(prior*probabilty_ratio_data_tissue_conservation$real_probability+(1-prior)*probabilty_ratio_data_tissue_conservation$random_probability)
	probabilty_ratio_data_tissue_conservation$posterior=posterior
	
	
	probabilty_ratio_data_tissue_conservation$posterior=posterior
	probabilty_ratio_data_tissue_conservation$num_tissues=1:(nrow(probabilty_ratio_data_tissue_conservation))
	probabilty_ratio_data_tissue_conservation=probabilty_ratio_data_tissue_conservation[complete.cases(probabilty_ratio_data_tissue_conservation),]
	x=probabilty_ratio_data_tissue_conservation$num_tissues
	y=probabilty_ratio_data_tissue_conservation$posterior
	probabilty_ratio_data_tissue_conservation$conditional_probability=probabilty_ratio_data_tissue_conservation$real_probability/(probabilty_ratio_data_tissue_conservation$real_probability+probabilty_ratio_data_tissue_conservation$random_probability)
	#y=probabilty_ratio_data_tissue_conservation$probability_ratio
	
	#fit to logistic function with asymptote at 1
	fit <- nls(y~1/(1+ exp(-b * (x-c))), start=list(b=.5,c=1))
	new_y=predict(fit,newdata=x)
	#for x's where we could not empirically estimate y, impute using this regression model.
	new_values=y[y<1]
	new_values=append(new_values,new_y[!y<1])
	probabilty_ratio_data_tissue_conservation$posterior_regression=new_values
	
	#do the same with conditional probability
	x=probabilty_ratio_data_tissue_conservation$num_tissues
	y=probabilty_ratio_data_tissue_conservation$conditional_probability
	#y=probabilty_ratio_data_tissue_conservation$probability_ratio
	
	#fit to logistic function with asymptote at 1
	fit <- nls(y~1/(1+ exp(-b * (x-c))), start=list(b=.5,c=1))
	new_y=predict(fit,newdata=x)
	#for x's where we could not empirically estimate y, impute using this regression model.
	new_values=y[y<1]
	new_values=append(new_values,new_y[!y<1])
	probabilty_ratio_data_tissue_conservation$conditional_probability=new_values
	probabilty_ratio_data_tissue_conservation$prior=prior
	probabilty_ratio_data_tissue_conservation$random_probability_regression=probabilty_ratio_data_tissue_conservation$real_probability/probabilty_ratio_data_tissue_conservation$conditional_probability-probabilty_ratio_data_tissue_conservation$real_probability
	
	probabilty_ratio_data_tissue_conservation$random_probability_regression=probabilty_ratio_data_tissue_conservation$random_probability_regression/sum(probabilty_ratio_data_tissue_conservation$random_probability_regression)
	
	plot(probabilty_ratio_data_tissue_conservation$num_tissues,probabilty_ratio_data_tissue_conservation$random_probability)
	lines(probabilty_ratio_data_tissue_conservation$num_tissues,probabilty_ratio_data_tissue_conservation$random_probability_regression,col='blue')
	plot(probabilty_ratio_data_tissue_conservation$num_tissues,probabilty_ratio_data_tissue_conservation$random_probability,ylim=c(0,.0001))
	lines(probabilty_ratio_data_tissue_conservation$num_tissues,probabilty_ratio_data_tissue_conservation$random_probability_regression,ylim=c(0,.0001),col='blue')
	probabilty_ratio_data_tissue_conservation=probabilty_ratio_data_tissue_conservation[,c('num_tissues','random_probability_regression','real_probability','prior','posterior','conditional_probability')]
	
	plot(probabilty_ratio_data_tissue_conservation$num_tissues,probabilty_ratio_data_tissue_conservation$conditional_probability,pch=20,col='green')
	save(probabilty_ratio_data_tissue_conservation,file="~/jb3401/scratch/RDA_files/probabilty_ratio_data_tissue_conservation_TCGA_only_2.rda")
	dev.off()
	
	
}

updateAracneLikelihoods <-function(){
	#This function reads in conservation probabilities and ARACNe networks, and updates the "likelihood" column for the Aracne networks
	#using Bayes Therom. THE UPDATE is done IN PLACE and written to disk, so be careful!!! 
	library(data.table)
	the_interactome=fread("/ifs/home/c2b2/ac_lab/jb3401/jb3401/scratch/Aracne_Networks/GTEx_TCGA_networks/tcga_luad_6cols.txt")
	load("/ifs/home/c2b2/ac_lab/jb3401/jb3401/scratch/RDA_files/probabilty_ratio_data_tissue_conservation_TCGA_only_2.rda")
	probabilty_ratio_data_tissue_conservation_TCGA=probabilty_ratio_data_tissue_conservation
	load("/ifs/home/c2b2/ac_lab/jb3401/jb3401/scratch/RDA_files/probabilty_ratio_data_tissue_conservation_GTEx_only_2.rda")
	probabilty_ratio_data_tissue_conservation_GTEx=probabilty_ratio_data_tissue_conservation
	
	#Test correlation between between GTEx and TCGA
	the_intersection=all_networks_together_TCGA[.(all_networks_together_GTEx$V1)]
	the_intersection_2=the_intersection[the_intersection$V2>=0,]
	the_intersection_3=all_networks_together_GTEx[.(the_intersection_2$V1)]
	cor.test(the_intersection_3$V2,the_intersection_2$V2)
	gtex_tcga_r_squared=cor.test(the_intersection_3$V2,the_intersection_2$V2)$estimate^2
	#png("Correlation_GTEx_TCGA.png")
	#plot(the_intersection_3$V2,the_intersection_2$V2)
	#dev.off()
	
	all_networks_together_TCGA=fread("/ifs/home/c2b2/ac_lab/jb3401/jb3401/scratch/GTEx_project/compare_to_noDPI/TCGA/temp_TCGA_networks_dir/output_file.txt",verbose=T)
	setkeyv(all_networks_together_TCGA,c("V1","V2"))
	
	all_networks_together_GTEx=fread("/ifs/home/c2b2/ac_lab/jb3401/jb3401/scratch/GTEx_project/compare_to_noDPI/GTEx/temp_GTEx_networks_dir/output_file.txt",verbose=T)
	setkeyv(all_networks_together_GTEx,c("V1","V2"))
	
	#now we are going to update the interactomes. THIS UPDATE IS DONE IN PLACE AND WRITTEN TO DISK!
	
	all_networks=read.delim("/ifs/scratch/c2b2/ac_lab/jb3401/Aracne_Networks/conserved_network_list.txt",header=F)
	
	for (the_network in all_networks$V1){
		the_network=as.character(the_network)
		the_interactome=fread(the_network)
		temp=paste(the_interactome$Hub,the_interactome$Target,sep="_")
		temp_GTEx=all_networks_together_GTEx[.(temp)]
		names(temp_reg_values)=probabilty_ratio_data_tissue_conservation_GTEx[,c("num_tissues")]
		old_probs=the_interactome$likelihood
		new_probs=updateProbabilities(old_probs,temp_GTEx$V2,probabilty_ratio_data_tissue_conservation_GTEx)
		
		
		
		#now perform this analysis with TCGA data
		
		temp=paste(the_interactome$Hub,the_interactome$Target,sep="_")
		temp_TCGA=all_networks_together_TCGA[.(temp)]
		names(temp_reg_values)=probabilty_ratio_data_tissue_conservation_TCGA[,c("num_tissues")]
		old_probs=new_probs
		new_probs_2=updateProbabilities(old_probs,temp_TCGA$V2,probabilty_ratio_data_tissue_conservation_TCGA)
		the_interactome$likelihood=new_probs_2
		write.table(the_interactome,row.names=F,file=the_network,col.names=TRUE,quote=FALSE,sep="\t")
		print(the_network)
	}
}

#update probablities using bayes theorem
updateProbabilities <- function(prior_probabilities,tissue_counts,prob_data_frame){
	#use bayes theorem
	A=prior_probabilities
	prior=prob_data_frame$prior[1]
	B_given_A=prob_data_frame[tissue_counts,'real_probability']
	B_give_not_A=prob_data_frame[tissue_counts,'random_probability_regression']
	#denominator=(B_given_A*prior)+(B_give_not_A)*(1-prior)
	denominator=(B_given_A*A)+(B_give_not_A)*(1-A)
	#if data is missing, replace with old probablility
	new_probability=(A*B_given_A)/denominator
	new_probability[is.na(new_probability)]=prior_probabilities[is.na(new_probability)]
	return(new_probability)
	
}

compareOriginalNewLikelihoodAracne <-function(){
	library(data.table)
	the_interactome=fread("/ifs/home/c2b2/ac_lab/jb3401/jb3401/scratch/Aracne_Networks/GTEx_TCGA_networks/tcga_luad_6cols.txt")
	the_interactome_conservation=fread("/ifs/home/c2b2/ac_lab/jb3401/jb3401/scratch/Aracne_Networks/GTEx_TCGA_tissue_conserved_networks/tcga_luad_6cols.txt")
	pdf("ARACNE_tissue_Conservation_correlation.pdf")
	plot(the_interactome$likelihood,the_interactome_conservation$likelihood,pch=20)
	dev.off()
	
	#Now read in raw values and plot
}

checkIndependenceofGTExTCGA <- function (){
	library(data.table)
	all_networks_together_TCGA=fread("/ifs/home/c2b2/ac_lab/jb3401/jb3401/scratch/GTEx_project/compare_to_noDPI/TCGA/temp_TCGA_networks_dir/output_file.txt",verbose=T)
	setkeyv(all_networks_together_TCGA,c("V1","V2"))
	
	all_networks_together_GTEx=fread("/ifs/home/c2b2/ac_lab/jb3401/jb3401/scratch/GTEx_project/compare_to_noDPI/GTEx/temp_GTEx_networks_dir/output_file.txt",verbose=T)
	setkeyv(all_networks_together_GTEx,c("V1","V2"))
	both=intersection(all_networks_together_GTEx$V1,all_networks_together_TCGA$V1)
	print("Interactions common to both")
	print(len(both))
	correlation_squared=cor.test(all_networks_together_GTEx[.(both)]$V2,all_networks_together_TCGA[.(both)]$V2)$estimate^2
	print(correlation_squared)
}

#This function compares how well the old and new versions (i.e. tissue conserved versions) of ARACNE capture the STRING database:
compareAracneToString <- function(){
	library(ROCR)
	source("~/jb3401/scripts/Utils/rFunctions.R")
	#load("/ifs/home/c2b2/ac_lab/jb3401/jb3401/scratch/GTEx_project/Compare_to_STRING/LUNG_original_aracne.R",verbose=T)
	#load("/ifs/home/c2b2/ac_lab/jb3401/jb3401/scratch/GTEx_project/Compare_to_STRING/LUNG_modified_aracne.R",verbose=T)
	string_data=read.delim("/ifs/scratch/c2b2/ac_lab/jb3401/GTEx_project/Compare_to_STRING/STRING_500_literature_2.txt",header=F)
	modified_aracne=read.delim("/ifs/home/c2b2/ac_lab/jb3401/jb3401/scratch/GTEx_project/Compare_to_STRING/Lung_modified_aracne_2.txt",header=F)
	rownames(modified_aracne)=modified_aracne$V1
	modified_aracne$V1=NULL
	old_aracne=read.delim("/ifs/home/c2b2/ac_lab/jb3401/jb3401/scratch/GTEx_project/Compare_to_STRING/Lung_original_aracne_2.txt",header=F)
	rownames(old_aracne)=old_aracne$V1
	old_aracne$V1=NULL
	
	both=intersection(string_data$V1,rownames(old_aracne))
	old_aracne[both,'V3']=1
	old_aracne$V3=0
	modified_aracne$V3=0
	modified_aracne[both,'V3']=1
	
	pred=prediction(modified_aracne$V2,modified_aracne$V3);
	perf_mod=performance(pred,measure = "tpr", x.measure = "fpr");
	pred=prediction(old_aracne$V2,old_aracne$V3);
	perf_orig=performance(pred,measure = "tpr", x.measure = "fpr");
	pdf("Old_aracne_new_aracne_compared_to_String.pdf")
	plot(perf_mod,col='red')
	par(new=T)
	plot(perf_orig,col='grey')
	abline(0,1)
	dev.off()
}

# *********************************************************************************
# Return the VST-normalized expression measurements for a query gene.
#
# This method assumes that the contents of the following have have been already loaded
# into the R environment:
#	/ifs/archive/shares/af_lab/GTEx/RDA_Files/all_64_expmat.rda
#
# ARGUMENTS:
# * gene:		The query gene.
# * mat_names:	Names of gene expression matrix variables. This must be a vector of
#		strings containing a subset of values from variable matNames. Or, it can
#		be NULL in which case it is assumed that all values in matNames are to be
#		used, i.e., mat_names = matNames
#
# RETURNS
# A list "res" with length(mat_names) entries, one each for the values in 
# mat_names. The entry:
#		res[[mat_names[i]]]
# is a vector containing the normalized expression measurements of "gene" in the 
# tissue identified by mat_names[i]. The list is named, specifically:
#		names(res) = mat_names.
# *********************************************************************************
getGeneExpression <- function(gene, mat_names = NULL){
	if(is.na(gene) || is.null(gene))
		return(NA)
	
	gene <- as.entrezId(gene)
	if (is.null(mat_names))
		mat_names = matNames
	res = lapply(mat_names, function(mat){
				mat = get(mat)
				return(as.numeric(mat[as.character(gene), ]))
			})
	names(res) = mat_names
	return(res)
}

#'****************************************************************************************
#' Generates tables listing the support and the maximum MI for every interaction.
#'
#' @title	computeSupportAndMaxMI
#' 
#' @param net	A subset of "varNames", i.e., a character string containing interactome 
#' 				variable names.
#' 
#' @return 		Returns a list with two members, each member being a N x M numeric matrix
#' 				where N is the total number of genes found in all interactomes specified 
#' 				in "nets" and M is the total number of hubs found in all interactomes in
#' 				"nets". In each matrix, rownames() and colnames() are the entrez Ids of 
#' 				the N and M genes, respectively. The entry [gene, hub] in the first 
#' 				matrix is the support of the interaction (hub, gene), i.e., the number of 
#' 				interactomes in "nets" where the interaction is present. In the second
#' 				matrix, the entry [gene, hub] is the maximum MI value for the interaction 
#' 				(hub, gene) across all interactomes in "nets".
#' 				NOTE: this is a time consuming method.
#'****************************************************************************************
computeSupportAndMaxMI <- function(nets = varNames){
  all_genes = getInteractomeGenes(nets, count = FALSE, hubs_only = FALSE, common = FALSE)
  all_hubs = getInteractomeGenes(nets, count = FALSE, hubs_only = TRUE, common = FALSE)
  sup = matrix(0, length(all_genes), length(all_hubs))
  maxmi = matrix(0, length(all_genes), length(all_hubs))
  rownames(sup) = rownames(maxmi) = all_genes
  colnames(sup) = colnames(maxmi) = all_hubs
  for (net in nets){
    writeLines(paste("Now processing --->", net))
    net = get(net)
    for (hub in names(net[[2]])){
      reg = net[[2]][[hub]]
      for (i in 1:nrow(reg)){
        if (reg[i, "Target"] > 0){
          sup[rownames(reg)[i], hub] = sup[rownames(reg)[i], hub] + 1
          maxmi[rownames(reg)[i], hub] = max(maxmi[rownames(reg)[i], hub], reg[i, "MI"])
        }
      }
    }
  }
  
  return(list(sup, maxmi))
}
