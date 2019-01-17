###
#   File name : PrepBeforeReComputeMI.R
#   Author    : Hyunjin Kim
#   Date      : Mar 12, 2017
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Make a bash script that does the following:
#
# ****************************************************************************************
# 1. Make a scripts that re-run the aracne command on the same dataset you used for
#    creating the original network using the following parameter settings:
# "aracne -e <exp_file> -o <output_dir> -t <hub_list> --pvalue 1.0 --nobootstrap --nodpi"
# 2. Create the output directory <output_dir>.
# 3. In the directory create a file named "miThreshold_p1E0_samples<N>.txt"
#    where <N> is the number of samples in the expression matrix that ARACNE will run on.
#    E.g., if the expression matrix has 200 samples the file name should be named
#    "miThreshold_p1E0_samples200.txt". In that file add the following single line: 0.0
# 4. You should combine the three hub gene lists (tf, co_tf, signaling) into on list
#    and add to the file <hub_list>
# ****************************************************************************************
#
#   Instruction
#               1. Source("PrepBeforeReComputeMI.R")
#               2. Run the function "prepScript()" - specify source path and output directory
#               3. The bash script will be generated in the output directory
#
#   Example
#               > source("The_directory_of_PrepBeforeReComputeMI.R/PrepBeforeReComputeMI.R")
#               > prepScript(examplePath="../ARACNe/GTEx/Aracne_recomp_mi_example.sh",
#                            inputPath="./results/aracne_ready/GTEx2/vst_clean/",
#                            outputPath="./results/scripts/GTEx2/MI/",
#                            hubPath="../ARACNe/hubGenes/",
#                            dos2unixPath="C:/cygwin64/bin/dos2unix.exe")
###

prepScript <- function(examplePath="../ARACNe/GTEx/Aracne_recomp_mi_example.sh",
                       inputPath="./results/aracne_ready/GTEx2/vst_clean/",
                       outputPath="./results/scripts/GTEx2/MI/",
                       hubPath="../ARACNe/hubGenes/",
                       dos2unixPath="C:/cygwin64/bin/dos2unix.exe") {
  
  ### STEP 1
  ### load example shell script
  example <- read.table(file=examplePath, sep='\n', header=FALSE, check.names = FALSE, stringsAsFactors = FALSE, comment.char = "", quote = "")
  
  ### collect files from the cntPath
  f <- list.files(inputPath)
  f <- f[which(endsWith(f, ".dat"))]
  
  ### master script's first two lines
  masterScript <- "#!/bin/bash\nqacct -o hk2990\n"
  
  ### make scripts for all the datasets
  for(i in 1: length(f)) {
    
    ### load dataset
    cnt <- read.table(paste0(inputPath, f[i]), sep="\t", row.names = 1, header = TRUE, check.names = FALSE)
    
    ### let's make a script for a specific file!
    scriptSample <- example
    
    scriptSample[2,1] <- paste0(scriptSample[2,1], "_", i)
    scriptSample[3,1] <- paste0("input=../data/", f[i])
    scriptSample[4,1] <- paste0("output=./results/", substr(f[i], 1, nchar(f[i])-4), "_mi/")
    
    ### save the result
    write.table(scriptSample, file=paste0(outputPath, substr(f[i], 1, nchar(f[i])-4), "_mi.sh"), sep="\n", quote=FALSE, row.names=FALSE, col.names=FALSE)
    
    ### dos2unix to the result
    system(paste(dos2unixPath, paste0(outputPath, substr(f[i], 1, nchar(f[i])-4), "_mi.sh")))
    
    ### STEP 2
    ### create a result directory for a specific file
    dir.create(paste0(outputPath, substr(f[i], 1, nchar(f[i])-4), "_mi"))
    
    ### STEP 3
    ### create a fixed mi threshold file for Aracne run
    temp <- "-99"
    write.table(temp, file=paste0(outputPath, substr(f[i], 1, nchar(f[i])-4), "_mi/", "miThreshold_p1E0_samples", ncol(cnt), ".txt"), sep="", quote=FALSE, row.names=FALSE, col.names=FALSE)
    
    ### dos2unix to the result
    system(paste(dos2unixPath, paste0(outputPath, substr(f[i], 1, nchar(f[i])-4), "_mi/", "miThreshold_p1E0_samples", ncol(cnt), ".txt")))
    
    ### append the master script
    masterScript <- paste0(masterScript, "qsub ./", substr(f[i], 1, nchar(f[i])-4), "_mi.sh\n")
    
  }
  
  ### save the master script
  write.table(masterScript, file=paste0(outputPath, "MASTER_SCRIPT_MI.sh"), sep="\n", quote=FALSE, row.names=FALSE, col.names=FALSE)
  
  ### dos2unix to the result
  system(paste(dos2unixPath, paste0(outputPath, "MASTER_SCRIPT_MI.sh")))
  
  ### STEP 4
  ### combining tf + cotf + signaling list
  cotfs <- read.table(paste0(hubPath, "cotf.txt"), sep='\n', header=FALSE, check.names = FALSE, stringsAsFactors = FALSE, comment.char = "", quote = "")
  tfs <- read.table(paste0(hubPath, "tf.txt"), sep='\n', header=FALSE, check.names = FALSE, stringsAsFactors = FALSE, comment.char = "", quote = "")
  signals <- read.table(paste0(hubPath, "signaling.txt"), sep='\n', header=FALSE, check.names = FALSE, stringsAsFactors = FALSE, comment.char = "", quote = "")

  hubs <- rbind(cotfs, tfs, signals)

  ### save the total hubs
  write.table(hubs, file=paste0(hubPath, "hubs.txt"), sep="\n", quote=FALSE, row.names=FALSE, col.names=FALSE)

  ### dos2unix to the result
  system(paste(dos2unixPath, paste0(hubPath, "hubs.txt")))
  
}
