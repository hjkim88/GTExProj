###
#   File name : CPNetworksOnly.R
#   Author    : Hyunjin Kim
#   Date      : Dec 14, 2017
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Network files are stored in server. I need to copy them to another server directory.
#
#   Instruction
#               1. Source("CPNetworksOnly.R")
#               2. Run the function "cp()" - specify source path and output directory
#               3. The bash copy script will be generated in the output directory
#
#   Example
#               > source("The_directory_of_CPNetworksOnly.R/CPNetworksOnly.R")
#               > cp(sourcePath="/ifs/scratch/c2b2/af_lab/hk2990/ARACNe_GTEx/results/",
#                    destinationPath="/ifs/archive/shares/af_lab/GTEx/Aracne/",
#                    outputDir="./results/scripts/CPNetworks.sh",
#                    dos2unixPath="C:/cygwin64/bin/dos2unix.exe")
###

cp <- function(sourcePath="/ifs/scratch/c2b2/af_lab/hk2990/ARACNe_GTEx/results/",
               destinationPath="/ifs/archive/shares/af_lab/GTEx/Aracne/",
               outputDir="./results/scripts/CPNetworks.sh",
               dos2unixPath="C:/cygwin64/bin/dos2unix.exe") {
  
  ### script's first line
  script <- "#!/bin/bash\n"
  
  ### move to the working directory
  script <- paste0(script, "cd ", sourcePath, "\n")
  
  ### get the directories
  script <- paste0(script, "DIRS='ls -d *'\n")
  
  ### mkdir and cp
  script <- paste0(script, "for i in $($DIRS); do\n")
  script <- paste0(script, "mkdir ", destinationPath, "${i}\n")
  script <- paste0(script, "cp ${i}/signal_run/finalNetwork_4col.tsv ", destinationPath, "${i}/${i}_signal.tsv\n")
  script <- paste0(script, "cp ${i}/cotf_run/finalNetwork_4col.tsv ", destinationPath, "${i}/${i}_cotf.tsv\n")
  script <- paste0(script, "cp ${i}/tf_run/finalNetwork_4col.tsv ", destinationPath, "${i}/${i}_tf.tsv\n")
  script <- paste0(script, "done\n")
  
  ### save the script
  write.table(script, file=outputDir, sep="\n", quote=FALSE, row.names=FALSE, col.names=FALSE)
  
  ### dos2unix to the result
  system(paste(dos2unixPath, outputDir))
  
}



