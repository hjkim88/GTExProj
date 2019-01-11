###
#   File name : ARACNeScripts.R
#   Author    : Hyunjin Kim
#   Date      : Dec 7, 2017
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Make Aracne scripts for HPC cluster submission
#
#   Instruction
#               1. Source("ARACNeScripts.R")
#               2. Run the function "makeScripts()" - specify the input file (Aracne-ready data) directory and output directory
#               3. The ARACNe scripts will be generated in the output directory
#
#   Example
#               > source("The_directory_of_ARACNeScripts.R/ARACNeScripts.R")
#               > makeScripts(examplePath="../../ARACNe/GTEx/Aracne_example.sh",
#                             inputPath="./results/aracne_ready/",
#                             outputPath="./results/scripts/",
#                             dos2unixPath="C:/cygwin64/bin/dos2unix.exe")
###

makeScripts <- function(examplePath="../ARACNe/GTEx/Aracne_example.sh",
                        inputPath="./results/aracne_ready/",
                        outputPath="./results/scripts/",
                        dos2unixPath="C:/cygwin64/bin/dos2unix.exe") {
  
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
    
    ### Aracne requires at least 100 samples
    if(ncol(cnt) >= 100) {
      
      ### let's make a script for a specific file!
      scriptSample <- example
      
      scriptSample[2,1] <- paste0(scriptSample[2,1], "_", i)
      scriptSample[3,1] <- paste0("input=./data/", f[i])
      scriptSample[4,1] <- paste0("output=./results/", substr(f[i], 1, nchar(f[i])-4), "/")
      
      ### save the result
      write.table(scriptSample, file=paste0(outputPath, substr(f[i], 1, nchar(f[i])-4), ".sh"), sep="\n", quote=FALSE, row.names=FALSE, col.names=FALSE)
      
      ### dos2unix to the result
      system(paste(dos2unixPath, paste0(outputPath, substr(f[i], 1, nchar(f[i])-4), ".sh")))
      
      ### create a result directory for a specific file
      dir.create(paste0(outputPath, substr(f[i], 1, nchar(f[i])-4)))
      
      ### append the master script
      masterScript <- paste0(masterScript, "qsub ./", substr(f[i], 1, nchar(f[i])-4), ".sh\n")
      
    }
  }
  
  ### save the master script
  write.table(masterScript, file=paste0(outputPath, "MASTER_SCRIPT.sh"), sep="\n", quote=FALSE, row.names=FALSE, col.names=FALSE)

  ### dos2unix to the result
  system(paste(dos2unixPath, paste0(outputPath, "MASTER_SCRIPT.sh")))
  
}


