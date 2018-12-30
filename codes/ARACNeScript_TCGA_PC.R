###
#   File name : ARACNeScript_TCGA_PC.R
#   Author    : Hyunjin Kim
#   Date      : Dec 30, 2018
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Make a shell script for Aracne runs with the TCGA data of our own (Will be run on Hyunjin's PC)
#
#   Instruction
#               1. Source("ARACNeScript_TCGA_PC.R")
#               2. Run the function "makeScript()" - specify the input file (Aracne-ready data) directory and output file path
#               3. The ARACNe script will be generated in the output file path
#
#   Example
#               > source("The_directory_of_ARACNeScript_TCGA_PC.R/ARACNeScript_TCGA_PC.R")
#               > makeScript(inputPath="./results/aracne_ready/TCGA_our_own/",
#                            outputFilePath="./results/scripts/TCGA_our_own/tcga_our_own_aracne_script.sh",
#                            dos2unixPath="C:/cygwin64/bin/dos2unix.exe")
###

makeScript <- function(inputPath="./results/aracne_ready/TCGA_our_own/",
                       outputFilePath="./results/scripts/TCGA_our_own/tcga_our_own_aracne_script.sh",
                       dos2unixPath="C:/cygwin64/bin/dos2unix.exe") {
  
  ### collect files from the inputPath
  f <- list.files(inputPath)
  
  
  ### only keep .dat files
  f <- f[which(endsWith(f, ".dat"))]
  
  
  ### set global variables in the script
  script <- "aracne=/cygdrive/c/Research/CUMC/ARACNe/aracne\n"
  script <- paste0(script, "tf=C:/Research/CUMC/ARACNe/hubGenes/tf.txt\n")
  script <- paste0(script, "cotf=C:/Research/CUMC/ARACNe/hubGenes/cotf.txt\n")
  script <- paste0(script, "signaling=C:/Research/CUMC/ARACNe/hubGenes/signaling.txt\n")
  script <- paste0(script, "pv=0.00000001\n")
  script <- paste0(script, "threads=4\n")
  script <- paste0(script, "bootstrap=100\n")
  script <- paste0(script, "\n")
  
  
  ### iteratively write Aracne run code for each TCGA tissue
  for(i in 1:length(f)) {
    ### set Aracne ready data path variable
    script <- paste0(script, "input=C:/Research/CUMC/GTExProj/results/aracne_ready/TCGA_our_own/",
                     f[i], "\n")
    
    ### make an output directory for the tissue
    dir <- strsplit(f[i], split = ".", fixed = TRUE)[[1]][1]
    script <- paste0(script, "mkdir /cygdrive/c/Research/CUMC/GTExProj/results/Aracne/TCGA_our_own/",
                     dir, "\n")
    
    ### set Aracne run output directory path variable
    script <- paste0(script, "output=C:/Research/CUMC/GTExProj/results/Aracne/TCGA_our_own/",
                     dir, "/\n")
    
    ### Aracne run with TF hubs
    script <- paste0(script, "/usr/bin/time -v ${aracne} -e ${input} -o ${output}tf_run -t",
                     " ${tf} --pvalue ${pv}  --threads ${threads} --bootstrap ${bootstrap}",
                     "> ${output}tf_log.txt 2>&1", "\n")
    
    ###Aracne run with co-TF hubs
    script <- paste0(script, "/usr/bin/time -v ${aracne} -e ${input} -o ${output}cotf_run -t",
                     " ${cotf} --pvalue ${pv}  --threads ${threads} --bootstrap ${bootstrap}",
                     "> ${output}cotf_log.txt 2>&1", "\n")
    
    ### Aracne run with signaling hubs
    script <- paste0(script, "/usr/bin/time -v ${aracne} -e ${input} -o ${output}signal_run -t",
                     " ${signaling} --pvalue ${pv}  --threads ${threads} --bootstrap ${bootstrap}",
                     "> ${output}signal_log.txt 2>&1", "\n")
    
    ### add line for each tissue run
    script <- paste0(script, "\n")
  }
  
  
  ### save the master script
  write.table(script, file = outputFilePath, sep="\n", quote=FALSE, row.names=FALSE, col.names=FALSE)
  
  
  ### dos2unix to the result
  system(paste(dos2unixPath, outputFilePath))
  
}


