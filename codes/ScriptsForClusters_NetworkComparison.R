###
#   File name : ScriptsForClusters_NetworkComparison.R
#   Author    : Hyunjin Kim
#   Date      : Dec 31, 2019
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : This will create R scripts and shell scripts for cluster and
#               a master shell script that will submitting all the shell scripts
#               for calculating similarities of all the Aracne network pairs
#
#   Instruction
#               1. Source("ScriptsForClusters_NetworkComparison.R")
#               2. Run the function "make_scripts" - specify the necessary input and output paths
#               3. The result file will be generated in the output path
#
#   Example
#               > source("The_directory_of_ScriptsForClusters_NetworkComparison.R/ScriptsForClusters_NetworkComparison.R")
#               > make_scripts(igraphPath="./data/RDA_Files/All_62_Aracne_igraphs.rda",
#                              sampleRPath="./results/network_comparison/cluster_scripts/sample_randomwalk.R",
#                              sampleSHPath="./results/network_comparison/cluster_scripts/sample_randomwalk.sh",
#                              outputDir="./results/network_comparison/cluster_scripts/",
#                              dos2unixPath="C:/cygwin64/bin/dos2unix.exe")
###

make_scripts <- function(igraphPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/GTEx/RDA_Files/All_62_Aracne_igraphs.rda",
                         sampleRPath="//isilon.c2b2.columbia.edu/ifs/scratch/c2b2/af_lab/hk2990/network_comparison/sample_randomwalk.R",
                         sampleSHPath="//isilon.c2b2.columbia.edu/ifs/scratch/c2b2/af_lab/hk2990/network_comparison/sample_randomwalk.sh",
                         outputDir="//isilon.c2b2.columbia.edu/ifs/scratch/c2b2/af_lab/hk2990/network_comparison/",
                         dos2unixPath="//isilon.c2b2.columbia.edu/ifs/home/c2b2/af_lab/hk2990/dos2unix.exe") {
  
  ### load example shell script
  sample_r <- read.table(file=sampleRPath, sep='\n', header=FALSE, check.names = FALSE, stringsAsFactors = FALSE, comment.char = "", quote = "")
  sample_sh <- read.table(file=sampleSHPath, sep='\n', header=FALSE, check.names = FALSE, stringsAsFactors = FALSE, comment.char = "", quote = "")
  
  ### load igraphs
  load(igraphPath)
  
  ### master script's first two lines
  masterScript <- "#!/bin/bash\nqacct -o hk2990\n"
  
  ### make R and Sh scripts for each tissue
  for(tissue in names(igs)) {
    
    ### let's make a R script for a specific tissue!
    scriptSample <- sample_r
    
    ### change the tissue name in the sample script
    scriptSample[4,1] <- paste0("given_tissue <- \"", tissue, "\"")
    
    ### save the R script
    write.table(scriptSample, file=paste0(outputDir, tissue, ".R"), sep="\n", quote=FALSE, row.names=FALSE, col.names=FALSE)
    
    ### dos2unix to the result
    system(paste(dos2unixPath, paste0(outputDir, tissue, ".R")))
    
    ### for Shell script
    scriptSample <- sample_sh
    
    ### change some things in the sample script
    scriptSample[2,1] <- paste0("#$ -l mem=12G,time=80:00:00 -S /bin/bash -cwd -j y -N hk_", tissue)
    scriptSample[4,1] <- paste0("R_code=/ifs/scratch/c2b2/af_lab/hk2990/network_comparison/", tissue, ".R")
    
    ### save the Shell script
    write.table(scriptSample, file=paste0(outputDir, tissue, ".sh"), sep="\n", quote=FALSE, row.names=FALSE, col.names=FALSE)
    
    ### dos2unix to the result
    system(paste(dos2unixPath, paste0(outputDir, tissue, ".sh")))
    
    ### append the master script
    masterScript <- paste0(masterScript, "qsub /ifs/scratch/c2b2/af_lab/hk2990/network_comparison/", tissue, ".sh\n")
    
  }
  
  ### save the master script
  write.table(masterScript, file=paste0(outputDir, "MASTER_SCRIPT.sh"), sep="\n", quote=FALSE, row.names=FALSE, col.names=FALSE)
  
  ### dos2unix to the result
  system(paste(dos2unixPath, paste0(outputDir, "MASTER_SCRIPT.sh")))
  
}
