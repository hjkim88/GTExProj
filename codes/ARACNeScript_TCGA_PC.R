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
  
  
  
  
}


