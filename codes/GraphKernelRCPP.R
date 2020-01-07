###
#   File name : GraphKernelRCPP.R
#   Author    : Hyunjin Kim
#   Date      : Dec 26, 2019
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : We would like to run graph kernel functions in the "graphkernels" package,
#               but we want to modify some of the CPP code of it. So, we are testing
#               the new function with the modified CPP code.
#
#   Instruction
#               1. Source("GraphKernelRCPP.R")
#               2. Run the function "geometricRandomWalk" - specify the necessary input and output paths
#               3. The result file will be generated in the output path
#
#   Example
#               > source("The_directory_of_GraphKernelRCPP.R/GraphKernelRCPP.R")
#               > geometricRandomWalk(cppPath="./codes/GeometricRandomWalk.cpp")
###

geometricRandomWalk <- function(cppPath="//isilon.c2b2.columbia.edu/ifs/scratch/c2b2/af_lab/hk2990/network_comparison/GeometricRandomWalk.cpp") {
  
  ### load libraries
  if(!require(Rcpp, quietly = TRUE)) {
    install.packages("Rcpp")
    require(Rcpp, quietly = TRUE)
  }
  if(!require(graphkernels, quietly = TRUE)) {
    install.packages("graphkernels")
    require(graphkernels, quietly = TRUE)
  }
  
  ### load the CPP function
  ### here the function name is already designated: CalculateKernelCpp()
  sourceCpp(file = cppPath)
  
  ### the wrapper function
  CalculateGeometricRandomWalkKernel <- function(G, par) {
    graph.info.list <- vector("list", length(G))
    for (i in 1:length(G))
      graph.info.list[[i]] <- GetGraphInfo(G[[i]])

    return(CalculateKernelCpp(graph.info.list, par, 8))
  }
  
  ### sample run
  data(mutag)
  writeLines(paste(CalculateGeometricRandomWalkKernel(list(mutag[[1]], mutag[[2]]), 0.1)))
  
}
