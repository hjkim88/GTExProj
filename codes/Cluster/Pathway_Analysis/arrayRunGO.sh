#!/bin/bash

# The code that will run on each CPU. Loads the R module on the cluster, changes to the parent
# directory where the results will be stored (the current directory is intended only for 
# saving the stdout and stderr files), and invokes the actual R script.

module load R/3.3.1
cd ..
Rscript ~/CPTAC/scripts/regulonGOenrichment_array.r $NET $ONT $CHUNK $SGE_TASK_ID
