This folder contains code for performing cluster-enabled GO term enrichment analysis 
on the regulons of all networks. The code assumes that root execution has been setup
on the "/ifs/scratch" partition with three subfolders (see file "arrayRunMaster.sh"
for details):
- "data": contains .rda files, one for each interactome to be processed.
- "runs": empty folder where the output from the cluster jobs will be stored.
- "scripts":  contains the script files needed for the run.

The file "arrayRunMaster.sh" is the master script, this is what should be run from
the commans line of a cluster node to submit the whole cluster computation. The 
CHUNK variable specifies how many regulons each cluster job will process. Assuming
that processing a single regulon takes about a minute, a value of 30 is reasonable.
Notice that this is setup as a so-called "array" cluster job. This is specified by
the "-t 1-200:1" argument which requires 200 CPUs for each interactome. The number 
200 is chosen based on the total number of hub genes, which is between 5000-6000. 
Assuming a CHUNK setting of 30, this means that there are enough available CPUS to
finish processing an entire interactome in about 30 minutes or so.
*************************************************************************************
---------------------------       !!!ATTENTION!!!          --------------------------
*****ATTENTION********: make sure to set the memory and time parameters as appropriate,
to ensure enough resources are available for a job to complete. This will require a 
bit of experimentation, i.e., run a couple of jobs with selected CHUNK setting, to 
test what are the resource requirements.
*************************************************************************************

Other relevant files are:
- "arrayRunGO.sh": this is the script file that will run on each of the 200 CPUs.
- "regulonGOenrichment_array.r": this is the R code invoked by arrayRunGO.sh. It 
	is the one that does all the work.
- "cleanupAll.r": this should be run after all the cluster jobs are finished, to 
	aggregate the individual CHUNKs generated by each array job into on final 
	object per interactome.