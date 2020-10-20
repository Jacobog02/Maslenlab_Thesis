#!/bin/bash

#SBATCH -p exacloud             # partition (queue) (exacloud, long_jobs)
#SBATCH -N 1                    # number of nodes
#SBATCH -n 8                    # number of cores
#SBATCH --mem 24000             # memory pool for all cores
#SBATCH -t 0-00:30               # time (D-HH:MM)

echo "SLURM_JOBID: " $SLURM_JOBID
#echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
#echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID

## This script is dependent on the coverage script being run sucessfully and will process the output of each file into one merged file named $expname_cov.bed. 
##### Checking Dependency Success #####
dependid=${SLURM_JOB_DEPENDENCY##*:}
deptasks=(`sacct -j $dependid --format State | grep -v [S/-] | tr -d ' ' `)

check_success ${deptasks[@]}


## Old! 4/15/20
# GO!
#Rscript -e 'Sys.setenv(RSTUDIO_PANDOC="/usr/lib/rstudio-server/bin/pandoc"); rmarkdown::render("probe_visuals.Rmd","html_document") ' $expname $wrkdir  
#Rscript -e 'rmarkdown::render("probe_visuals.Rmd","html_document") ' $expname $wrkdir $notation 

#mv probe_visuals.html $wrkdir/${expname}_visuals.html


## NEW
# Updated to have a wrapper R script to catch command line arguements
Rscript $script_dir/vis_wrapper.R $expname $wrkdir $script_dir


# Remap output to be comparable.
# if [ -n "$notation" ]; then
# #echo Im doing something > beep.boop
# # Convert the bedfile to UCSC if needed
# cat $wrkdir/${expname}_failed_probes.bed | $notation $mapping flip > $wrkdir/${expname}_temp.bed
# # Overwrite the old file bc this one is standard
# mv $wrkdir/${expname}_temp.bed $wrkdir/${expname}_failed_probes.bed
# # Clear the tmp file
# rm $wrkdir/${expname}_temp.bed
# 
# fi
