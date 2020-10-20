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

# GO!
if [ ! -f $wrkdir/${expname}_cov.bed ]; then
    Rscript --vanilla $script_dir/cov_merge.R $args_file $wrkdir $expname
fi

