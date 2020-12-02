#!/bin/bash

#SBATCH -p exacloud                # partition (queue) (exacloud, long_jobs)
#SBATCH -N 2                       # number of nodes
#SBATCH -n 8                       # number of cores
#SBATCH --mem 20000              # 36 GB memory pool for all cores
#SBATCH -t 00-01:00                 # time (D-HH:MM)

##### Checking Dependency Success #####
dependid=${SLURM_JOB_DEPENDENCY##*:}

deptasks=(`sacct -j $dependid --format State | grep -v [S/-] | tr -d ' ' `)

check_success ${deptasks[@]}

echo "SLURM_JOBID: " $SLURM_JOBID



## set outdir path
#$wrkdir outdir=/home/groups/hoolock/u1/jg/replicate/qc/output


#	TASKS:
#echo $INFILE
# Do something to each input file:
#$exafastqc -t 2 -o $wrkdir/qc $sourcedir/$INFILE*.fastq.gz 
## Check expected output: If it doesnt exist execute otherwise skip. 
$exa_mutliqc -d $wrkdir/ -n ${expname}_report_multiqc


