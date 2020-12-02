#!/bin/bash

#SBATCH -p exacloud             # partition (queue) (exacloud, long_jobs)
#SBATCH -N 1                    # number of nodes
#SBATCH -n 8                    # number of cores
#SBATCH --mem 64000             # memory pool for all cores
#SBATCH -t 0-1:00               # time (D-HH:MM)

##### Checking Dependency Success #####
dependid=${SLURM_JOB_DEPENDENCY##*:}

deptasks=(`sacct -j $dependid --format State | grep -v [S/-] | tr -d ' ' `)

check_success ${deptasks[@]}

# Start script
echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID


# set paths 
#wrkdir=/home/groups/hoolock/u1/jg/probe_val/ECP4wd
#args_file=/home/groups/hoolock/u1/jg/probe_val/ECP4_cov_args
# epic_seq renamed to probe_bed
#epic_seq=/home/groups/hoolock/u1/jg/probe_val/hg38-truseq-methyl-capture-epic-manifest-file.bed 

# set variable "ARGS" to be the output of the correct line from args_file (matches SLURM_ARRAY_TASK_ID)
ARGS=`head -$SLURM_ARRAY_TASK_ID $args_file | tail -1`

# set variables based on the output of the line ARGS
#IFS=: read ORIGBAM  INFILE<<< $ARGS

## JG 10/29/20: Adding Single/Paired End Sequencing Formatting. 
# set variables based on the output of the line ARGS
IFS=: read INFILE INFILE2 SEQ <<< $ARGS

## JG 10/29/20: modified script to take in the user supplied sample name. 
if [ ! -f $wrkdir/cov.$INFILE2.bed ]; then
    bedtools coverage -sorted \
    -g $wrkdir/${expname}_genome.txt \
    -a $wrkdir/srt-${probe_name} \
    -b $sourcedir/reads/srt.$INFILE2.bam > $wrkdir/cov.$INFILE2.bed
fi
