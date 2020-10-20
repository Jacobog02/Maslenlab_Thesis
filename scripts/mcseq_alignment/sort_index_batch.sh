#!/bin/bash

#SBATCH -p exacloud                # partition (queue) (exacloud, long_jobs)
#SBATCH -N 2                      # number of nodes
#SBATCH -n 24                       # number of cores
#SBATCH --mem 64000                # 36 GB memory pool for all cores
#SBATCH -t 0-02:00                 # time (D-HH:MM)

##### Checking Dependency Success #####
dependid=${SLURM_JOB_DEPENDENCY##*:}

deptasks=(`sacct -j $dependid --format State | grep -v [S/-] | tr -d ' ' `)

check_success ${deptasks[@]}

echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID


# Set Variables But when generalized the variables are in enviroment. 
#sourcedir=/home/groups/hoolock/u1/bd/Projects/ECP4/bismark/output3_min20max75
#wrkdir=/home/groups/hoolock/u1/jg/probe_val/ECP4wd
#args_file=/home/groups/hoolock/u1/jg/probe_val/raw_ECP4_args.txt.fin

# set variable "ARGS" to be the output of the correct line from args_file (matches SLURM_ARRAY_TASK_ID)
ARGS=`head -$SLURM_ARRAY_TASK_ID $args_file | tail -1`
## JG 10/11/20: Adding Single/Paired End Sequencing Formatting. 
# set variables based on the output of the line ARGS
#IFS=: read INFILE INFILE2 <<< $ARGS
IFS=: read INFILE INFILE2 SEQ <<< $ARGS

## JG 10/12/20: Single/Paired End mod
## I will do this by changing the suffix for file. 
## Taken from bismark_post_align_batch.sh
if [ $SEQ == "SE" ]; then
  
  suffix=".bam"
  #param="-s"
  
elif [ $SEQ == "PE" ]; then 
  
  suffix="_pe.bam"
  #param="-p"
  
else
   ## This field is not set correctly
   echo $SEQ "Field is not PE or SE"
   exit 0 ## leave
fi

#	TASKS:
# Do something to each input file:
catch=`ls $wrkdir/reads | grep ".bam$" | grep "^srt.${INFILE2}[_|\.]"`


if [ ! -f $wrkdir/reads/$catch ]; then 

    echo "Sorting and Indexing" ${INFILE2}

    $exa_sam sort -o $wrkdir/reads/srt.${INFILE2}.bam -@ 12 -m 8G $wrkdir/align/${INFILE2}${suffix} #_pe.bam

    $exa_sam index -@ 8 $wrkdir/reads/srt.${INFILE2}.bam
    
    else
    echo "Previous output detected now skipping"
fi

