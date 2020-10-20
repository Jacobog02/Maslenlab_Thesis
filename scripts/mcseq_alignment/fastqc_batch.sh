#!/bin/bash

#SBATCH -p exacloud                # partition (queue) (exacloud, long_jobs)
#SBATCH -N 2                       # number of nodes
#SBATCH -n 16                       # number of cores
#SBATCH --mem 2400              # 36 GB memory pool for all cores
#SBATCH -t 00-02:00                 # time (D-HH:MM)

echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID


## set outdir path
#$wrkdir outdir=/home/groups/hoolock/u1/jg/replicate/qc/output

# set infile directory path
#inpath=/home/groups/hoolock/u0/rawdata/Epigenetics_Core/ECP15

# set path to arguments file
#args_file=/home/groups/hoolock/u1/jg/replicate/qc/ECP15.txt

## JG 10/11/20: Added Single/Paired End input. I believe this script is agnostic to this fact. 

# set variable "ARGS" to be the output of the correct line from args_file (matches SLURM_ARRAY_TASK_ID)
ARGS=`head -$SLURM_ARRAY_TASK_ID $args_file | tail -1`

# set variables based on the output of the line ARGS
IFS=: read INFILE INFILE2 SEQ <<< $ARGS


#	TASKS:
#echo $INFILE
# Do something to each input file:
#$exafastqc -t 2 -o $wrkdir/qc $sourcedir/$INFILE*.fastq.gz 
## Check expected output: If it doesnt exist execute otherwise skip. 
if [ ! -f $wrkdir/qc/${INFILE}*R1*fastqc.html  ]; then 
    echo perform fastqc on $INFILE
    $exafastqc -t 8 -o $wrkdir/qc $sourcedir/$INFILE*.f*q.gz 
else
  echo output detected skipping
fi


