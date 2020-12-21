#!/bin/bash
## JG 12/15/20: Code was refactored to be derived from the BLANK_template_batch.sh script. 
## Was orignally sourced from the previous version without deduplication. 

#SBATCH -p exacloud                # partition (queue) (exacloud, long_jobs)
#SBATCH -N 1                       # number of nodes
#SBATCH -n 16                       # number of cores
#SBATCH -c 1                      # Request that ncpus be allocated per process. (--cpus-per-task=1 (or -c) tells slurm to make each )
#SBATCH --mem 2000              # 36 GB memory pool for all cores
#SBATCH -t 00-02:00                 # time (D-HH:MM)

echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID


## JG 12/16/20: I believe this will work even in nodes where they have no dependencies.... This is because check_success actually looked for failed states lol... I need to rename it. 
## Just double checked and check_success searches for dependency ids and makes sure each id is "COMPLETED". If it is empty it skips this check and says it is fine! 
## I Think this can be generically used for my purposes. 
##### Checking Dependency Success #####

dependid=${SLURM_JOB_DEPENDENCY##*:} ## extract dependency if available

## Request the dependency task results from slurm! 
deptasks=(`sacct -j $dependid --format State | grep -v [S/-] | tr -d ' ' `)

## Use custom function to check if they succeeded. 
## The way it is coded now is that dependent jobs will still execute regardless of completion. 
## I read online this is the best way currently to make sure you are doing the analysis correctly. 
check_success ${deptasks[@]} ## Exits the job if any previous tasks failed. 


## JG 10/11/20: Added Single/Paired End input. I believe this script is agnostic to this fact. 

# set variable "ARGS" to be the output of the correct line from args_file (matches SLURM_ARRAY_TASK_ID)
ARGS=`head -$SLURM_ARRAY_TASK_ID $args_file | tail -1`

# set variables based on the output of the line ARGS
readin_ARG $ARGS
echo File Prefix: $INFILE
echo Sample Name: $INFILE2
echo Seq Type: $SEQ
echo Lib Type: $LIB

## Initalize seq_suffix for sequence type specific calls. 
echo SeqType Suffix Parse
parse_SeqType $SEQ ## Declare seq_suffix and Get the appropriate seq suffix "_pe" or nothing for paired end or single end. 
echo $seq_suffix ##  to be used for sequence type specific calls. 


## Initalize dedup_suffix using $LIB for MethyCapture vs. WGBS. (WBGS needs deduplication after alignment) 
echo deduplication parse
parse_Libtype $LIB
echo $dedup_suffix
echo $DEDUP


#	TASKS:

## Poorly optimized fastqc call. It just does the same analysis with same resources for both SE or PE.... again could be optimized. 

## Check expected output: If it doesnt exist execute otherwise skip. 
if [ ! -f $wrkdir/qc/${INFILE}*R1*fastqc.html  ]; then 
    echo perform fastqc on $INFILE
    $exafastqc -t 8 -o $wrkdir/qc $sourcedir/$INFILE*.f*q.gz 
else
  echo output detected skipping
fi


