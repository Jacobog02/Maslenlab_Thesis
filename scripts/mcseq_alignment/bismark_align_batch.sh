#!/bin/bash

#SBATCH -p exacloud                # partition (queue) (exacloud, light, gpu) OLD see QOS long_jobs
#SBATCH -N 2                       # number of nodes
#SBATCH -n 24                       # number of cores
#SBATCH --mem 96000                # memory pool for all cores
#SBATCH -t 3-00:00                 # time (D-HH:MM)


##### Checking Dependency Success #####
dependid=${SLURM_JOB_DEPENDENCY##*:}

deptasks=(`sacct -j $dependid --format State | grep -v [S/-] | tr -d ' ' `)

check_success ${deptasks[@]}

echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID

## Tool paths:
#exa_opt=/opt/installed/

# set outdir path
#outdir=/home/groups/hoolock/u1/jg/replicate/align/output

# set infile directory path
#inpath=/home/groups/hoolock/u1/jg/replicate/trim/data

# set path to arguments file
#args_file=/home/groups/hoolock/u1/jg/replicate/align/al_arg


## JG 10/11/20: Adding Single/Paired End Sequencing Formatting. 

# set variable "ARGS" to be the output of the correct line from args_file (matches SLURM_ARRAY_TASK_ID)
ARGS=`head -$SLURM_ARRAY_TASK_ID $args_file | tail -1`

## JG 10/11/20: Adding Single/Paired End Sequencing Formatting. 
# set variables based on the output of the line ARGS
#IFS=: read INFILE INFILE2 <<< $ARGS
IFS=: read INFILE INFILE2 SEQ <<< $ARGS

#	TASKS:
# Do something to each input file:

## Run bismark with default parameters. 
## Add file detection
## Make it look for a unique file ending
catch=`ls $wrkdir/align | grep ".bam$" | grep "^${INFILE2}[_|\.]"`

## If output doesnt exist run the alignment
if [ ! -f $wrkdir/align/$catch  ]; then 
  ## JG 10/10/20: Needing to add a Single End or Paired End Options: Check 
  ## Check Paired End
  if [ $SEQ == "PE" ]; then
    ## JG 10/11/20: Adding Single/Paired End Conditionals. 
    echo "Perform PAIRED end read alignment on" $INFILE
    ## Get Readnames
    ## Final grep for unique file changed 7/10/20 from grep -w "^${INFILE}" to grep "^${INFILE}[_|\.]" to match a _ or a . \.
    ## THis was done to stop it match S1 to S1,S10,S11 etc... 
    ## The -w option does not count underscore as a non-word character and doesnt match NAME_SAMPLE1_SOMETHING
    R1=`ls $wrkdir/trim/ | grep "_val_1.fq.gz" | grep "^${INFILE}[_|\.]"`
    R2=`ls $wrkdir/trim/ | grep "_val_2.fq.gz" | grep "^${INFILE}[_|\.]"`
    echo $R1
    echo $R2
    
    $bispath/bismark \
    -p 8 \
    $hg38 \
    --temp_dir $wrkdir/align \
    -1 $wrkdir/trim/$R1 \
    -2 $wrkdir/trim/$R2 \
    -o $wrkdir/align \
    --basename ${INFILE2} #\
    #--bam
  ## Else check Single End  
  elif [ $SEQ == "SE" ]; then 
    #echo "Perform SINGLE end read alignment on" $INFILE
    ## Get Readnames
    ## Final grep for unique file changed 7/10/20 from grep -w "^${INFILE}" to grep "^${INFILE}[_|\.]" to match a _ or a . \.
    ## THis was done to stop it match S1 to S1,S10,S11 etc... 
    ## The -w option does not count underscore as a non-word character and doesnt match NAME_SAMPLE1_SOMETHING
    R1=`ls $wrkdir/trim/ | grep "trimmed.fq.gz" | grep "^${INFILE}[_|\.]"`
    #R2=`ls $wrkdir/trim/ | grep "_val_2.fq.gz" | grep "^${INFILE}[_|\.]"`
    #echo $R1
    echo "Perform SINGLE end read alignment on" $R1
    #echo $R2
    
    $bispath/bismark \
    -p 8 \
    $hg38 \
    --temp_dir $wrkdir/align \
    -o $wrkdir/align \
    --basename ${INFILE2} \
    $wrkdir/trim/$R1
  
  else 
      echo ${INFILE} "Does not have proper Seq_Type only PE or SE accepted."
      exit 1 
    fi ## End of PE/SE Check
## REPORT ALIGNMENT DONE
else 
echo "Previous Alignment Output Detected now skipping" $INFILE

fi ## END OF CATCH CHECK


########## END TASKS #########
## Perform methyl extractor to obtain coverage report! 
## Just comment out for this alignment. 
