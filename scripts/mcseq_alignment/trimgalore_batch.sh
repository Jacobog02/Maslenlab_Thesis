#!/bin/bash

#SBATCH -p exacloud                # partition (queue) (exacloud, long_jobs)
#SBATCH -N 2                       # number of nodes
#SBATCH -n 15                       # number of cores
#SBATCH --mem 24000                # memory pool for all cores
#SBATCH -t 00-06:00                  # time (D-HH:MM)

##### Checking Dependency Success #####
dependid=${SLURM_JOB_DEPENDENCY##*:}

deptasks=(`sacct -j $dependid --format State | grep -v [S/-] | tr -d ' ' `)

check_success ${deptasks[@]}

##### RUN INFO #####
echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID

########## Processs Array Job Input #########
## JG 10/11/20: Added Single/Paired End sequencing parameters within input.
## A single input file contains the seq type within the file prefix. This is bc a batch could contain samples with both types. 
## I would rather have the pipeline built to be agnostic of the file types per batch. 

# set variable "ARGS" to be the output of the correct line from args_file (matches SLURM_ARRAY_TASK_ID)
ARGS=`head -$SLURM_ARRAY_TASK_ID $args_file | tail -1`

# set variables based on the output of the line ARGS
#IFS=: read INFILE INFILE2 <<< $ARGS
IFS=: read INFILE INFILE2 SEQ <<< $ARGS


########## Processs Array Job Input #########

########## TASKS #########
# Do something to each input file:
## Check expected output: If it doesnt exist execute otherwise skip. 
## Check for trimmed reads
#catch=`ls $wrkdir/trim | grep "_val_1.fq.gz$" | grep "^${INFILE}[_|\.]"`

######## CHANGED 7/10/20 ######## 
## Check for trimmed fastqc report! This is done in the fastqc check as well. 
catch=`ls $wrkdir/trim | grep "_fastqc.html$" | grep "^${INFILE}[_|\.]"`



if [ ! -f $wrkdir/trim/${catch} ]; then 

  ## JG 10/10/20: Needing to add a Single End or Paired End Options: Check 
  ## Check Paired End
  if [ $SEQ == "PE" ]; then
    echo "Perform PAIRED end read trimming on" $INFILE
    ## Get Readnames
     ## Final grep for unique file changed 7/10/20 from grep -w "^${INFILE}" to grep "^${INFILE}[_|\.]" to match a _ or a . \.
    ## THis was done to stop it match S1 to S1,S10,S11 etc... 
    ## The -w option does not count underscore as a non-word character and doesnt match NAME_SAMPLE1_SOMETHING
    R1=`ls $sourcedir | grep ".f*q.gz$" | grep "R1" | grep "^${INFILE}[_|\.]"`
    R2=`ls $sourcedir | grep ".f*q.gz$" | grep "R2" | grep "^${INFILE}[_|\.]"`
    echo $R1
    echo $R2
    
    
    $exa_trimg \
    -j 4 \
    -o $wrkdir/trim \
    --fastqc \
    --fastqc_args "-t 4 -o ${wrkdir}/trim" \
    --paired  ${sourcedir}/${R1} ${sourcedir}/${R2}
    
  ## Else check Single End  
  elif [ $SEQ == "SE" ]; then 
    echo "Perform SINGLE end read trimming on" $INFILE
    ## Get Readnames
     ## Final grep for unique file changed 7/10/20 from grep -w "^${INFILE}" to grep "^${INFILE}[_|\.]" to match a _ or a . \.
    ## THis was done to stop it match S1 to S1,S10,S11 etc... 
    ## The -w option does not count underscore as a non-word character and doesnt match NAME_SAMPLE1_SOMETHING
    R1=`ls $sourcedir | grep ".f*q.gz$" | grep "R1" | grep "^${INFILE}[_|\.]"`
    #R2=`ls $sourcedir | grep ".f*q.gz$" | grep "R2" | grep "^${INFILE}[_|\.]"`
    echo $R1
    #echo $R2
    
    
    $exa_trimg \
    -j 4 \
    -o $wrkdir/trim \
    --fastqc \
    --fastqc_args "-t 4 -o ${wrkdir}/trim" \
    ${sourcedir}/${R1}
    #--paired  ${sourcedir}/${R1} ${sourcedir}/${R2}
    ## or else just report the mapping didnt work and exit badly
    else 
      echo ${INFILE} "Does not have proper Seq_Type only PE or SE accepted."
      exit 1 
    fi ## End of PE/SE Check
    
  else 
  echo "Previous Trimming Output Detected now skipping" $INFILE

fi

########## END TASKS #########

