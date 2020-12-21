#!/bin/bash
## JG 12/16/20: 

#SBATCH -p exacloud                # partition (queue) (exacloud, long_jobs)
#SBATCH -N 2                       # number of nodes
#SBATCH -n 16                       # number of cores
#SBATCH -c 1                      # Request that ncpus be allocated per process. (--cpus-per-task=1 (or -c) tells slurm to make each )
#SBATCH --mem 24000                # memory pool for all cores (right now 24000MB = 24G) 12/15/20: modified to 
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


########## Processs Array Job Input #########

########## TASKS #########

######## CHANGED 7/10/20 ######## 
## Check for trimmed fastqc report! This is done in the fastqc check as well. 
catch=`ls $wrkdir/trim | grep "_fastqc.html$" | grep "^${INFILE}[_|\.]" |  tr ' ' '\n'`

## JG 10/23/20: modify script to only check one input... Consider checking all results. 
## See here: https://stackoverflow.com/questions/10586153/split-string-into-an-array-in-bash
catch=($(echo "$catch"))


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

