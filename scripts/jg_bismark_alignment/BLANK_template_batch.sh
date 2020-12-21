#!/bin/bash
## JG 12/15/20: Template slurm script. I figured I should make a generic template that I could copy to make new parts of the pipeline. 

#SBATCH -p exacloud               # partition (queue) (exacloud, light, gpu) OLD see QOS long_jobs
#SBATCH -N 1                      # number of nodes (when N == 1 all cores are on 1 machines. )
#SBATCH -n 2                      # number of cores
#SBATCH -c 1                      # Request that ncpus be allocated per process. (--cpus-per-task=1 (or -c) tells slurm to make each )
#SBATCH --mem 2000                # memory pool for all cores ( 2000MB = 2 G )
#SBATCH -t 0-00:10                # time (D-HH:MM)


##### Checking Dependency Success #####
dependid=${SLURM_JOB_DEPENDENCY##*:} ## extract dependency if available

## Request the dependency task results from slurm! 
deptasks=(`sacct -j $dependid --format State | grep -v [S/-] | tr -d ' ' `)

## Use custom function to check if they succeeded. 
## The way it is coded now is that dependent jobs will still execute regardless of completion. 
## I read online this is the best way currently to make sure you are doing the analysis correctly. 
check_success ${deptasks[@]} ## Exits the job if any previous tasks failed. 

echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID

## Tool paths:
## JG 12/15/20: Historically you could define the tools you needed directly in each script. 
## I used the arguably risky approach to throw these tool variables exported into the global environment. 

## JG 10/11/20: Adding Single/Paired End Sequencing Formatting. 

# set variable "ARGS" to be the output of the correct line from args_file (matches SLURM_ARRAY_TASK_ID)
ARGS=`head -$SLURM_ARRAY_TASK_ID $args_file | tail -1`

## JG 10/11/20: Adding Single/Paired End Sequencing Formatting. 
# set variables based on the output of the line ARGS
#IFS=: read INFILE INFILE2 <<< $ARGS
#IFS=: read INFILE INFILE2 SEQ <<< $ARGS

#	TASKS:
# Do something to each input file:

## Parse the colon delimited mapping files with unified call. 
## Standardize how my workers get the data based on the jobID == row of mapping file. 
## This is a naive approach but worked for my intial pipeline and I am committed. 
readin_ARG $ARGS
echo $INFILE
echo $INFILE2
echo $SEQ
echo $LIB


## Initalize seq_suffix for sequence type specific calls. 
echo SeqType Suffix Parse
parse_SeqType $SEQ ## Declare seq_suffix and Get the appropriate seq suffix "_pe" or nothing for paired end or single end. 
echo $seq_suffix ##  to be used for sequence type specific calls. 


## Initalize dedup_suffix using $LIB for MethyCapture vs. WGBS. (WBGS needs deduplication after alignment) 
echo deduplication parse
parse_Libtype $LIB
echo $dedup_suffix
echo $DEDUP


##### Run Code of Interest Below ######

exit 0 ## end forcefully for debugging 


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
    
    ## 10/22/20: Added ambiguous call. Consider parameterizing to allow multiple kinds of behavior. 
    $bispath/bismark \
    -p 8 \
    $hg38 \
    --temp_dir $wrkdir/align \
    -1 $wrkdir/trim/$R1 \
    -2 $wrkdir/trim/$R2 \
    -o $wrkdir/align \
    --basename ${INFILE2} \
    --ambig_bam --ambiguous 
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
    $wrkdir/trim/$R1 \
    --ambig_bam --ambiguous 
  
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
