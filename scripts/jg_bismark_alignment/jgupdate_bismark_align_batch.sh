#!/bin/bash
## JG 12/17/20: Modifiying to apply bismark deduplication (FILE COPIED FROM bismark_align_batch.sh)
## Key change is requiring SAM output at this step. I may want to just make the fastqc/trim step happen within this script.
## Looking at the tutorial: https://genefish.wordpress.com/2019/09/23/code-bin-bash-job-name-sbatch-12/ (I am blessed I found a slurm one that is doing my approach)
## I will still use the bam output from the version I have. It is very parallelized and works like a charm :)  SAM output breaks everything so figure it out!
## The person in the tutorial just gave the bam file into the deduplication function so... Ill try? ls
## Alter this script (or copy) to run for 6 days and do all 3 steps in one call. 
## I can check this worked and then do the post alignment call. 


#SBATCH -p exacloud                # partition (queue) (exacloud, light, gpu) OLD see QOS long_jobs
#SBATCH -N 2                       # number of nodes
#SBATCH -n 24                       # number of cores
#SBATCH --mem 120000                # memory pool for all cores ## 120G memory to use for alignment (and my other processes. Make it beefy)
#SBATCH -t 5-00:00                 # time (D-HH:MM) ## JG 12/17/20: Now runs for 5 days although only 4 are needed probably (some alignments take may days)


##### Checking Dependency Success #####
dependid=${SLURM_JOB_DEPENDENCY##*:}

deptasks=(`sacct -j $dependid --format State | grep -v [S/-] | tr -d ' ' `)

check_success ${deptasks[@]}

echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID


## JG 10/22/10: Adding Ambiguous Reads AND Alignment output. This might be integrated for SA1.
## I will use the random sample of reads to see if there is enrichment for targeted regions etc. 

## JG 10/11/20: Adding Single/Paired End Sequencing Formatting. 

# set variable "ARGS" to be the output of the correct line from args_file (matches SLURM_ARRAY_TASK_ID)
ARGS=`head -$SLURM_ARRAY_TASK_ID $args_file | tail -1`

## JG 10/11/20: Adding Single/Paired End Sequencing Formatting. 
# set variables based on the output of the line ARGS
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
## JG 12/17/20: I will now perform all the tasks in 1 large script. This will allow me to use maximum resource I need for all the tasks. 


#### Raw FASTQC
## Check to see if the html file exists for each file. If it doesnt run my hodgepodge code that runs all the INFILE matches.

#if [ ! -f $wrkdir/qc/${INFILE}*R1*fastqc.html  ]; then 
    #echo perform fastqc on $INFILE
    #$exafastqc -t 8 -o $wrkdir/qc $sourcedir/$INFILE*.f*q.gz 
#else
  #echo output detected skipping
#fi


######## CHANGED 7/10/20 ######## 
## Check for trimmed fastqc report! This is done in the fastqc check as well. 
catch=`ls $wrkdir/qc | grep "_fastqc.html$" | grep "^${INFILE}[_|\.]" |  tr ' ' '\n'`

## JG 10/23/20: modify script to only check one input... Consider checking all results. 
## See here: https://stackoverflow.com/questions/10586153/split-string-into-an-array-in-bash
catch=($(echo "$catch"))
if [ ! -f $wrkdir/qc/${catch}  ]; then 
    echo perform fastqc on $INFILE
    $exafastqc -t 8 -o $wrkdir/qc $sourcedir/$INFILE*.f*q.gz 
else
  echo output detected skipping
fi



#### Read Trimming & QC
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


## Run bismark with default parameters. 
## Add file detection
## Make it look for a unique file ending
## OLD CATCH!!
#catch=`ls $wrkdir/align | grep ".bam$" | grep "^${INFILE2}${seq_suffix}[_|\.]"`
## NEW CATCH I MADE THIS IT EXACTLY MATCHES the base name no ambig or deduplicated for now. 
catch=`ls $wrkdir/align | grep ".bam$" | grep "^${INFILE2}${seq_suffix}[_|\.]*bam$"`


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
