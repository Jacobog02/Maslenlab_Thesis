#!/bin/bash
## JG 12/17/20: Modifiying to apply bismark deduplication (FILE COPIED FROM bismark_align_batch.sh)
## Key change is requiring SAM output at this step. I may want to just make the fastqc/trim step happen within this script.
## Looking at the tutorial: https://genefish.wordpress.com/2019/09/23/code-bin-bash-job-name-sbatch-12/ (I am blessed I found a slurm one that is doing my approach)
## I will still use the bam output from the version I have. It is very parallelized and works like a charm :)  SAM output breaks everything so figure it out!
## The person in the tutorial just gave the bam file into the deduplication function so... Ill try? ls
## Alter this script (or copy) to run for 6 days and do all 3 steps in one call. 
## I can check this worked and then do the post alignment call. 

## JG 12/20/20: Modifying to create needed directories on the fly as one way to monitor the success of the job. 
## I will additionally create a monitor directory so each job is pumping out my custom made updates to monitor individual sample progress. 


#SBATCH -p exacloud               # partition (queue) (exacloud, light, gpu) OLD see QOS long_jobs this is set in the wrapper script. 
#SBATCH --qos long_jobs           # This is generally called in the wrapper but I put it here to be safe. 
#SBATCH -N 2                      # number of nodes
#SBATCH -n 24                     # number of cores
#SBATCH --mem 120000              # memory pool for all cores ## 120G memory to use for alignment (and my other processes. Make it beefy)
#SBATCH -t 7-00:00                # time (D-HH:MM) ## JG 12/17/20: Now runs for 5 days although only 4 are needed probably (some alignments take may days)


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
#echo File Prefix: $INFILE
#echo Sample Name: $INFILE2
#echo Seq Type: $SEQ
#echo Lib Type: $LIB

## Initalize seq_suffix for sequence type specific calls. 
#echo SeqType Suffix Parse
parse_SeqType $SEQ ## Declare seq_suffix and Get the appropriate seq suffix "_pe" or nothing for paired end or single end. 
#echo $seq_suffix ##  to be used for sequence type specific calls. 


## Initalize dedup_suffix using $LIB for MethyCapture vs. WGBS. (WBGS needs deduplication after alignment) 
#echo deduplication parse
parse_Libtype $LIB
#echo $dedup_suffix
#echo $DEDUP

## PRINT THE DATA TO THE OUTFILE. 
print_sample_data

###########################################################	TASKS ###########################################################	
## JG 12/17/20: I will now perform all the tasks in 1 large script. This will allow me to use maximum resource I need for all the tasks. 

## setup monitoring script
mkdir -p monitor

sample_monitor_f=${wrkdir}/monitor/${INFILE2}.$SLURM_JOBID

## Start the monitor file!!!
## give it the same output we see here
print_sample_data > $sample_monitor_f



## Adding header first. 
printf "%s\t%s\t%s\t%s\n" "Step_Name" "start_stop_runtime" "Datetime" "Message" >> $sample_monitor_f

## Starting slurmjob time for each sample
cur_step=SLURM_JOB
cur_strstop=START
#a_message=`echo -e "" "new slurm job for bisulfite sequencing mapping for ${INFILE2} ${expname}"`
a_message="Bisulfite sequencing mapping for slurmjob: ${INFILE2} ${expname}"
jg_format_monitor $cur_step $cur_strstop "$a_message" >> $sample_monitor_f



#################### Raw FASTQC ############################
## Check to see if the html file exists for each file. If it doesnt run my hodgepodge code that runs all the INFILE matches.

#if [ ! -f $wrkdir/qc/${INFILE}*R1*fastqc.html  ]; then 
    #echo perform fastqc on $INFILE
    #$exafastqc -t 8 -o $wrkdir/qc $sourcedir/$INFILE*.f*q.gz 
#else
  #echo output detected skipping
#fi


################ CHANGED 12/20/20 ################
## Adding mkdir call here
mkdir -p qc

## Adding monitor call here. 
cur_step=QC
jg_format_monitor $cur_step START "Raw Fastq reads FASTQC" >> $sample_monitor_f

######## CHANGED 7/10/20 ######## 
## Check for trimmed fastqc report! This is done in the fastqc check as well. 
catch=`ls $wrkdir/qc | grep "_fastqc.html$" | grep "^${INFILE}[_|\.]" |  tr ' ' '\n'`

## JG 10/23/20: modify script to only check one input... Consider checking all results. 
## See here: https://stackoverflow.com/questions/10586153/split-string-into-an-array-in-bash
catch=($(echo "$catch"))
if [ ! -f $wrkdir/qc/${catch}  ]; then 
    echo perform fastqc on $INFILE
    $exafastqc -t 8 -o $wrkdir/qc $sourcedir/$INFILE*.f*q.gz 
    jg_format_monitor $cur_step STOP "RUN" >> $sample_monitor_f
else
  echo output detected skipping
  jg_format_monitor $cur_step STOP "OUT DETECTED" >> $sample_monitor_f
fi



################################ Read Trimming & QC################################################
######## CHANGED 12/20/20 ########
## adding mkdir call 
mkdir -p trim

## Adding monitor call
cur_step=TRIM
jg_format_monitor $cur_step START "Raw Read Trimming TRIMEGALORE" >> $sample_monitor_f


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
    
    ## JG 12/20/20: monitor out
    jg_format_monitor $cur_step STOP "RUN" >> $sample_monitor_f
    
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
    
    jg_format_monitor $cur_step STOP "RUN" >> $sample_monitor_f

    
    else 
      echo ${INFILE} "Does not have proper Seq_Type only PE or SE accepted."
      exit 1 
    fi ## End of PE/SE Check
    
  else 
  echo "Previous Trimming Output Detected now skipping" $INFILE
    jg_format_monitor $cur_step STOP "OUT DETECTED" >> $sample_monitor_f


fi



################################ BISMARK ALIGMENT ################################################
## Run bismark with default parameters. 
## Add file detection
## Make it look for a unique file ending
######## CHANGED 12/20/20 ########
## adding mkdir call 
mkdir -p align

## Adding monitor call
cur_step=ALIGN
jg_format_monitor $cur_step START "Read Alignment w/ bismark" >> $sample_monitor_f

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
    
    ## JG 12/20/20: monitor
    jg_format_monitor $cur_step STOP "RUN" >> $sample_monitor_f

    
    
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
    
    jg_format_monitor $cur_step STOP "RUN" >> $sample_monitor_f

  
  else 
      echo ${INFILE} "Does not have proper Seq_Type only PE or SE accepted."
      exit 1 
    fi ## End of PE/SE Check
## REPORT ALIGNMENT DONE
else 
echo "Previous Alignment Output Detected now skipping" $INFILE
jg_format_monitor $cur_step STOP "OUT DETECTED" >> $sample_monitor_f


fi ## END OF CATCH CHECK



################################ DEDUPLICATION PROCESS ################################################
######## CHANGED 12/20/20 ########
## adding mkdir call 
#mkdir -p trim # not needed for dedup uses align. 

## Adding monitor call
cur_step=DEDUP
jg_format_monitor $cur_step START "checking deduplication status" >> $sample_monitor_f

## Check if we need to do deduplication
## JG 12/17/20: Modifying the if to check if the file exists as well. 
if [ "$DEDUP" = true ]; then
    echo $INFILE2 : $LIB "Library Detected " 


    ## Now check if file exists. 
    catch=`ls $wrkdir/align | grep ".deduplicated.bam$" | grep "^${INFILE2}[_|\.]"`
    
    if [ ! -f $wrkdir/meth/$catch  ]; then 
    
    echo "No Output Detected, Now performing deduplication"
    
    ## RUN DEDUPLICATION!
    ## TEST 3! Passing bam directly as input
    ## JG 12/17/20: THIS WORKED! GENERAL SO I DO NOT NEED PE VS SE. 
    $bispath/deduplicate_bismark --bam --output_dir $wrkdir/align $wrkdir/align/${INFILE2}${seq_suffix}.bam 
    ## Now all downstream calls use `${INFILE2}${seq_suffix}${dedup_suffix}.bam`
    
    jg_format_monitor $cur_step STOP "RUN" >> $sample_monitor_f

    
    else 
      echo "Deduplicated output detected Now skipping"
      jg_format_monitor $cur_step STOP "OUT DETECTED" >> $sample_monitor_f

      
    fi 
    
    
## IT is methyl capture    
 else
  echo $INFILE2 : $LIB "Library Detected"
  echo "Deduplication is not not applied to targeted sequencing."
  jg_format_monitor $cur_step STOP "MethylCapture detected skipping duplication" >> $sample_monitor_f

  
fi



################################ Methylation Extractor ################################################
######## CHANGED 12/20/20 ########
## adding mkdir call 
mkdir -p meth 

## Adding monitor call
cur_step=METH_EXTRACT
jg_format_monitor $cur_step START "getting CpG calling w/ bismark methylation extractor" >> $sample_monitor_f


## Add file detection
## Make it look for a unique file ending
#catch=`ls $wrkdir/meth | grep "_pe.bismark.cov.gz" | grep "^${INFILE2}[_|\.]"`
catch=`ls $wrkdir/meth | grep "\.bismark.cov.gz" | grep "^${INFILE2}[_|\.]"`

if [ ! -f $wrkdir/meth/$catch  ]; then 

  echo "Perform Methylation Extraction & Summary Report Generation" $INFILE2
  
  ## Parallel is set to 8 bc there are 3 threads in this process so 8*3 is 24 cores. 
  $bispath/bismark_methylation_extractor \
  $param --gzip --bedGraph \
  -o $wrkdir/meth/ \
  --parallel 8 \
  $wrkdir/align/${INFILE2}${seq_suffix}${dedup_suffix}.bam
  jg_format_monitor $cur_step STOP "RUN" >> $sample_monitor_f

  
else
  echo "Previously Made Bismark Coverage Reports Detected skipping step"
  jg_format_monitor $cur_step STOP "OUT DETECTED" >> $sample_monitor_f

  
fi


################################ NucleoTide Stats ################################################
######## CHANGED 12/20/20 ########
## adding mkdir call 
#mkdir -p meth # not needed for NUC_STATS uses meth. 

## Adding monitor call
cur_step=NUC_STATS
jg_format_monitor $cur_step START "getting nucleotide statistics for qc" >> $sample_monitor_f



catch=`ls $wrkdir/meth | grep ".nucleotide_stats.txt" | grep "^${INFILE2}[_|\.]"`

if [ ! -f $wrkdir/meth/$catch  ]; then 
  echo "Perform bam2nuc cytosine *.nucleotide_stats.txt"
  ## Create Nucleotide report for multiqc
  $bispath/bam2nuc \
  --dir $wrkdir/meth/ \
  --genome_folder $hg38 \
  $wrkdir/align/${INFILE2}${seq_suffix}${dedup_suffix}.bam
  
  jg_format_monitor $cur_step STOP "RUN" >> $sample_monitor_f

  
  else
  echo "Previously Made Nucleotide Stats Detected skipping step"
  jg_format_monitor $cur_step STOP "OUT DETECTED" >> $sample_monitor_f
  
fi



################################ Sorting Indexing ################################################
######## CHANGED 12/20/20 ########
## adding mkdir call 
mkdir -p reads # not needed for NUC_STATS uses meth. 

## Adding monitor call
cur_step=SORT_INDEX
jg_format_monitor $cur_step START "Sorting indexing appropriate reads" >> $sample_monitor_f

### JG 12/17/20: Now adding sorting and indexing calls HERE
### START OF COPIED sort_index_batch.sh I WILL MODIFY TO FIT INTO THIS WRAPPER #####

## I will know as a human based on the WGBS status if the sorted indexed bam files were deduplicated or not. 
catch=`ls $wrkdir/reads | grep ".bam$" | grep "^srt.${INFILE2}[_|\.]"`


if [ ! -f $wrkdir/reads/$catch ]; then 

    echo "Sorting and Indexing Aligned" ${INFILE2}

    $exa_sam sort -o $wrkdir/reads/srt.${INFILE2}.bam -@ 12 -m 8G $wrkdir/align/${INFILE2}${seq_suffix}${dedup_suffix}.bam

    $exa_sam index -@ 8 $wrkdir/reads/srt.${INFILE2}.bam
    
    jg_format_monitor $cur_step STOP "RUN" >> $sample_monitor_f

    else
    echo "Previous Algined output detected now skipping"
      jg_format_monitor $cur_step STOP "OUT DETECTED" >> $sample_monitor_f

fi


################################ Ambiguous Data Cleaning ################################################
######## CHANGED 12/20/20 ########
## adding mkdir call 
mkdir -p ambig # not needed for NUC_STATS uses meth. 

## Adding monitor call
cur_step=AMBIG
jg_format_monitor $cur_step START "Cleaning Ambiguous reads/alignments. sorting indexing appropriate reads" >> $sample_monitor_f


## JG 10/24/20: Adding Ambiguous Data Moving and Cleaning. 
## This section will check to see if the data is in the final location. If not move it (alignments AND input reads) to the right place and index. 
## JG 11/17/20: My desired output file is a bed version of the ambig.bam files. 
catch=`ls $wrkdir/ambig | grep ".ambig.bed$" | grep "^srt.${INFILE2}[_|\.]"`


if [ ! -f $wrkdir/ambig/$catch ]; then 

    echo "Sorting and Indexing Ambiguous Alignments and moving reads. " ${INFILE2}
    
    ## Copy the ambiguous reads to the ambig directory. 
    cp $wrkdir/align/${INFILE2}*ambiguous_reads*  $wrkdir/ambig/

    ## Sort and Index the Ambiguous Reads. 
    $exa_sam sort -o $wrkdir/ambig/srt.${INFILE2}.ambig.bam -@ 12 -m 8G $wrkdir/align/${INFILE2}${suffix}.ambig.bam #_pe.bam

    $exa_sam index -@ 8 $wrkdir/ambig/srt.${INFILE2}.ambig.bam
    
    ## JG 11/17/20: bamtobed conversion
    $methyl_bed bamtobed -i $wrkdir/ambig/srt.${INFILE2}.ambig.bam > $wrkdir/ambig/${INFILE2}.ambig.bed
    
    jg_format_monitor $cur_step STOP "RUN" >> $sample_monitor_f

    
    else
    echo "Previous Ambiguous output detected now skipping"
    jg_format_monitor $cur_step STOP "OUT DETECTED" >> $sample_monitor_f

fi



### END OF jgupdate_bismark_post_align_batch.sh script ####
## JG 12/17/20 Note: I think I can copy and paste the executables of the directory into one file to be run on one node for many days. 

