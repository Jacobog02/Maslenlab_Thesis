#!/bin/bash
## 10/22/20: Added ambiguous call. Consider parameterizing to allow multiple kinds of behavior. 
## JG 12/16/20: Modifying post_align_batch.sh to handle deduplication. 
## IN hindsight I must rewrite this script to do the sorting/indexing of both the reads and the ambiguous data. 
## Those functions need to know if deduplciation must occur so I will unite them as one. 
## JG 12/17/20: Looking at tutorial: https://genefish.wordpress.com/2019/09/23/code-bin-bash-job-name-sbatch-12/ 
## I will pass the bam output directly into deduplciate 

#SBATCH -p exacloud                # partition (queue) (exacloud, light, gpu) OLD see QOS long_jobs
#SBATCH -N 2                       # number of nodes
#SBATCH -n 24                       # number of cores
#SBATCH --mem 96000                # memory pool for all cores
#SBATCH -t 1-00:00                 # time (D-HH:MM)


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

# set variable "ARGS" to be the output of the correct line from args_file (matches SLURM_ARRAY_TASK_ID)
ARGS=`head -$SLURM_ARRAY_TASK_ID $args_file | tail -1`

## JG 10/11/20: Adding Single/Paired End Sequencing Formatting. 
# set variables based on the output of the line ARGS
#IFS=: read INFILE INFILE2 <<< $ARGS
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

### DEDUPLICATION PROCESS #####

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
    
    else 
      echo "Deduplicated output detected Now skipping"
      
    fi 
    
    
## IT is methyl capture    
 else
  echo $INFILE2 : $LIB "Library Detected"
  echo "Deduplication is not not applied to targeted sequencing."
  
fi   


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
  
else
  echo "Previously Made Bismark Coverage Reports Detected skipping step"
  
fi

catch=`ls $wrkdir/meth | grep ".nucleotide_stats.txt" | grep "^${INFILE2}[_|\.]"`

if [ ! -f $wrkdir/meth/$catch  ]; then 
  echo "Perform bam2nuc cytosine *.nucleotide_stats.txt"
  ## Create Nucleotide report for multiqc
  $bispath/bam2nuc \
  --dir $wrkdir/meth/ \
  --genome_folder $hg38 \
  $wrkdir/align/${INFILE2}${seq_suffix}${dedup_suffix}.bam
  
  
  #else
  #echo "Previously Made Bismark Coverage Reports Detected skipping step"
  
fi


### JG 12/17/20: Now adding sorting and indexing calls HERE
### START OF COPIED sort_index_batch.sh I WILL MODIFY TO FIT INTO THIS WRAPPER #####

## I will know as a human based on the WGBS status if the sorted indexed bam files were deduplicated or not. 
catch=`ls $wrkdir/reads | grep ".bam$" | grep "^srt.${INFILE2}[_|\.]"`


if [ ! -f $wrkdir/reads/$catch ]; then 

    echo "Sorting and Indexing Aligned" ${INFILE2}

    $exa_sam sort -o $wrkdir/reads/srt.${INFILE2}.bam -@ 12 -m 8G $wrkdir/align/${INFILE2}${seq_suffix}${dedup_suffix}.bam

    $exa_sam index -@ 8 $wrkdir/reads/srt.${INFILE2}.bam
    
    else
    echo "Previous Algined output detected now skipping"
fi


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
    
    else
    echo "Previous Ambiguous output detected now skipping"
fi



### END OF jgupdate_bismark_post_align_batch.sh script ####
## JG 12/17/20 Note: I think I can copy and paste the executables of the directory into one file to be run on one node for many days. 