#!/bin/bash
## JG 12/16/20: Modifying post_align_batch.sh to handle deduplication. 
## IN hindsight I must rewrite this script to do the sorting/indexing of both the reads and the ambiguous data. 
## Those functions need to know if deduplciation must occur so I will unite them as one. 

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
  $wrkdir/align/${INFILE2}${suffix}
  
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
  $wrkdir/align/${INFILE2}${suffix}
  
  
  #else
  #echo "Previously Made Bismark Coverage Reports Detected skipping step"
  
fi

## Perform methyl extractor to obtain coverage report! 
## Just comment out for this alignment. 
