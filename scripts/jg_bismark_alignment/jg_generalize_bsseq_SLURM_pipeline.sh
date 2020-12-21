#!/bin/bash

##### PURPOSE #####
## Jacob Gutierrez gutierja@ohsu.edu
## JG 11/17/20: Now Adding Output file check 
## Reused 3/29/20
## The purpose of this file is to intitialize the probe validation program. PREVIOUSLY
# This script is now the dependecy file for the preprocessing of the raw methyl-capture sequencing data. 


## generlized_bsseq_SLURM_pipeline.sh Issues ###
## The main issue with my approach is that it doesnt check for the output until it is time to do it. 
## I realized I Should just deploy only 1 (ONE!!!) larger worker node for each file which does qc > indexing... 
## This is not what you want to deploy if your code has risk to fail. Because of this I just made individual scripts that spawn jobs. 
## BUT because I made them all use the sample template all of the <compute> sections can be directly computed in serial by each worker. 

## JG 12/20/20: I made the modification this script now calls only 1 massive file to do the full analysis from raw qc > read sorting and indexing. 
## The job is expected to take 6-7 (depending on the amount of data) days to finish per sample. 
## Modifying pipeline to make this a shallow wrapper to initalize the slurm environment and to allocate array jobs for each sample. 
## The jobs themselves will make the directories and have a new monitor script. 

# Analysis Pipeline:
## 0) Create Mapping File from excel: Python pandas
## 1) raw sequencing : Fastqc
## 2) read triming : Trimegalore (cutadapt + fastqc)
## 3) Alignment : Bismark (bowtie2)
## 4) Sorting & Indexing : Samtools

## JG 12/15/20: Converted this function from mcseq_SLURM_pipeline.sh >>> generalized_bsseq_SLURM_pipeline.sh in another directory. 
## I save the old one as a working example. 

## JG 12/14/20: Now adding deduplication features. But I will do this by adding another column in the input excel sheet.
## Additionally I will modify this pipeline to have a template BLANK_slurm_batch.sh so I can create new analyses easily. 
## I also need to harmonize the IFS:: call so that I only need to modify that once. I could even have a python program make dynamic memory passed into bash instead of a flat file. 

## Requires the Experiment setup portion to be specified (perhaps consider making command line 

## AND! a formatted file of of the source bams and the new file names to simplify. Manually formated. 

# Run by sbatch probe_pipeline.sh in a wrapper script that includes the experiment specific information.  

##### SBATCH INITALIZATION #####
#SBATCH -p exacloud                # partition (queue) (exacloud, long_jobs)
#SBATCH -N 1                       # number of nodes 
#SBATCH -n 2                      # number of cores
#SBATCH --mem 2400                # memory pool for all cores
#SBATCH -t 0-02:00                 # time (D-HH:MM)
#SBATCH -o ./slurmout/source.%A_%a.out           # Standard output
#SBATCH -e ./slurmerr/source.%A_%a.err           # Standard error


#### Experiment Set up #####
echo Depedency Script

#### Set Job limit for parallelized jobs#####
joblim=12
#### Exported Variables that live in the void #####

# Excel mapping file for all experiments
# excel #Excel mapping file.

# common name for experiment
# expname=ECP4 # String for the experiment name to be used for standard output notation no spaces andbest if its short. 

# Formatted Arguements file can be used for all processes just make it focus on new name 
## This is generated from the excel file!
# args_file=/home/groups/hoolock2/u0/jg/thesis/sa1/${expname}/${expname}_mapping.txt # Absolute path to argument files (colon delimited)

# Where the source raw fastq files live
# sourcedir=/home/groups/hoolock2/u0/archive/Epigenetics_Core_rawdata/ECP4 # Absolute path to the directory where the original fastq files are contained

# Where the output be placed 
#wrkdir=/home/groups/hoolock2/u0/jg/thesis/sa1/${expname}/ # Absolute path to the directory to be used for output writing 

## Where the neccesary scripts to run the pipeline live
# script_dir=/home/groups/hoolock2/u0/jg/thesis/scripts/ # Absolute path to directory for the pipeline scripts



##### Export Software and Functions #####

## JG 4/1/20: I started using a unified conda environment to simply function calls.

#export exa_opt=/opt/installed/ # Exacloud stuff but I will use my conda enviroment
#1) fastqc
#export exafastqc=${exa_opt}/fastqc/fastqc  
export exafastqc=/home/groups/hoolock2/u0/jg/thesis/env/methylseq/bin/fastqc # I can just say fastqc bc it is on the path

#2) Trimegalore
#export exa_trimg=/home/exacloud/lustre1/BioCoders/Applications/trim_galore
export exa_trimg=/home/groups/hoolock2/u0/jg/thesis/env/methylseq/bin/trim_galore
#export exa_cut=~/.local/bin/cutadapt

#3.1) MultiQC report for raw and trim
export exa_mutliqc=~/.local/bin/multiqc

#3) Bismark Alignment
# Bismark 
#export bispath=/home/groups/hoolock2/u0/jg/bin/bismark/Bismark-0.22.3 # Newly installed bismark
## Consider putting bismark in the environment 
export bispath=/home/groups/hoolock2/u0/jg/thesis/env/methylseq/bin/ #bismark

# Human Genome (CHECK TO SEE IF BISULFITE CONVERTED)
export hg38=/home/groups/hoolock2/u0/jg/thesis/human_genome


#4) Sort and Index 
## JG 11/17/20: AND BEDTOOLS!!! I already have it in the environment. 
#export exa_sam=/opt/installed/samtools-1.6/bin/samtools
export exa_sam=/home/groups/hoolock2/u0/jg/thesis/env/methylseq/bin/samtools
export methyl_bed=/home/groups/hoolock2/u0/jg/thesis/env/methylseq/bin/bedtools

## Adding Bigwig Conversion 
## JG 10/30/20 
## Please run install_bedGraphToBigWig.sh in scripts dir. 
export bedGraphToBigWig=/home/groups/hoolock2/u0/jg/thesis/scripts/mcseq_alignment/bedGraphToBigWig
#export fetchChromSizes=/home/groups/hoolock2/u0/jg/thesis/scripts/mcseq_alignment/fetchChromSizes
#export ChromSizes=/home/groups/hoolock2/u0/jg/thesis/annotation_info/hg38.chrom.size
export ChromSizes=/home/groups/hoolock2/u0/jg/thesis/sa1/ECP15/ECP15_genome.txt

##### Functions ####

## Check success function queries slurm based on passed jobids and ends tasks if a previous failure was detected.
## JG 10/14/20: Updating the function to echo that the task failed. Look into if can capture the job name.... 
check_success () {
#echo inside check success

#echo "$@"
## Loop over
for task in "$@"; do
    #echo $task 
    if [ $task != "COMPLETED" ]; then
    #echo stopping script here
    echo exiting: previous job failed
    exit 1 # fail here!!! Or MAYBE try scancel 
    fi
done # end of task checking
## Print Success
echo success: Job OK

}
export -f check_success ## export to env for job nodes to access


## JG 12/14/20: Now adding readin_ARG function to change read in command only once.
## This command will allow each worker node 
readin_ARG () {
  
#ARGS=$1 ## Catch args variable. 
#IFS=: read INFILE INFILE2 SEQ <<< $ARGS # double check if these variables are available to the worker nodes
#IFS=: read INFILE INFILE2 SEQ <<< $1

IFS=: read INFILE INFILE2 SEQ LIB <<< $1 ## JG 12/15/20: modified for LIB. 

## Be sure the template uses this call instead of the historic one. Simplfies many things. 

}
export -f readin_ARG ## export to env for job nodes to access

## JG 12/14/20: Parse_Seqtype function. Accepts a argument which is the SeqType column and then returns the corresponding suffix.
parse_SeqType () {
  #echo $1
  local THESEQ=$1 ## This is expected to be $SEQ from readin_ARG
  if [ $THESEQ == "SE" ]; then
    
    #suffix=".bam"
    seq_suffix=""
    echo Single end SEQ $seq_suffix
    #param="-s"
    
  elif [ $THESEQ == "PE" ]; then 
    
    #suffix="_pe.bam"
    seq_suffix="_pe"
    echo Paired end SEQ $seq_suffix
    #param="-p"
    
  else
     ## This field is not set correctly
     echo $THESEQ "Field is not PE or SE"
     exit 1 ## leave with anger
  fi

}
export -f parse_SeqType ## export to env for job nodes to access

## JG 12/14/20: Parse_Libtype function. Accepts an argument which is the LibType column (MethylCapture or WGBS) and returns the appropriate deduplication suffix.
parse_Libtype () {
  echo $1
  local THELIB=$1 ## This is expected to be $SEQ from readin_ARG
  ## If it is methylcapture do not do the deduplication step. 
  if [ $THELIB == "MethylCapture" ]; then
    
    ## Set 'global' (will be set on each worker node but should not overwrite any other worker node environment variables)
    DEDUP=false
    #suffix=".bam"
    dedup_suffix=""
    echo MethylCapture Sample $dedup_suffix
    #param="-s"
    
  elif [ $THELIB == "WGBS" ]; then 
    
    DEDUP=true
    #suffix="_pe.bam"
    dedup_suffix=".deduplicated"
    echo Whole Genome Bisfultie Seq $seq_suffix
    #param="-p"
    
  else
     ## This field is not set correctly
     echo $THELIB "Field is not MethylCapture or WGBS. Only types supported as of 12/15/20"
     exit 1 ## leave with anger 
  fi

}
export -f parse_Libtype ## export to env for job nodes to access

## JG 12/20/20: Making display of the extracted data. 
## This should be used within the out and the monitor. 
print_sample_data () {
  
  ## Sample data 
  printf "%s\n" "File Prefix: $INFILE"
  printf "%s\n" "Sample Name: $INFILE2"
  printf "%s\n" "Seq Type: $SEQ"
  printf "%s\n" "Lib Type: $LIB"
  
  
  ## Initalize seq_suffix for sequence type specific calls. 
  printf "%s\n" "Seq Suffix: $seq_suffix" ##  to be used for sequence type specific calls. 
  
  ## Initalize dedup_suffix using $LIB for MethyCapture vs. WGBS. (WBGS needs deduplication after alignment) 
  printf "%s\n" "Dedup Suffix: $dedup_suffix"
  printf "%s\n" "Dedup Flag: $DEDUP"
  
}
export -f print_sample_data 


jg_datetime (){
  
  #nowtime=`date "+%x %T"`
  #printf "%s" $nowtime
  date "+%Y-%m-%d %H:%M:%S"
}
export -f jg_datetime 


## JG 12/20/20: functional monitor output. Gives standard output in TSV format to quickly analyze results if I wanted. 
## I really just want this so I can see where each node is at individually. 
jg_format_monitor (){
  
  ## Format the monitor output.
  MYCUSTOMTAB='   '
  step_name=$1
  start_stop_runtimemess=$2
  the_time=`jg_datetime`
  #message=$3
  
  
  ## Print formatted message
  ## Print the first parts
  printf "%s\t%s\t%s %s\t%s\n" $step_name $start_stop_runtimemess $the_time "$3"
  #printf "%s\n" "$3"
  
}
export -f jg_format_monitor 

## Build sbatch command 
## JG 12/15/20: Trying to make this format_sbatch command be build incrementally. 
format_sbatch() {

#outerrsbatch='sbatch --output=./out/'${expname}'.'$1'.%A_%a.out --error=./err/'${expname}'.'$1'.%A_%a.err '
  #echo 'sbatch --account '$slurm_acct' --job-name='${expname}'_'$1' --output=./slurmout/'${expname}'.'$1'.%A_%a.out --error=./slurmerr/'${expname}'.'$1'.%A_%a.err'
  #echo 'sbatch --account '$slurm_acct' --job-name='$1' --output=./slurmout/'${expname}'.'$1'.%A_%a.out --error=./slurmerr/'${expname}'.'$1'.%A_%a.err'
  
  ## Attempting to make list that get appended to then printed.
  local -a sbatch_print_list=()
  
  ## Sbatch command plus account 
  local sbatchandaccount=`echo 'sbatch --account='$slurm_acct ` ##Format item 
  sbatch_print_list+=$sbatchandaccount ## Append to list
  
  
  local jobname=`echo ' --job-name='$1` ## Be sure to add spaces before each arg
  sbatch_print_list+=$jobname ## Append to list
  
  local outerr=`echo ' --output=./slurmout/'${expname}'.'$1'.%A_%a.out --error=./slurmerr/'${expname}'.'$1'.%A_%a.err'`
  sbatch_print_list+=$outerr
  
  echo $sbatch_print_list ## Echoing the list gives the correctly formated sbatch command. Gives more flexibilty to create smart pipelines. 
  
  
}


## JG 12/14/20: add_dependency function. Takes in an argument which can be empty or a dependency. It will output the proper formatted sbatch argument and append it to command. 



####### END OF FUNCTIONS #####


####### Preprocessing Input Parameters #######
## JG 12/15/20: Modified to build the args_file directly
## JG 12/17/20: GAVE ME A DAMN HEADACHE OVER HERE PUTTING BACK IN THE run_ file. 
#export args_file=${wrkdir}/${expname}_mapping # Absolute path to argument files (colon delimited)


##### Create Mapping File and Find Number of File (nodes) TO call####
## JG 12/17/20 expect args_file to be called in the run_ file. 
if [ ! -f $args_file ]; then 
    #echo INSIDE ARGS BUILD
    #echo args file $args_file
    #echo wrkdir $wrkdir
    $scriptdir/create_args.py $args_file $excel $expname
    
else
## CHECK THE SECOND VARIABLE FOR IT TO EXIST
ARGS=`head -1 $args_file | tail -1`

# set variables based on the output of the line ARGS
IFS=: read CHECK <<< $ARGS

fi


# Find number of files 
filenum=`cat $args_file | wc -l`


####### END OF Preprocessing Input Parameters #######


######### START BATCH JOBS SEQUENTIALLY #############


#### Testing and Debugging Script. 
## JG 12/16/20: I made this section of code to allow for incremental debugging of generic slurm array scripts. 
## The template can be copied into another file for a specific analysis which can be debugged uniquly from there. 
#echo JG 12-15-20 Testing

#base_sbatch=`format_sbatch test_slurm` # set prefix
#echo $base_sbatch
#slurm_test=$($base_sbatch \
#--array=1-${filenum}%${joblim} \
#$scriptdir/BLANK_template_batch.sh)

#slurm_test=${slurm_test##* }
#echo $slurm_test

#echo JG 12 15 20 Debug ending
#exit 0 ## JG 1l2/15/20 for debugging. 

#### Raw FastQC #####
## JG 12/15/20: Updated to match the template. 
# echo Raw FastQC
# base_sbatch=`format_sbatch rawqc` # set prefix
# rawqc=$($base_sbatch \
# --array=1-${filenum}%${joblim} \
# $scriptdir/fastqc_batch.sh)
# 
# rawqc=${rawqc##* }
# echo $rawqc
# 
# 
# sleep 1 ## Be kind to the scheduler!!!! 


#### TimeGalore #####
# echo Read Trimming
# base_sbatch=`format_sbatch trim` # set prefix
# trim=$($base_sbatch \
# --dependency=afterany:${rawqc} \
# --array=1-${filenum}%${joblim} \
# $scriptdir/trimgalore_batch.sh)
# 
# 
# trim=${trim##* }
# echo $trim
# 
# sleep 1 ## Be kind to the scheduler!!!! 

#exit 0 ## JG 10/11/20: Added this to debug the triming command. 

#### Alignment #####
# echo Aligning to hg38 assembly 
# base_sbatch=`format_sbatch align` # set prefix
# align=$($base_sbatch \
# --dependency=afterany:${trim} \
# --array=1-${filenum}%${joblim} \
# --qos=long_jobs \
# $scriptdir/bismark_align_batch.sh)
# 
# 
# align=${align##* }
# echo $align
# 
# sleep 1 ## Be kind to the scheduler!!!! 

#exit 0 ## JG 10/11/20: Added this to debug the align command. 


#### Post Alignment Processing #####
# echo Post Alignment Processing 
# base_sbatch=`format_sbatch post_align` # set prefix
# post_align=$($base_sbatch \
# --dependency=afterany:${align} \
# --array=1-${filenum}%${joblim} \
# --qos=long_jobs \
# $scriptdir/bismark_post_align_batch.sh)
# 
# 
# post_align=${post_align##* }
# echo $post_align
# 
# sleep 1 ## Be kind to the scheduler!!!! 


#### Sorting and Indexing #####
# echo Sorting and Indexing aligned output 
# base_sbatch=`format_sbatch sortid` # set prefix
# sortid=$($base_sbatch \
# --dependency=afterany:${post_align} \
# --array=1-${filenum}%${joblim} \
# $scriptdir/sort_index_batch.sh)
# 
# 
# sortid=${sortid##* }
# echo $sortid
# 
# sleep 1 ## Be kind to the scheduler!!!! 


#### MultiQC report #####
# echo MutltiQC report
# base_sbatch=`format_sbatch multiqc` # set prefix
# #report_name="preprocessing_report"
# mqc=$($base_sbatch \
# --dependency=afterany:${sortid} \
# $scriptdir/multiqc_batch.sh)
# 
# mqc=${mqc##* }
# echo $mqc

#### Clean Excess sequencing files #####
# Remove large files from 1) trimming and 2) raw alignments
#echo Remove excess reads
#base_sbatch=`format_sbatch rmreads` # set prefix

#rmreads=$($base_sbatch \
#--dependency=afterany:${mqc} \
#$script_dir/remove_excess_reads.sh)


##### JG 12/17/20 MegaJob Implementation. 
## To optimize resource utilization I just request a node with the maximum power I need as an array for files of interest. 
## Each node is responsible for only 1 file and uses the same initalization to do the complete analysis. 
## I left the incremental code as single callable files so you CAN DEVELOP NEW ANALYSES AND THEN MERGE THEM IN THE END. 
## I may want to look into a way to code this directly so I can merge a list of files to make a mega sbatch call. 


## Post Alignment Process 
# echo Post Aligment Processing: deduplication,methyl extractor, nucleotide stats, sorting indexing, ambiguous processing
# base_sbatch=`format_sbatch fullpostalign` # set prefix
# fullpostalign=$($base_sbatch \
# --array=1-${filenum}%${joblim} \
# $scriptdir/jgupdate_bismark_post_align_batch.sh)
# 
# fullpostalign=${fullpostalign##* }
# echo $fullpostalign


## FULL ALIGNMENT PROCESS
echo "Full QC,  Aligment, post alignment processing:"
echo "raw fasta qc, trimming, alignment, deduplication (if neccesary),methyl extractor, nucleotide stats, sorting indexing, ambiguous processing"
echo Aligning to hg38 assembly
base_sbatch=`format_sbatch full_bismark_align` # set prefix
full_bismark=$($base_sbatch \
--array=1-${filenum}%${joblim} \
--qos=long_jobs \
$scriptdir/jg_full_bismark_alignpostprocess_batch.sh)


full_bismark=${full_bismark##* }
echo $full_bismark

sleep 1 ## Be kind to the scheduler!!!!


#### MultiQC report #####
echo MutltiQC report
base_sbatch=`format_sbatch multiqc` # set prefix
#report_name="preprocessing_report"
mqc=$($base_sbatch \
--qos=light \
--dependency=afterany:${full_bismark} \
$scriptdir/multiqc_batch.sh)

mqc=${mqc##* }
echo $mqc



######### END OF START BATCH JOBS SEQUENTIALLY #############

echo All Job Submitted