#!/bin/bash

##### PURPOSE #####
## Jacob Gutierrez gutierja@ohsu.edu
## JG 11/17/20: Now Adding Output file check 
## Reused 3/29/20
## The purpose of this file is to intitialize the probe validation program. PREVIOUSLY
# This script is now the dependecy file for the preprocessing of the raw methyl-capture sequencing data. 
# Analysis Pipeline:
## 0) Create Mapping File from excel: Python pandas
## 1) raw sequencing : Fastqc
## 2) read triming : Trimegalore (cutadapt + fastqc)
## 3) Alignment : Bismark (bowtie2)
## 4) Sorting & Indexing : Samtools

# Could consider manually setting the probe_bed file

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
    exit 1 # fail here? Or MAYBE try scancel 
    fi
done # end of task checking
## Print Success
echo success: Job OK

}
export -f check_success


## Build sbatch command 
format_sbatch() {

#outerrsbatch='sbatch --output=./out/'${expname}'.'$1'.%A_%a.out --error=./err/'${expname}'.'$1'.%A_%a.err '
  #echo 'sbatch --account '$slurm_acct' --job-name='${expname}'_'$1' --output=./slurmout/'${expname}'.'$1'.%A_%a.out --error=./slurmerr/'${expname}'.'$1'.%A_%a.err'
  echo 'sbatch --account '$slurm_acct' --job-name='$1' --output=./slurmout/'${expname}'.'$1'.%A_%a.out --error=./slurmerr/'${expname}'.'$1'.%A_%a.err'
  
}


####### END OF FUNCTIONS #####


##### Create Mapping File and Find Number of File (nodes) TO call#####
if [ ! -f $args_file ]; then 
    $script_dir/create_args.py $args_file $excel $expname
    
else
## CHECK THE SECOND VARIABLE FOR IT TO EXIST
ARGS=`head -1 $args_file | tail -1`

# set variables based on the output of the line ARGS
IFS=: read CHECK <<< $ARGS

fi


# Find number of files 
filenum=`cat $args_file | wc -l`



#### Raw FastQC #####
echo Raw FastQC
base_sbatch=`format_sbatch rawqc` # set prefix
rawqc=$($base_sbatch \
--array=1-${filenum}%${joblim} \
$script_dir/fastqc_batch.sh)

rawqc=${rawqc##* }
echo $rawqc




#### TimeGalore #####
echo Read Trimming
base_sbatch=`format_sbatch trim` # set prefix
trim=$($base_sbatch \
--dependency=afterany:${rawqc} \
--array=1-${filenum}%${joblim} \
$script_dir/trimgalore_batch.sh)


trim=${trim##* }
echo $trim

#exit 0 ## JG 10/11/20: Added this to debug the triming command. 

#### Alignment #####
echo Aligning to hg38 assembly 
base_sbatch=`format_sbatch align` # set prefix
align=$($base_sbatch \
--dependency=afterany:${trim} \
--array=1-${filenum}%${joblim} \
--qos=long_jobs \
$script_dir/bismark_align_batch.sh)


align=${align##* }
echo $align

#exit 0 ## JG 10/11/20: Added this to debug the align command. 


#### Post Alignment Processing #####
echo Post Alignment Processing 
base_sbatch=`format_sbatch post_align` # set prefix
post_align=$($base_sbatch \
--dependency=afterany:${align} \
--array=1-${filenum}%${joblim} \
--qos=long_jobs \
$script_dir/bismark_post_align_batch.sh)


post_align=${post_align##* }
echo $post_align


#### Sorting and Indexing #####
echo Sorting and Indexing aligned output 
base_sbatch=`format_sbatch sortid` # set prefix
sortid=$($base_sbatch \
--dependency=afterany:${post_align} \
--array=1-${filenum}%${joblim} \
$script_dir/sort_index_batch.sh)


sortid=${sortid##* }
echo $sortid


#### MultiQC report #####
echo MutltiQC report
base_sbatch=`format_sbatch multiqc` # set prefix
#report_name="preprocessing_report"
mqc=$($base_sbatch \
--dependency=afterany:${sortid} \
$script_dir/multiqc_batch.sh)

mqc=${mqc##* }
echo $mqc

#### Clean Excess sequencing files #####
# Remove large files from 1) trimming and 2) raw alignments
#echo Remove excess reads
#base_sbatch=`format_sbatch rmreads` # set prefix

#rmreads=$($base_sbatch \
#--dependency=afterany:${mqc} \
#$script_dir/remove_excess_reads.sh)




echo All Job Submitted