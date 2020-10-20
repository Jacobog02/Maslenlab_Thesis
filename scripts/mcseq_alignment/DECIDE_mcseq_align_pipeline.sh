#!/bin/bash

##### PURPOSE #####
## Jacob Gutierrez guteirja@ohsu.edu
## Reused 3/29/20
## The purpose of this file is to intitialize the probe validation program. PREVIOUSLY
# This script is now the dependecy file for the preprocessing of the raw methyl-capture sequencing data. 
# Analysis Pipeline:
## 0) Create Mappign File from excel: Python pandas
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
#SBATCH -o ./out/source.%A_%a.out           # Standard output
#SBATCH -e ./err/source.%A_%a.err           # Standard error


#### Experiment Set up #####
echo Depedency Script

#### Set Job limit for parallelized jobs#####
joblim=24
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

#3) Bismark Alignment
# Bismark 
export bispath=/home/groups/hoolock2/u0/jg/thesis/env/methylseq/bin/bismark # Newly installed bismark
## Consider putting bismark in the environment 

# Human Genome (CHECK TO SEE IF BISULFITE CONVERTED)
export hg38=/home/groups/hoolock2/u0/jg/thesis/human_genome


#4) Sort and Index.
export exa_sam=/opt/installed/samtools-1.6/bin/samtools

##### Functions ####

## Check success function queries slurm based on passed jobids and ends tasks if a previous failure was detected.
check_success () {
echo inside check success

#echo "$@"

for task in "$@"; do
    #echo $task 
    if [ $task != "COMPLETED" ]; then
    #echo stopping script here
    exit 1 # fail here? Or MAYBE try scancel 
    fi
done 
}
export -f check_success


## Build sbatch command 
format_sbatch() {

#outerrsbatch='sbatch --output=./out/'${expname}'.'$1'.%A_%a.out --error=./err/'${expname}'.'$1'.%A_%a.err '
  echo 'sbatch --account '$slurm_acct' --job-name='${expname}'_'$1'%'${joblim}' --output=./slurmout/'${expname}'.'$1'.%A_%a.out --error=./slurmerr/'${expname}'.'$1'.%A_%a.err'
  
}

#format_sbatch rawrxd


####### END OF FUNCTIONS #####



###### Base SBATCH function #####
#outerrsbatch='sbatch --output=./slurmout/'${expname}'.'`${delim}`'.%A_%a.out --error=./slurmerr/'${expname}'.'`${delim}`'.%A_%a.err '


##### Create Mapping File and Find Number of File (nodes) TO call#####
if [ ! -f $args_file ]; then 
    $script_dir/create_args.py $args_file $excel $expname
    
else
## CHECK THE SECOND VARIABLE FOR IT TO EXIST
ARGS=`head -2 $args_file | tail -1`

# set variables based on the output of the line ARGS
IFS=: read CHECK <<< $ARGS

fi


# Find number of files 
filenum=`cat $args_file | wc -l`



#### Raw FastQC #####
echo Raw FastQC
base_sbatch=`format_sbatch rawqc` # set prefix
rawqc=$($base_sbatch \
--array=1-${filenum} \
$script_dir/fastqc_batch.sh)



rawqc=${rawqc##* }
echo $rawqc

#### TimeGalore #####
echo Read Trimming
base_sbatch=`format_sbatch trim` # set prefix
trim=$($base_sbatch \
--dependency=afterany:${rawqc} \
--array=1-${filenum} \
$script_dir/trimgalore_batch.sh)


trim=${trim##* }
echo $trim


#### Alginment #####
echo Aligning to hg38 assembly 


#### Sorting and Indexing #####
## This is a main dependency no other process can occur if this doesnt happen successfully. The secret is to query $? which is the exit status of the previous program using the following code repeated for each id
#sacct -j 10576550 --format State | grep -v [S/-] | tr -d ' ' gives FAILED or SUCCESSFUL check for the latter 


## Start first job here requires working directory and the arguments filei
## I will save this into a special command
#outerrsbatch='sbatch --output=./out/'${expname}'.srt.%A_%a.out --error=./err/'${expname}'.srt.%A_%a.err '
#sortid=$($outerrsbatch \
#--array=1-${filenum} \
#--job-name=${expname}_sort \
#sort_index_batch.sh)
#testing_1.sh

# format id
#sortid=${sortid##* }

#echo $sortid

##### Genome Sorting #####


#gen_srt_id=$($outerrsbatch \
#--dependency=afterany:${sortid} \
#--job-name=${expname}_gen \
#create_genometxt.sh )

#--dependency=afterany:${sortid} \

#gen_srt_id=${gen_srt_id##* }



##### Coverage #####

# Batch coverage command goes here 

#cov_id=$($outerrsbatch \
#--dependency=afterany:${gen_srt_id} \
#--array=1-${filenum} \
#--job-name=${expname}_cov \
#coverage_batch.sh)


#covid=$($outerrsbatch --dependency=afterany:${sortid} `\
#testing_2.sh)

#cov_id=${cov_id##* }

# Coverage processing to generate a processed bed file

#pro_cov_id=$($outerrsbatch \
#--dependency=afterany:${cov_id} \
#--job-name=${expname}_p_cov \
#cov_processing.sh)

#pro_cov_id=${pro_cov_id##* } # Make sure that this is for the rendering. 





echo All Job Submitted