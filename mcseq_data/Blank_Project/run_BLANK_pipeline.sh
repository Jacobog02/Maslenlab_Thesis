#!/bin/bash

##### PURPOSE #####
## Jacob Gutierrez guteirja@ohsu.edu
## The purpose of this file is to intitialize the probe validation program.
## Basically You declare the required parameters within this script and the named variables are then inherited to the enviroment to be used within the cluster. 

## Requires the Experiment setup portion to be specified including:
# Excel mapping file to generate the project specific mappings to be deployed on the cluster. 

# Run by bash run_EXPERIMENTNAME.sh which sets up the variables and calls sbatch sending the job to the cluster. 

#### Experiment Set up #####
# Excel mapping file for all experiments
## JG 10/11/20: Consider making this an argument but it does make more sense to manually declare paths with a wrapper. 
## Exacloud is easier to work with when accessing files using explicit paths from /home/groups/....
#export excel=/home/groups/hoolock2/u0/jg/thesis/mcseq_data/mc-seq_mapping.xlsx  #Excel mapping file.
export excel=/home/groups/hoolock2/u0/jg/thesis/mcseq_data/mc-seq_mapping_plusECP42.xlsx  #Excel mapping file.
# common name for experiment
export expname=ECP4 # String for the experiment name to be used for standard output notation no spaces andbest if its short. 
# Where the output be placed 
export wrkdir=/home/groups/hoolock2/u0/jg/thesis/mcseq_data/${expname}/ # Absolute path to the directory to be used for output writing 

# Formatted Arguements file can be used for all processes just make it focus on new name 
## This is generated from the excel file!
export args_file=${wrkdir}/${expname}_mapping # Absolute path to argument files (colon delimited)

# Where the source raw fastq files live
#export sourcedir=/home/groups/hoolock2/u0/archive/Epigenetics_Core_rawdata/ECP4 # Absolute path to the directory where the original fastq files are contained
#export sourcedir=/home/groups/hoolock2/u0/archive/Epigenetics_Core_rawdata/ECP10_clean/methylseq
#export sourcedir=/home/groups/hoolock2/u0/archive/Epigenetics_Core_rawdata/ECP15
#export sourcedir=/home/groups/hoolock2/u0/archive/Epigenetics_Core_rawdata/ECP16/re_demultiplex
export sourcedir=/home/groups/hoolock2/u0/archive/Epigenetics_Core_rawdata/ECP42 # Absolute path to the directory where the original fastq files are contained

## Where the neccesary scripts to run the pipeline live
export script_dir=/home/groups/hoolock2/u0/jg/thesis/scripts/mcseq_alignment # Absolute path to directory for the pipeline scripts

## Set Slurm Account
export slurm_acct=maslenlab


##### Check or Create all Needed Directoires
mkdir -p slurmout slurmerr qc trim align reads meth

##### Launch the Slurm Process! #####

## Alignment and Methylation Calling
sbatch --account $slurm_acct --job-name=${expname}_dependency $script_dir/mcseq_SLURM_pipeline.sh


## Cleanup Big Files. 
#bash $script_dir/remove_excess_reads.sh

