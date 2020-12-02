#!/bin/bash

##### PURPOSE #####
## Jacob Gutierrez guteirja@ohsu.edu
## The purpose of this file is to intitialize the probe validation program.
## Basically You declare the required parameters within this script and the named variables are then inherited to the enviroment to be used within the cluster. 

## Requires the Experiment setup portion to be specified including:
# formatted file of of the source bams names colon delimited with the and the new file names to simplify manually formated. 

# Run by bash run_EXPERIMENTNAME.sh which sets up the variables and calls sbatch sending the job to the cluster. 

#### Experiment Set up #####

# common name for experiment
export expname=ECP4

# Where the source bam live
export sourcedir=/home/groups/hoolock2/u0/jg/thesis/mcseq_data/${expname}/

# Where the bams will be converted to bed output
export wrkdir=/home/groups/hoolock2/u0/jg/thesis/sa1/${expname}/

# Formatted Arguements file can be used for all processes just make it focus on new name
export args_file=/home/groups/hoolock2/u0/jg/thesis/mcseq_data/${expname}/${expname}_mapping

# Probe file 
export probe_bed=/home/groups/hoolock2/u0/jg/thesis/probe_info/hg38-ensemble-truseq-manifest.bed

## Where the neccesary scripts to run the pipeline live
export script_dir=/home/groups/hoolock2/u0/jg/thesis/scripts/mcseq_validation # Absolute path to directory for the pipeline scripts

## Set Slurm Account
export slurm_acct=maslenlab

##### Check or Create all Needed Directoires
mkdir -p slurmout slurmerr

sbatch --account $slurm_acct --job-name=${expname}_val_dependency $script_dir/mcseq_probe_pipeline.sh

