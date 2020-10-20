#!/bin/bash

##### PURPOSE #####
## Jacob Gutierrez guteirja@ohsu.edu
## The purpose of this file is to intitialize the probe validation program. Additionally it will ensure all process are run with the correct dependencies and finally generating HTML output. 

# Could consider manually setting the probe_bed file

## Requires the Experiment setup portion to be specified (perhaps consider making command line 

## AND! a formatted file of of the source bams and the new file names to simplify. Manually formated. 

# Run by sbatch probe_pipeline.sh in a wrapper script that includes the experiment specific information.  

##### SBATCH INITALIZATION #####
#SBATCH -p exacloud                # partition (queue) (exacloud, long_jobs)
#SBATCH -N 1                       # number of nodes
#SBATCH -n 2                      # number of cores
#SBATCH --mem 2400                # memory pool for all cores
#SBATCH -t 0-00:30                 # time (D-HH:MM)
#SBATCH -o ./slurmout/source.%A_%a.out           # Standard output
#SBATCH -e ./slurmerr/source.%A_%a.err           # Standard error


#### Experiment Set up #####
echo Depedency Script

# common name for experiment
#expname=ECP4

# Where the source bam live
#sourcedir=/home/groups/hoolock/u1/bd/Projects/ECP4/bismark/output3_min20max75

# Where the bams will be converted to bed output
#wrkdir=/home/groups/hoolock/u1/jg/val_pipe/ECP4_wd

# Formatted Arguements file can be used for all processes just make it focus on new name
#args_file=/home/groups/hoolock/u1/jg/probe_val/raw_ECP4_args.txt.fin
#args_file=/home/groups/hoolock/u1/jg/val_pipe/sleep_args

# Find number of files 
filenum=`cat $args_file | wc -l`
#filenum=3

# Probe file 
#probe_bed=/home/groups/hoolock/u1/jg/val_pipe/new-srt-hg38-truseq.bed
#genome=/home/groups/hoolock/u1/jg/val_pipe/new_genome.txt

export probe_name=`basename $probe_bed`

##### Software and Functions #####
#exa_sam=/opt/installed/samtools-1.6/bin/samtools

## Environment specific samtools
export exa_sam=/home/groups/hoolock2/u0/jg/thesis/env/methylseq/bin/samtools
## Enviroment specific bedtools
export exa_bed=/home/groups/hoolock2/u0/jg/thesis/env/methylseq/bin/bedtools

check_success () {
echo inside check success

#echo "$@"

for task in "$@"; do
    #echo $task 
    if [ $task != "COMPLETED" ]; then
    #echo stopping script here
    echo Failure Detected Stopping Job!
    exit 1 # fail here? Or MAYBE try scancel 
    fi
done 
}

export -f check_success

## Build sbatch command 
format_sbatch() {

#outerrsbatch='sbatch --output=./out/'${expname}'.'$1'.%A_%a.out --error=./err/'${expname}'.'$1'.%A_%a.err '
  echo 'sbatch --account '$slurm_acct' --job-name='${expname}'_'$1' --output=./slurmout/'${expname}'.'$1'.%A_%a.out --error=./slurmerr/'${expname}'.'$1'.%A_%a.err'
  
}


##### Export Commands #####
## Export sends the variable into the enviroment to be inherited by all child processes. Like turtles all the way down. 

## Commented commands are done in the wrapper
#export expname
#export probe_name
#export exa_sam
#export sourcedir
#export wrkdir
#export args_file
#export filenum
#export probe_bed
#export probe_name
#export exa_sam
#export -f check_success

#### Sorting and Indexing #####
## This is a main dependency no other process can occur if this doesnt happen successfully. The secret is to query $? which is the exit status of the previous program using the following code repeated for each id
#sacct -j 10576550 --format State | grep -v [S/-] | tr -d ' ' gives FAILED or SUCCESSFUL check for the latter 



##### Genome Sorting #####
echo Genome Sorting
base_sbatch=`format_sbatch gensrt` # set prefix
gen_srt_id=$($base_sbatch \
$script_dir/create_genometxt.sh )

#--dependency=afterany:${sortid} \

gen_srt_id=${gen_srt_id##* }

# exit 1 ## Debugged 4/11/20 Couldnt get samtools to load see (https://github.com/sunbeam-labs/sunbeam/issues/181) 


##### Coverage #####

# Batch coverage command goes here 
base_sbatch=`format_sbatch cov`
cov_id=$($base_sbatch \
--dependency=afterany:${gen_srt_id} \
--array=1-${filenum} \
$script_dir/coverage_batch.sh)


#covid=$($outerrsbatch --dependency=afterany:${sortid} `\
#testing_2.sh)

cov_id=${cov_id##* }

# Coverage processing to generate a processed bed file
base_sbatch=`format_sbatch p_cov`
pro_cov_id=$($base_sbatch \
--dependency=afterany:${cov_id} \
$script_dir/cov_processing.sh)

pro_cov_id=${pro_cov_id##* } # Make sure that this is for the rendering. 






##### Interval ######

# Batch interval command goes here
base_sbatch=`format_sbatch int`
int_id=$($base_sbatch \
--dependency=afterany:${gen_srt_id} \
--array=1-${filenum} \
$script_dir/intersect_batch.sh)

int_id=${int_id##* }

# Interval processing here 1) clean the files making int_out.txt 2) Make granges object
base_sbatch=`format_sbatch p_int`
pro_int_id=$($base_sbatch \
--dependency=afterany:${int_id} \
$script_dir/int_processing.sh)

pro_int_id=${pro_int_id##* } # Make sure that this is for the rendering. 




###### Visualizations #####
## MAKE FINAL OUTPUT DATA START WITH EXPERIMENT NAME!!! THIS LETS THE MARKDOWN FIND EVERYTHING STANDARD.

# Parameters experiment name, working dir to get data files.  
# Render command goes here with experiment name and working dir 
base_sbatch=`format_sbatch vis`
vis_id=$($base_sbatch \
--dependency=afterany:${pro_int_id}:${pro_cov_id} \
$script_dir/vis_processing.sh)




