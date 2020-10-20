#!/bin/bash
## JG 8/16/20: This script can be called as apart of the pipeline or directly by the user. 
## Checks to see if it is running on a node then checks pipeline otherwise it just removes excess large files. 

#SBATCH -p exacloud                # partition (queue) (exacloud, long_jobs)
#SBATCH -N 1                       # number of nodes
#SBATCH -n 4                       # number of cores
#SBATCH --mem 2400              # 36 GB memory pool for all cores
#SBATCH -t 00-01:00                 # time (D-HH:MM)


## JG 8/16/20: Check if on slurm node id is set.
if [ ! -z "$SLURM_JOBID"]; then 

  echo "SLURM_JOBID: " $SLURM_JOBID

  ##### Checking Dependency Success #####
  dependid=${SLURM_JOB_DEPENDENCY##*:}

  deptasks=(`sacct -j $dependid --format State | grep -v [S/-] | tr -d ' ' `)

  check_success ${deptasks[@]}

fi
## This script deletes extra files that take up space and are not used downstream. 


## JG 8/18/20: Too many files I found some code to do this here: 
## https://stackoverflow.com/questions/14765569/test-if-multiple-files-exist

####### JG 8/17/20: Added Disc Usage Before and after. 
echo "####### Current Disk Usage #######"

du $wrkdir -h

echo "###################################"

## set outdir path
#$wrkdir outdir=/home/groups/hoolock/u1/jg/replicate/qc/output

echo Now Removing Excess Files! 

# Remove large files from 1) trimming and 2) raw alignments
## Trimming: trimmed reads are large
# catch=`ls -1 $wrkdir/trim/*.fq.gz`
# 
# if [ -f $catch ]; then 
# echo Removing Trimmed Data
# #rm $wrkdir/trim/*.fq.gz
# 
# else
# echo No Trim Bam Files Detected: Skipping
# fi

### OLD STUFF

if [ `ls -1 $wrkdir/trim/*.fq.gz 2>/dev/null | wc -l ` -gt 0 ];
then
    echo "Files Detected Now Removing Trimmed"
    rm $wrkdir/trim/*.fq.gz
else
    echo "No Trimming Files Detected"
fi



## Bismark Alignment: unsorted bam files
#rm $wrkdir/align/*.bam
#catch=`ls ${wrkdir}align/*.bam`

if [ `ls -1 ${wrkdir}align/*.bam 2>/dev/null | wc -l ` -gt 0 ];
then
    echo "Files Detected Now Removing Aligned"
    rm $wrkdir/align/*.bam

else
echo No Alig Bam Files Detected: Skipping
fi


## Bismark coverage: 
#rm $wrkdir/meth/C*.txt.gz
#catch=`ls $wrkdir/meth/C*.txt.gz`

#if [ -f $catch ]; then 
if [ `ls -1 $wrkdir/meth/C*.txt.gz 2>/dev/null | wc -l ` -gt 0 ];
then 
  echo Removing CpG Coverage Data
  rm $wrkdir/meth/C*.txt.gz
  
  else
  echo No coverage Files Detected: Skipping
fi
## Bismark nucleotide report: 


#### Cleaned Directory Disc Usage
echo "####### Current Disk Usage #######"

du $wrkdir -h

echo "###################################"




