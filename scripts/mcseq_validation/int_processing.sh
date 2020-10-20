#!/bin/bash

#SBATCH -p exacloud             # partition (queue) (exacloud, long_jobs)
#SBATCH -N 1                    # number of nodes
#SBATCH -n 8                    # number of cores
#SBATCH --mem 24000             # memory pool for all cores
#SBATCH -t 0-00:30               # time (D-HH:MM)

echo "SLURM_JOBID: " $SLURM_JOBID
#echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
#echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID

## This script is dependent on the coverage script being run sucessfully and will process the output of each file into one merged file named $expname_cov.bed. 
##### Checking Dependency Success #####
dependid=${SLURM_JOB_DEPENDENCY##*:}
deptasks=(`sacct -j $dependid --format State | grep -v [S/-] | tr -d ' ' `)

check_success ${deptasks[@]}

# GO!
output=$wrkdir/${expname}_int_output.csv 
buffer=buffer.txt
files=` ls $wrkdir/int.* `
#echo $files

# Only run this if the output doesnt exist yet! 
## 9/4/19 BUG: If the output exists in the pipeline it will skip this and a parsing error in int_merge.R will happen! 

if [ ! -f "$output" ]; then
   # Add column titles
echo  'filename,total reads,reads out!' >  $buffer

# Add info file name , read total, read outof frame. 

for file in $files;  do 

  echo `echo $file | awk '{n=split($0,a,"/"); print a[n];}'` ',' >> $buffer ; 
  
  echo `head -1 $file` ',' >> $buffer ;
  
  echo `cat $file | wc -l`  '!' >> $buffer ;
  
  # Remove first line of data here so each file is in bed format. 
  sed '1d' $file > tmpfile; mv tmpfile $file

done

#buffer=`cat $output | tr -d '\n' | tr -d ' ' | tr '!' '\n'` 
## Modify the buffer to be the correct output
cat $buffer | tr -d '\n' | tr -d ' ' | tr '!' '\n' > $output
rm $buffer # goodbye buffer!
fi ## End of output csv

if [ ! -f $wrkdir/${expname}_int.Rdata ]; then

    Rscript --vanilla $script_dir/int_merge.R $args_file $wrkdir $expname
fi
