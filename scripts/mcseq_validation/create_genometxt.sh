#! /bin/bash

#SBATCH -p exacloud             # partition (queue) (exacloud, long_jobs)
#SBATCH -N 1                    # number of nodes
#SBATCH -n 8                    # number of cores
#SBATCH --mem 24000             # memory pool for all cores
#SBATCH -t 0-00:30               # time (D-HH:MM)

echo "SLURM_JOBID: " $SLURM_JOBID



if [ ! -f $wrkdir/srt-${probe_name} ]; then

    $exa_sam view -H $sourcedir/reads/srt.*.bam  | grep @SQ |sed 's/@SQ\tSN:\|LN://g' > $wrkdir/${expname}_genome.txt


    #### Make sure the chr has the string 'chr' in it!!! 

    ## Try #1
    # the -i option operates on the file in place
    #if [ "$add_chr_genome" = true] ; then
    #sed -i -e 's/^/chr/' $wrkdir/${expname}_genome.txt # THIS DIDNT WORK :(
    #fi


    ## Try #2 8/30/19 updated to work with any genome file
    ## 9/4/19 Ok i used the mapping file to remap the chr on the bed file which has reduced representation. 
    #Rscript fix_genome.R $wrkdir/${expname}_genome.txt
      

    # Now check if there is a genome file we would like to merge on. Normally this isnt a deal but if the Coverage breaks try this 

    ## JG 9/3/19 coverage breaks bc the bam header has 1 instead of chr1 it is an issue with the UCSC to ensembl mapping not currently resolved. 

    #if [ -f "$merge_genome" ]; then

    #Rscript merge_genome.R $wrkdir/${expname}_genome.txt $merge_genome

    #fi
      
    $exa_bed sort -g $wrkdir/${expname}_genome.txt -i $probe_bed \
    > $wrkdir/srt-${probe_name}

fi
