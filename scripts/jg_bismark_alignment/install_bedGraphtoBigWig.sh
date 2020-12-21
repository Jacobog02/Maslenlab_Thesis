#! /bin/bash
## Jacob Gutierrez gutierja@ohsu.edu
## 2/5/20
## Purpose ##
## This script is to reproducibly generate a bed file from the ensemble gtf file.

## JG 10/30/20: I am usurping this script to download the bigwig converter. 
## I need this to rapidly visualize methyaltion data. 
if [ -f "bedGraphToBigWig" ]; then
  ## File exists
  echo bedGraphToBigWig has been downloaded
else ##Download and unzip gtf
  echo Downloading bedGraphToBigWig
  wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig
  
fi

if [ -f "fetchChromSizes" ]; then
  ## File exists
  echo fetchChromSizes has been downloaded
else ##Download and unzip gtf
  echo Downloading fetchChromSizes
  wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/fetchChromSizes
  
fi

### Obtain gtf file
## JG 8/21/20: I reproducibly copied the gtf file from hoolock. I will use this as input. 
#wget ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/Homo_sapiens.GRCh38.99.gtf.gz
#gunzip Homo_sapiens.GRCh38.99.gtf.gz # Into plain text

## JG 8/21/20: I reproducibly copied the gtf file from hoolock. I will use this as input.
# 
# if [ -f "Homo_sapiens.GRCh38.99.gtf" ]; then
#   ## File exists
#   echo GTF file previously downloaded
# else ##Download and unzip gtf
#   wget ftp://ftp.ensembl.org/pub/release-98/gtf/homo_sapiens/Homo_sapiens.GRCh38.98.gtf.gz
#   gunzip Homo_sapiens.GRCh38.98.gtf.gz
# fi

## Obtain UCSC Scripts
## Note: I modified the code from this link: https://gist.github.com/gireeshkbogu/f478ad8495dca56545746cd391615b93#file-convert_gtf_to_bed12-sh-L8

## Check to see if it exists. 
# if [ -f "gtfToGenePred" ]; then
#   ## File exists
#   echo gtfToGenePred has been downloaded
# else ##Download and unzip gtf
#   echo Downloading gtftoGenePred
#   wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred
#   
# fi
#wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred

# if [ -f "genePredToBed" ]; then
#   ## File exists
#   echo genePredToBed has been downloaded
# else ##Download and unzip gtf
# echo Downloading genePredToBed
#   wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/genePredToBed
#   
# fi
#wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/genePredToBed

#chmod +x gtfToGenePred genePredToBed # Make executable.

chmod +x bedGraphToBigWig fetchChromSizes # Make executable.


## Parallel Examples! Considering doing the man tutorial
#https://gist.github.com/Brainiarc7/7af2ab5e88ef238da2d9f36b4be203c0

#https://github.com/FelixKrueger/Bismark/issues/274

export bedGraphToBigWig=/home/groups/hoolock2/u0/jg/thesis/scripts/mcseq_alignment/bedGraphToBigWig
#export fetchChromSizes=/home/groups/hoolock2/u0/jg/thesis/scripts/mcseq_alignment/fetchChromSizes
#export ChromSizes=/home/groups/hoolock2/u0/jg/thesis/annotation_info/hg38.chrom.size
export ChromSizes=/home/groups/hoolock2/u0/jg/thesis/sa1/ECP15/ECP15_genome.txt

## HERE is how I parallelized cleaning .bedgraph bismark output to be bigwigs. 
parallel gunzip {} ::: *.bedGraph.gz ## Unzip the files. I cant figureout how to pass it into the ensembl function

parallel sed -i '1d' {} ::: *.bedGraph ## Clean the Title track, apperatly the function doesnt like it. 

parallel 'bedtools sort -i {} > {.}.srt' ::: *.bedGraph ## Sort the data. 

## I need to use the genome.txt file! Consider when I should do this. I placed it in the SLURM pipeline. 
## The chrome sizes are computed from the bam data for SA1.... 
parallel $bedGraphToBigWig {} $ChromSizes {.}.bw ::: *.srt



# ## Convert to bed
# ./gtfToGenePred Homo_sapiens.GRCh38.98.gtf ensembl_ann.pred
# ./genePredToBed ensembl_ann.pred ensembl_GRCh38.98_ann.bed
# 
# ## Cleanup 
# rm ensembl_ann.pred
