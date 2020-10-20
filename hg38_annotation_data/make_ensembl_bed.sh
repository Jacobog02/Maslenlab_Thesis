#! /bin/bash
## Jacob Gutierrez gutierja@ohsu.edu
## 2/5/20
## Purpose ##
## This script is to reproducibly generate a bed file from the ensemble gtf file.


### Obtain gtf file
## JG 8/21/20: I reproducibly copied the gtf file from hoolock. I will use this as input. 
#wget ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/Homo_sapiens.GRCh38.99.gtf.gz
#gunzip Homo_sapiens.GRCh38.99.gtf.gz # Into plain text

## JG 8/21/20: I reproducibly copied the gtf file from hoolock. I will use this as input.

if [ -f "Homo_sapiens.GRCh38.99.gtf" ]; then
  ## File exists
  echo GTF file previously downloaded
else ##Download and unzip gtf
  wget ftp://ftp.ensembl.org/pub/release-98/gtf/homo_sapiens/Homo_sapiens.GRCh38.98.gtf.gz
  gunzip Homo_sapiens.GRCh38.98.gtf.gz
fi

## Obtain UCSC Scripts
## Note: I modified the code from this link: https://gist.github.com/gireeshkbogu/f478ad8495dca56545746cd391615b93#file-convert_gtf_to_bed12-sh-L8

## Check to see if it exists. 
if [ -f "gtfToGenePred" ]; then
  ## File exists
  echo gtfToGenePred has been downloaded
else ##Download and unzip gtf
  echo Downloading gtftoGenePred
  wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred
  
fi
#wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred

if [ -f "genePredToBed" ]; then
  ## File exists
  echo genePredToBed has been downloaded
else ##Download and unzip gtf
echo Downloading genePredToBed
  wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/genePredToBed
  
fi
#wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/genePredToBed

chmod +x gtfToGenePred genePredToBed # Make executable.


## Convert to bed
./gtfToGenePred Homo_sapiens.GRCh38.98.gtf ensembl_ann.pred
./genePredToBed ensembl_ann.pred ensembl_GRCh38.98_ann.bed

## Cleanup 
rm ensembl_ann.pred
