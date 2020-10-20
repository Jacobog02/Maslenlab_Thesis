#!/bin/bash
## Jacob Gutierrez gutierja@ohsu.edu
## 2/6/2020 
## THIS IS NOT MY SCRIPT THANK YOU BRETT!
## Copied from /home/groups/hoolock2/u0/bd/Projects/ECP15/annotate/DMC/GR1case_v_GR1ctrl/bedclosest.sh
## 2/10/20 I got reference annotation files from brett in /home/groups/hoolock2/u0/bd/annot_ref_files/hg38/

## ON HOOLOCK WE MUST FIRST LOAD THE BEDTOOLS MODULE IN SHELL BEFORE RUNNING POTENTIALLY
ml BEDTools


## JG 8/21/20: Adding command line parameters to improve the workflow. 
input_file=$1
#ann_dir=$2
outdir=`dirname $1` 
file_name=`basename $1`

Ann_dir=/home/groups/hoolock2/u0/jg/thesis/hg38_annotation_data/
Ann_name=ensembl_GRCh38.98_ann.bed
Prom_name=ensembl_GRCh38.98_promoter.bed


## NOTE!! All files must be bed tools sorted to work!

## Gene Annotation Sorting Check
if [ -f ${Ann_dir}/srt.${Ann_name} ]; then
  ## Already exists
  echo "Gene Annotation already sorted"
else
  ## Doesnt exist make it.
  echo "Sorting Gene Annotations"
  # #sort -k1,1n -k2,2n  hg38_annotation_data/ensembl_GRCh38.99_ann.bed >  hg38_annotation_data/ensembl_GRCh38.99_ann.srt.bed
  #bedtools sort -i hg38_annotation_data/ensembl_GRCh38.98_ann.bed > hg38_annotation_data/ensembl_GRCh38.99_ann.srt.bed
  bedtools sort -i ${Ann_dir}/${Ann_name} > ${Ann_dir}/srt.${Ann_name}

  
fi

## Promoter Annotation Sorting Check
if [ -f ${Ann_dir}/srt.${Prom_name} ]; then
  ## Already exists
  echo "Promoter Annotation already sorted"
else
  ## Doesnt exist make it.
  echo "Sorting Promoter Annotations"
  # #sort -k1,1n -k2,2n  hg38_annotation_data/ensembl_GRCh38.99_promoter.bed >  hg38_annotation_data/ensembl_GRCh38.99_promoter.srt.bed
  #bedtools sort -i hg38_annotation_data/ensembl_GRCh38.98_promoter.bed >  hg38_annotation_data/ensembl_GRCh38.99_promoter.srt.bed
  bedtools sort -i ${Ann_dir}/${Prom_name} >  ${Ann_dir}/srt.${Prom_name}

  
fi


#sort -k1,1n -k2,2n  data/sig_methdiff_DMR_results.bed >  data/sig_methdiff_DMR_results.srt.bed
#bedtools sort -i data/sig_methdiff_DMR_results.bed >  data/sig_methdiff_DMR_results.srt.bed
## Input File Sorti
## 8/21/20 Accept parameters 
if [ -f ${outdir}/srt.${file_name} ]; then
  ## Already exists
  echo "Output already sorted"
else
  ## Doesnt exist make it.
  echo "Sorting Input File"
  bedtools sort -i ${input_file} > ${outdir}/srt.${file_name}
fi

# gene
## NOTE: input bed has chr, start, end, id
## NOTE NOTE:  chr is just an integer doesnt have chr in front of it....
if [ -f ${outdir}/"bedtools.closest.gene.output" ]; then
  ## File exists
  echo "Closest Genes Previously Done"
else
  echo "Finding Closest Gene Regions"
  bedtools closest -a ${outdir}/srt.${file_name} \
  -b ${Ann_dir}/srt.${Ann_name} \
  -d > ${outdir}/bedtools.closest.gene.output
fi

# promoters

if [ -f ${outdir}/"bedtools.closest.promoters.output" ]; then
  ## File exists
  echo "Closest Promoters Previously Done"
else
  echo "Finding Closest Promoter Regions"
  bedtools closest -a ${outdir}/srt.${file_name} \
  -b ${Ann_dir}/srt.${Prom_name} \
  -d > ${outdir}/bedtools.closest.promoters.output
fi


echo "Complete!"
## IT RAN :DDDD 2/7/20
