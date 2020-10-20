#!/bin/bash
## Jacob Gutierrez gutierja@ohsu.edu
## 2/8/20

## I DID NOT MAKE THIS I SOURCED THIS FROM BRETT.

## JG 8/21/20: Adding command line parameters to improve the workflow. 
## I did add this bit
#wd=/home/groups/hoolock2/u0/jg/methyl-seq-replicate/pipeline_dev
#ann_dir=/home/groups/hoolock2/u0/bd/annot_ref_files/hg38
## JG 8/21/20: Adding command line parameters to improve the workflow. 
input_file=$1
#ann_dir=$2
wd=`dirname $1` 
file_name=`basename $1`

dmr_file=$2

ann_dir=/home/groups/hoolock2/u0/jg/thesis/hg38_annotation_data/
ann_name=ensembl_GRCh38.98_ann.bed
#prom_name=ensembl_GRCh38.98_promoter.bed
t2g_name=ensembl_biomart_t2g.tsv
biomart_name=ensembl_biomart_annotations.tsv

script_dir=/home/groups/hoolock2/u0/jg/thesis/scripts/Methylation_SA2/

## See function for docutmentation of inputs

Rscript ${script_dir}/xlsx_genomation.R \
-i $wd/${file_name} \
-d $dmr_file \
-a $ann_dir/srt.${ann_name} \
-t $ann_dir/${t2g_name} \
-g $wd/bedtools.closest.gene.output \
-p $wd/bedtools.closest.promoters.output \
-b $ann_dir/${biomart_name} \
-n TS_BAV_DMR
