cp ma #!/bin/bash
## Jacob Gutierrez gutierja@ohsu.edu
## 2/6/2020 
## THIS IS NOT MY SCRIPT THANK YOU BRETT!
## Copied from /home/groups/hoolock2/u0/bd/Projects/ECP15/annotate/DMC/GR1case_v_GR1ctrl/bedclosest.sh

## ON HOOLOCK WE MUST FIRST LOAD THE BEDTOOLS MODULE IN SHELL BEFORE RUNNING POTENTIALLY
ml BEDTools


## NOTE!! All files must be bed tools sorted to work!

## USE UNIX SORT INSTEAD!!!! MAYBE :/ MAYEB DONT.... I dont know how to use the sort utility I might stick to bedtools. 
#bedtools sort -i ensembl_GRCh38.99_ann.bed > ensembl_GRCh38.99_ann.srt.bed
#bedtools sort -i data/sig_methdiff_DMR_results.bed  > data/sig_methdiff_DMR_results.srt.bed ## THIS FAILED BC MY PROMOTER BED FILE WAS BROKEN!!!

echo "Sorting Gene Annotations"
#sort -k1,1n -k2,2n  hg38_annotation_data/ensembl_GRCh38.99_ann.bed >  hg38_annotation_data/ensembl_GRCh38.99_ann.srt.bed
bedtools sort -i hg38_annotation_data/ensembl_GRCh38.99_ann.bed > hg38_annotation_data/ensembl_GRCh38.99_ann.srt.bed


echo "Sorting Promoter Annotations"
#sort -k1,1n -k2,2n  hg38_annotation_data/ensembl_GRCh38.99_promoter.bed >  hg38_annotation_data/ensembl_GRCh38.99_promoter.srt.bed
bedtools sort -i hg38_annotation_data/ensembl_GRCh38.99_promoter.bed >  hg38_annotation_data/ensembl_GRCh38.99_promoter.srt.bed

echo "Sorting Input File"
#sort -k1,1n -k2,2n  data/sig_methdiff_DMR_results.bed >  data/sig_methdiff_DMR_results.srt.bed
bedtools sort -i data/sig_methdiff_DMR_results.bed >  data/sig_methdiff_DMR_results.srt.bed


# gene
## NOTE: input bed has chr, start, end, id
## NOTE NOTE:  chr is just an integer doesnt have chr in front of it....
bedtools closest -a /home/groups/hoolock2/u0/jg/methyl-seq-replicate/pipeline_dev/data/sig_methdiff_DMR_results.srt.bed \
-b /home/groups/hoolock2/u0/jg/methyl-seq-replicate/pipeline_dev/hg38_annotation_data/ensembl_GRCh38.99_ann.srt.bed \
-d > /home/groups/hoolock2/u0/jg/methyl-seq-replicate/pipeline_dev/data/bedtools.closest.gene.output

# promoters
bedtools closest -a /home/groups/hoolock2/u0/jg/methyl-seq-replicate/pipeline_dev/data/sig_methdiff_DMR_results.srt.bed \
-b /home/groups/hoolock2/u0/jg/methyl-seq-replicate/pipeline_dev/hg38_annotation_data/ensembl_GRCh38.99_promoter.srt.bed \
-d > /home/groups/hoolock2/u0/jg/methyl-seq-replicate/pipeline_dev/data/bedtools.closest.promoters.output

## IT RAN :DDDD 2/7/20
