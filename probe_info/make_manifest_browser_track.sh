#!/bin/bash

##### PURPOSE #####
## Jacob Gutierrez guteirja@ohsu.edu
## The purpose of this file is to construct a compressed formatted bed track file from the Methyl Capture Probe Validation. 
## To Run: 
## 1) be in the working directory with all the files
## 2) bash make_browser_track.sh

##### UPDATE 7/21/20 #######
This became to complicated to run please see manifest_browsertrack_generation.Rmd
The code found here is incorportated into that file. 

## Cat all the files into one file
awk '{$1="chr"$1 ; print $0}' hg38-ensemble-truseq-manifest.bed > ucsc_format.bed

cat manifest_track_header.txt ucsc_format.bed > hg38-mcseq-truseq-track.bed

gzip hg38-mcseq-truseq-track.bed

## Clear buffer
rm ucsc_format.bed