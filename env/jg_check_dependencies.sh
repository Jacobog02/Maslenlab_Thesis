#!/bin/bash

## Jacob Gutierrez
## 4/15/20
## gutierja@ohsu.edu


##### PURPOSE #####

#The purpose of this script is to make sure all dependencies are callable.


##########

## Check samtools Known to sometimes be buggy
samtools --help


## Check python: pandas 
python check_pandas.py
#import pandas as pd  
#pd.read_excel("../mcseq_data/mc-seq_mapping.xlsx") 
#quit()



## Check R: rtracklayer, tidyverse
Rscript check_R.R
#require(rtracklayer)
#require(tidyverse)
#q("no")



## End of dependency check
echo End of Dependency check


