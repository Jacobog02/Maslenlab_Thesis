#!/home/groups/hoolock2/u0/jg/thesis/env/methylseq/bin/python

## Jacob Gutierrez gutierja@ohsu.edu
### 3/29/20

## This script takes in two arguements 1) Working Directory1) Excel file 2) Tab name from the command line


### libraries
import sys # Catch arguments
import pandas as pd # process excel
#import numpy as np

## Catch!
args_file = sys.argv[1]
#print("HELLO IN PYTHON")
#print("Args Type:" + type(args_file))
## Test to see if this is empty
#if args_file is None: stop("Arguments File Missing")
mappings = sys.argv[2]
project_name = sys.argv[3]


## Read df
df = pd.read_excel(mappings,project_name)

## Write colon delimited file
#df[["File_Prefix","Sample_Name"]].to_csv(args_file, sep = ":", header = False, index=False)

## JG 10/10/20: Adding Paired/Single End processes with same script.
#df[["File_Prefix","Sample_Name","Seq_Type"]].to_csv(args_file, sep = ":", header = False, index=False)

## JG 12/15/20: Modified to include Lib_Type for deduplication
df[["File_Prefix","Sample_Name","Seq_Type","Lib_Type"]].to_csv(args_file, sep = ":", header = False, index=False)





