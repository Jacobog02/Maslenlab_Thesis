#!/usr/bin/env Rscript
## Jacob Gutierrez 6/27/19 
## gutierja@ohsu.edu
## This file needs to be modified but it will use arguements to take in a file that has 
##  the paths of bed files to join into one dataframe and write out with another parameter

## JG 8/16/19 Modifying script in order to work in pipeline. 
## Specifically new arguments are : 
## 1) make it read in the file names from the arguments document  
## 2) Working directory to place new file generated 
## 3) Experiment name to name the new file

## OLD::::::Rscript --vanilla draft_cov_merge.R ../ECP4_bed.txt ../ECP4wd/ECP4_multicov.bed

# see the dependency file to see what files these reference but the experiment name is a string! 
#Rscript --vanilla draft_cov_merge.R $args_file $wrkdir $expname

## ------------------------Arguements-------------------------------------
#https://www.r-bloggers.com/passing-arguments-to-an-r-script-from-command-lines/
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)!=3) {
  stop("Parameters needed to process bed files not present. See documentation", call.=FALSE)
} 
#else if (length(args)==1) {
  # default output file
  #args[2] = "multi_out.bed"
#}

# Maybe check file is real
args_file = args[1] 

# Maybe check for valid directory
wrkdir = args[2]

exp_name = args[3]



## ------------------------Libraries-------------------------------------
require(tidyverse)


## ------------------------Get File Paths and Names --------------------------------
all_paths <- read_delim(args_file, delim=':',col_names = c('source_bams','Path'))

#print(all_paths)

# Extracting only the file names
all_paths <- all_paths %>% select(Path)

# Extract unique ID's for appending
#all_paths <-  all_paths %>% mutate(Names = str_split_fixed(all_paths$Path,'/',n=2)[,2]) %>%  mutate(Names = str_split_fixed(Names,'\\.',n=4)[,3])
all_paths <-  all_paths %>% mutate(Names = str_split_fixed(all_paths$Path,'_',n=2)[,2])

# Now formating the names to be the correct file extentions this is just the base name
all_paths <- all_paths %>% mutate(Path = sprintf('cov.%s.bed', Path)) %>% mutate(Path = paste(wrkdir,Path,sep='/'))


# Modify paths to go back once  NOW NOT NEEDED
#all_paths  <- all_paths %>% mutate(Path = paste('..',Path,sep='/')) ## DEBUGGING WONT BE NEEDED


# Get vectors for looping
paths <- all_paths %>% pull(Path)
names <- all_paths %>% pull(Names)

#print(names)

## ------------------------Make Dataframe-------------------------------------

for (i in seq(length(paths))){
  file <-  paths[i]
  id <- names[i]
  cov_header <- c('chr','start','stop', 'probe_id' , paste('read_num',id,sep='_'), paste('base_overlap',id,sep='_'),'probe_length',paste('frac_cov',id,sep='_'))

  ## Check if its #1 then start the main_frame
  if (i == 1){
    main_frame <- read_delim(file,delim='\t',col_names = cov_header,col_types = cols()) %>% select(chr,start,stop,probe_id,probe_length,everything())
  }
  else{

    new_file <-  read_delim(file,delim='\t',col_names = cov_header,col_types = cols())

    main_frame <- left_join(main_frame, new_file, by = c('chr','start','stop','probe_id','probe_length'))
  }

  #dim(main_frame)
}



## ------------------------Write Out DataFrame-------------------------------------
#print(dim(main_frame))
# Formatting outout file! This allows me to standardize the inputs. 
outfile = sprintf('%s_cov.bed', exp_name)
outfile = paste(wrkdir,outfile,sep='/')
write_delim(main_frame,path=outfile)

