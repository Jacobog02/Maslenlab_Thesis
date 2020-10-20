#!/usr/bin/env Rscript
## ----information,eval = F------------------------------------------------
## # Jacob Gutierrez gutierja@ohsu.edu 7/29/19
## # This file needs to be modified but it will use arguements to take in a file that has
## #  the paths of bed files to join into one dataframe and write out with another parameter
##

## JG 8/16/19 Modifying script in order to work in pipeline. 
## Specifically new arguments are : 
## 1) make it read in the file names from the arguments document  
## 2) Working directory to place new file generated 
## 3) Experiment name to name the new file


## OLD::::: # Rscript --vanilla draft_cov_merge.R ../ECP4_bed.txt ../ECP4wd/int_reads.Rdata
# see the dependency file to see what files these reference but the experiment name is a string! 
#Rscript --vanilla int_merge.R $args_file $wrkdir $expname

#wrkdir='/home/groups/hoolock/u1/jg/val_pipe/ECP4_wd'


#args_file='/home/groups/hoolock/u1/jg/probe_val/raw_ECP4_args.txt.fin'
#expname='ECP4'

## ------------------------Arguements-------------------------------------
#https://www.r-bloggers.com/passing-arguments-to-an-r-script-from-command-lines/
args = commandArgs(trailingOnly=TRUE)

## test if there is at least one argument: if not, return an error
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

expname = args[3]


## ----libraries-----------------------------------------------------------
suppressMessages(require(tidyverse))
suppressMessages(require(rtracklayer))

## ----functions-----------------------------------------------------------
format_paths <- function(dat){
  n <- nrow(dat)
  output <- c()
  for (i in seq(n)){
      buffer <- dat[i,]
      buffer <- paste(buffer,collapse = '/')
      output <- c(output,buffer)
  }
  return(output)
}



## ----manual override to debug, eval = F----------------------------------
## args <- c('../ECP4_beds.txt', '../ECP4wd/int_grange.Rdata')


## ----format the input----------------------------------------------------
all_paths <- read_delim(args_file, delim=':',col_names = c('source_bams','Path'))
all_paths <- all_paths %>% select(Path)


# Now modify the file names adding int.
#all_paths <-  all_paths %>% mutate(Names = str_split_fixed(all_paths$Path,'\\.',n=2)[,2]) ## Period delimited
all_paths <-  all_paths %>% mutate(Names = str_split_fixed(all_paths$Path,'_',n=2)[,2]) ##_ delimiteed 


# Now formating the names to be the correct file extentions this is just the base name
all_paths <- all_paths %>% mutate(Path = sprintf('int.%s.bed', Path)) %>% mutate(Path = paste(wrkdir,Path,sep='/'))

#path_buffer <- str_split_fixed(paths, '/' ,n =2) # Format into the matrix
#path_buffer[,2] <- paste('int.',path_buffer[,2], sep = '') # Add the int.
#all_paths$Path <-  format_paths(path_buffer) # Now modify the original dataframe. 


# Extract unique ID's for appending
# THIS NEEDS TO BE GENERALIZED I had to modify the n in Names and the id to extract. 
#all_paths <-  all_paths %>% mutate(Names = str_split_fixed(all_paths$Path,'/',n=2)[,2]) %>%  mutate(Names = str_split_fixed(Names,'\\.',n=5)[,4])

# Modify paths to go back once 
#all_paths  <- all_paths %>% mutate(Path = paste('..',Path,sep='/')) ## DEBUGGING WONT BE NEEDED


# Get vectors for looping
paths <- all_paths %>% pull(Path)
names <- all_paths %>% pull(Names)


#first <- paths[1]

#second <- paste('..',paths[2],sep='/')




## ----generate grangelist-------------------------------------------------
# 9/4/19 Added bc the 
int_cols <- cols(chr = col_character(),
       start = col_double(),
       stop = col_double(),
       id_10 = col_character())


for (i in seq(length(paths))){
  file <-  paths[i] # get file path
  the_id <- names[i] # get id 
  formatted_col <- paste('id',the_id,sep='_')
  cov_header <- c('chr','start','stop', formatted_col )
  
  ## 4/13/20 Adding empty sample comprehension check if no rows if so say that sample is empty. The lack of data will be understood
  #buffer <- read_delim(file,delim=' ', col_names = cov_header,col_types = int_cols) %>% mutate(id = the_id)
  #print(paste('Processing sample:', the_id))
  
  
  ## Check if its #1 then start the main_frame
  if (i == 1){
   #main_frame <- makeGRangesFromDataFrame(read_delim(file,delim=' ', col_names = cov_header,col_types = cols(), skip =1) %>% mutate(id = the_id))
    buffer <- read_delim(file,delim=' ', col_names = cov_header,col_types = int_cols) %>% mutate(id = the_id) %>% select(-formatted_col)
    main_grange <- makeGRangesListFromDataFrame( buffer,  split.field = 'id' , ignore.strand = T, keep.extra.columns = T)
    main_grange[[1]]$id = the_id
    #main_grange <- GRangesList(buffer) # buffer was just a Grange object NOT a list
  }else{
    #print(file)
    
    buffer <- read_delim(file,delim=' ', col_names = cov_header,col_types = int_cols) %>% mutate(id = the_id) %>% select(-formatted_col)
    buffer <- makeGRangesListFromDataFrame( buffer,  split.field = 'id' , ignore.strand = T, keep.extra.columns = T)
    buffer[[1]]$id = the_id
    main_grange <- suppressWarnings(append(main_grange, buffer))
    #length(main_grange)
    #main_frame <- left_join(main_frame, new_file, by = c('chr','start','stop'))
  }
  
  #print(dim(main_frame))
}



## ------------------------------------------------------------------------
buffer <- main_grange %>% unlist() #%>% as.data.frame() %>% head() # This commented stuff doesnt work 
names(buffer) <-NULL # need this to coherse into dataframe
# Issue with function conflict must explicityly say dplyr::
out_of_probe_reads <- buffer %>% as.data.frame() %>% dplyr::rename('chr' = seqnames) %>% select(-strand) # make it long! 8/20/19 renamed to chr and remove strand



## ----save output---------------------------------------------------------
# NEED THE OUTFILE MAYBE MAYBE IT AN ARGUMENT
#getwd()
outfile = sprintf('%s_int.Rdata', expname)
outfile = paste(wrkdir,outfile,sep='/')
save(main_grange, out_of_probe_reads, file = outfile)

