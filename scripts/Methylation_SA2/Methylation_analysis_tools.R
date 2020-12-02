
### Libraries
require(matrixStats)


## Custom Row Invariance removal function 
## JG 10/30/20: 
remove_invariance <- function(mat, return_idx = F){
  # compute row variance
  all_rows <- matrixStats::rowVars(mat, na.rm = TRUE)
  
  
  ## Find invariant posoitons where == 0
  invariant_idxs <- all_rows == 0   
  
  # Now say some things. 
  nrows <-  nrow(mat)
  to_remove <- sum(invariant_idxs)
  perc <- round(to_remove/nrows, 3) * 100 
  str1 <- sprintf("Total Number of CpGs : %i", nrows)
  str2 <- sprintf("Number of CpGs that are invariant across all samples (removed): %i\nAs percentage of total complete NA's: %f\n",to_remove, perc)
  #str3 <- sprintf("As percentage of total complete NA's: %f", perc)
  str3 <- sprintf("Final Number of CpGs: %i", nrows-to_remove)
  print(str1)
  cat(str2)#print(str2)
  print(str3)
  
  if (return_idx){
    return(invariant_idxs) ## JG 10/30/20: Making the function return the index, 
  } else {
  new_mat <- mat[!invariant_idxs,]
  return(new_mat)
  }
  #return(invariant_idxs) ## JG 10/30/20: Making the function return the index, 
} ## end of remove invariance Function. 



## Custom function to make pretty output for analysis
assess_diffmethobj <- function(obj, total_cpgs, diff = 10, thres = .05){
  
  
  print(sprintf("Subsetting to Signficaint sites with methylation difference greater than %i and adjusted pvalues less than %.1f", diff, thres))
  
  print(summary(obj))
  
  myDiff10p.hyper=getMethylDiff(obj,difference=diff,qvalue=thres,type="hyper")
  myDiff10p.hypo=getMethylDiff(obj,difference=diff,qvalue=thres,type="hypo")
  myDiff10p=getMethylDiff(obj,difference=diff,qvalue=thres)
  
  #total_cpgs <- dim(meth)[1]
  
  total_size <- dim(myDiff10p)[1]
  hyper <- dim(myDiff10p.hyper)[1]
  hypo <- dim(myDiff10p.hypo)[1]
  not_sig <- total_cpgs - total_size
  
  
  display <- data.frame(bav_tav_sig = c(hyper, not_sig, hypo))
  rownames(display) <-  c("hyper","not_sig","hypo")
  print(display)
  
}

## END OF Methylation_analysis_tools.R


jg_methyl_subset_canonical <- function(methylbase){
  ### SUbset based on canonical 
  canonical <- c(seq(22), "X","Y")
  canonical_pos <- methylbase$chr %in% canonical
  
  ## Subset methylbase object
  cleaned_meth <- methylbase[canonical_pos, ]
  return(cleaned_meth)
}

## JG 10/29/20: The methylkit tool was not cutting it I needed a custom main. 
jg_coverageplot <- function(one, title = "bp"){
  tmp <- getData(one)[,c("coverage")]
  plt <- hist(log10(tmp),plot = FALSE)
  my.labs <- as.character(round(100*plt$counts/length(tmp),1))
  
  #mtext(one@sample.id,side = 3)
  
  hist(log10(tmp), plot = TRUE,col = "chartreuse4", 
       xlab=paste("log10 of read coverage per",one@resolution),
       #main=paste("Histogram of", one@context, "coverage"), 
       main = sprintf( "%s Coverage Histogram", title),
       labels = my.labs)
  mtext(one@sample.id,side = 3)
  
  
  ## Report plot to list
  plt <- recordPlot()
  
  return(plt)
}


## JG 10/1/20 UPDATE!
## I will now functionalize this function to accept 2 matrices. First being CpG matrix OR grange object with a 1bp size. Second is a tiled Methylation matrix with various ranges to fill NA value with based on Grange getoverlap function. 

## Start of Imputation Function Is called below. 
## This accepts the output of percMethylation() OR a two matrices with the same n columns that have rowids of chr.start.end format to coherse into granges. 

JG_1kb_imputation <- function(raw, tile_mat){
  
  ## Now I need to query the raw data for a single row and extract the methylation value. 
  imputed <- raw ## Declare Output Space to modify
  buffer <- is.na(imputed) ## Create NA logical matrix
  dims <- dim(imputed) #dims[1 ] = row , dims [2] = Col ## Good reminder and here I am getting the input size.
  ## Converting rownames from . format to bed. 
  imp_rows <- rownames(imputed) %>% stringr::str_split_fixed(.,'\\.',n=3) %>% as.data.frame() %>% mutate_all(as.character)
  colnames(imp_rows) <- c('chr','start','end') ## Unfold to bed
  imp_grange <- makeGRangesFromDataFrame(imp_rows) ## Do the thing! (make grange)
  
  ## Tile Data into grange
  ## Converting rownames from . format to bed. 
  tile_rows <- rownames(tile_mat) %>% stringr::str_split_fixed(.,'\\.',n=3) %>% as.data.frame() %>% mutate_all(as.character) ## Unfold to bed
  colnames(tile_rows) <- c('chr','start','end')  
  rownames(tile_rows) <- rownames(tile_mat)
  final <- cbind(tile_rows,tile_mat) ## Make bed df. 
  fil_grange <- makeGRangesFromDataFrame(final, keep.extra.columns = TRUE) ## Convert to grange with extra data. 
  
  
  ## Get overlaps 
  hits <-  findOverlaps(imp_grange, fil_grange, select = "arbitrary")
  #tile_idx <- subjectHits(hits)
  ## Convert hits into  of methylation value matrix!
  til_dat <- mcols(fil_grange)[hits,] %>% as.matrix()
  
  
  
  #sum(is.na(imputed))
  
  ## Replace NA's with Tile Value. 
  ## For each row
  for (i in seq(dims[1])){
    # Get one row of data
    #print(i)
    a_row <- imputed[i,]
    ## Vector of NA logicals
    to_add <- is.na(a_row)
    
    if(any(to_add)){
      imp_dat <- til_dat[i,]
      a_row[to_add] <- imp_dat[to_add]
      imputed[i,] <- a_row
      #print("IMPUTED!")
    }
    
  }
  
  return (imputed)
} ## End of Tile_mat function. 

## Convert Methylkit Rowids into grange df. 
rows_to_bed <- function(meth_df){
  
  rowsids <- rownames(meth_df)
  
  bed_df <- data.frame(id = rowsids)
  
  col_names <- c("chr","start","end")
  buffer <- bed_df %>% tidyr::separate(id, col_names, "\\.", remove = F) %>% dplyr::select(all_of(col_names),everything())
  
  return(buffer)
  
  
}
