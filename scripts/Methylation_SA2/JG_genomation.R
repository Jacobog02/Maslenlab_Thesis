#!/usr/bin/Rscript

# NOTE: For 'weird' genomes with scaffolds and no chromosomes (killifish), made the following changes:
# use 'remove.unusual=FALSE' on line 112 when reading in the bed12 file 

# SANITY CHECKS to add:
# make sure there are no E's in promoter reference bed file used for bedtools closest?
# make sure all end - start for promoter regions are <= 3000
# columns 6 and 7 from bedtools.closest.promoter.output

# EDITS
# do you really want to remove duplicate unqIDs, multiple overlaps, etc
# answer: at the end, we want to remove rows with duplicate:
#  <unqID> <gene.overlap.GeneID> <promoter.GeneID> <tss.GeneID>
#  meaning, if all 4 of these values are the same, remove dups and keep 1 row

# SCRIPT IMPROVEMENTS:
# modularize?
# run bedtools closest from within this script?
#  if you do this, have an option to create the promoter bed file from the bed12 file, or supply the promoter bed file
# clean up in general


## JG 8/24/20: Updated script to directly query description data from biomart. 
## Issue is annotation data downloaded from biomart has source (ensembl, havana, mirbase, etc) limitations can only get one at a time. 
## Will also attempt to improve workflow and automate certain workflows such as computing the closest gene and closest promoters. 
## I will now remove the bedclosest step and just use genomation to get geneIDs of overlaping regions. 


#-------------------------------------------------------------------------------------------------------#
# HOW TO RUN:												#
#													#
# Command Line Arguments										#
#													#
# 1. -i, --infile: input regions BED file in format <chr><start><end><region_id>				#
#                   must contain canonical contigs only, no header					#
# 2. -a, --annoBED: path to ensembl annotation BED12 file   						#
# 3. -t, --t2g: path to t2g reference file: <GeneID> <TranscriptID>, no header				#
# 5. -g, --bcog: bedtools closest gene output file   		     					#
# 6. -p, --bcop: bedtools closest promoter output file							#
# 7. -b, --mart: path to bioMart gene info ref file:							#
#        <GeneID><gene.name><source.of.gene.name><gene.type><gene.description><Ensembl.family.desc>	#
# 8. -n, --name: name at the beginning of the output table						#
# 9. -x, --xlsx: Output the results in XLSX format #
#-------------------------------------------------------------------------------------------------------#

require("optparse")

#----------------------#
# COMMAND LINE OPTIONS #
#----------------------#
option_list = list(
  make_option(c("-i", "--infile"), type="character", default=NULL, 
              help="path to input regions bed file: <chr> <start> <end> <region_id>", metavar="character"),
  make_option(c("-d", "--dmr"), type="character", default=TRUE,
              help="Path to DMR results.txt file", metavar="character"),
  make_option(c("-a", "--annoBED"), type="character", default=NULL, 
	      help="path to ensembl annotation BED12 file", metavar="character"),
  make_option(c("-t", "--t2g"), type="character", default=NULL,
  	      help="path to t2g reference file: <GeneID> <TranscriptID>, no header", metavar="character"),
  make_option(c("-n", "--name"), type="character", default=NULL,
              help="name at the beginning of the output table", metavar="character") 
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


#-----------------------------------#
# COMMAND LINE OPTION SANITY CHECKS #
#-----------------------------------#
if (is.null(opt$infile)){
  print_help(opt_parser)
  stop("Missing input file.n", call.=FALSE)
}
if (is.null(opt$annoBED)){
  print_help(opt_parser)
  stop("Missing annotation BED file.n", call.=FALSE)
}
if (is.null(opt$t2g)){
  print_help(opt_parser)
  stop("Missing TranscriptID to GeneID file.n", call.=FALSE)
}
if (is.null(opt$name)){
  print_help(opt_parser)
  stop("Missing output prefix name argument.n", call.=FALSE)
}

## Debuging
# saveRDS(opt, 'parser.rds')
# stop("Parser saved for debugging")
# opt <- readRDS('parser.rds')

# BEGINNING OF SCRIPT

require(genomation)
#require(GenomicRanges) ## rtracklayer is used 
require(rtracklayer)
require(dplyr)
require(writexl)
require(biomaRt)
select <- dplyr::select


################################ Read in Bed file to annotate ################################ 

# read in bed file of input regions ( -i )
## USE THIS LATER 
myreg <- read.table(opt$infile, header=FALSE)
## JG 8/23/20: Changed the uniID to region_id to match the new differential methylation output
names(myreg) <- c("chr","start","end","region_id")


# convert to GRanges object
#myreg.gr <- makeGRangesFromDataFrame(myreg, seqnames.field="V1", start.field="V2", end.field="V3", ignore.strand=TRUE, starts.in.df.are.0based=TRUE, keep.extra.columns=TRUE)
myreg.gr <- makeGRangesFromDataFrame(myreg, ignore.strand=TRUE, starts.in.df.are.0based=TRUE, keep.extra.columns=TRUE)


## JG 8/21/20: Getting Directory to write to 
out_dir <- dirname(opt$infile)

# read in the ensembl annotation bed12 file ( -a )
#  (use unique.prom = FALSE to get the promoter region for each Transcript ID with the ID name included)
#  (use 'remove.unusual=FALSE' for 'weird' genomes such as killifish)
gene.obj <- readTranscriptFeatures(opt$annoBED, up.flank=3000, down.flank=0, unique.prom=FALSE)

# read in the TranscriptID to GeneID df ( -t )
#  used to match TranscriptIDs to GeneIDs, only want one entry per gene promoter
t2g <- read.delim(opt$t2g, header=TRUE)
colnames(t2g) <- c("GeneID", "TranscriptID")

######################################################

# get overlaps, use function from genomation
# strand info is false (default)
foo <- annotateWithGeneParts(myreg.gr, gene.obj)
# get TSS associations
# one entry per cpg, closest TSS matched to specific TranscriptID
mytss <- getAssociationWithTSS(foo)

# get Members (intron, exon, promoter overlaps)
mymems <- getMembers(foo)
# convert the members matrix to data.frame
mymems.df <- as.data.frame(mymems)
# add an "intergenic" column and set all values equal to "0"
mymems.df$intergenic <- 0
# if the prom, exon, and intron cols are all 0, set intergenic col equal to 1
mymems.df <- within(mymems.df, intergenic[prom == 0 & exon == 0 & intron == 0] <- 1)


## Updated



######################################################
## JG 8/24/20: Updating this script to perform the nearest gene and promoter. I dont need sorted files to do this. 
## I will use rtracklayer in order to find the nearest gene. 

# The point of this section seems to be finding the "Gene.overlap.GeneID" for each unqID
#  we end up with mytab1, which is <unqID> <gene.overlap.geneID>

## JG 8/24/20: I will now get the gene overlap with genomation
gene_grs <- c(gene.obj$exons,gene.obj$introns)
gene.hits <- findOverlaps(myreg.gr,gene_grs)

## Getting a dataframe of region_ids and promoter geneids. 
bco <- data.frame(region_id = myreg.gr$region_id[queryHits(gene.hits)], 
                  "TranscriptID" = data.frame(gene_grs)$name[subjectHits(gene.hits)])

# read in the bedtools closest gene overlap output file ( -g )
#bco <- read.delim(opt$bcog, header=FALSE)

# clean up bco, keep columns (1,2,3,4),(5,6,7,8),10,17
# input.coords, gene.coords, strand, distance
#bco <- bco[ ,c(1,2,3,4,5,6,7,8,10,17)]
#colnames(bco) <- c("foo.chr","foo.start","foo.end","region_id","gene.chr","gene.start","gene.end","TranscriptID","gene.strand","Distance")

# reference the t2g table to match TranscriptID to GeneID
bco2 <- merge(bco, t2g, by="TranscriptID", all.x=TRUE)

## JG 8/24/20: Check to see if regions overlap more than one gene 
## There will be more than one overlap it will happen. I must take the first instance for each. 
#bco2 %>% group_by(region_id) %>% count(GeneID) %>% arrange(desc(n))

# Function: Check to see if any input regions overlap more than one gene
# input: df: dataframe (bco2)
#       ID1: column number of unqID
#       ID2: column number of GeneID
#       dis: column number of "distance"
#genecheck <- function(df, ID1, ID2, dis){
# genecheck <- function(df, ID1, ID2){
#     
#       #test <- df[df[ ,dis] == "0", ]				          # grab rows that overlap a Gene THIS IS DONE BY FIND OVERLAP!
#     test <- df
#     test2 <- test[ ,c(ID1,ID2)]	 	       	      	   	          # keep only the unqID and GeneID
#     test3 <- test2[!duplicated(test2), ]             	       	          # remove duplicate rows (where both unqID and GeneID is the same, keeps one)
#     check <- test3[test3[,1]==test3[,1][duplicated(test3[,1])],]         # first try to do below, worked on test but not on real data
#     check <- as.data.frame(test3 %>% group_by(region_id) %>% filter(n() > 1)) # keep rows where unqID appears more than once (and must have a different GeneID per above)
#     return(nrow(check))						          # return the number of rows
# }



# genecheck <- function(df, ID1, ID2){
#   cols <- c(ID1,ID2)
#   test <- df %>% select(all_of(cols)) %>% add_count(region_id) %>% arrange(region_id,desc(n))
#   print(head(test))
#   return(test)						          # return the number of rows
# }

# check to see how many unqIDs overlap more than one GeneID (want it to be 0)
#zero_check <- genecheck(bco2, 5, 11, 10)
#zero_check <- genecheck(bco2, 2, 3)

#noquote(sprintf("Number of unique IDs with repeated GeneIDs: %i", zero_check)) 


# ## JG 2/10/20 Adding automated cleaning feature
# if(zero_check == 0){
#   # If we have 0 
#   bco3 <- bco2
# } else{
#   # Results in zero, so we can remove duplicate unqIDs leaving us with one row per unqID
#   noquote("Removing Duplicate IDS")
#   bco3 <- bco2[!duplicated(bco2$region_id), ]
# }


zero_check <- bco2 %>% count(region_id,GeneID) #%>% count(region_id) %>% arrange(desc(n))
noquote(sprintf("Number of unique IDs with repeated GeneIDs: %i", zero_check %>% count(region_id) %>% filter(n >1) %>% nrow)) 

#bco3 <- zero_check %>% group_by(region_id) %>% slice_max(order_by = region_id,n=1)#%>% top_n(n=1) %>% select(-n) 
bco3 <- zero_check %>% group_by(region_id) %>% slice_max(order_by = region_id,n=1, with_ties = F) %>% dplyr::select(-n) %>% ungroup()
zero_check <- bco3 %>% count(region_id) %>% filter(n >1) %>% nrow
noquote(sprintf("Repeated GeneIDs after processing (should be zero): %i", zero_check))
if(zero_check != 0){
  stop("Gene Overlap Unique DMR Annotation Failed!")
}


# remove unnecessary columns, keep: <unqID><GeneID>
#bco3 <- bco3[ ,c("region_id","GeneID")]
colnames(bco3) <- c("region_id","gene.overlap.GeneID")

# merge with results info

# should probably delete this too, old
######mytab1 <- merge(myres, bco3, by="unqID", all.x=TRUE)

mytab1 <- bco3


##########################################################

# 2) add members info (make this step 1 at some point?)

# this section merges our membership info with our unqIDs based on row number (it works)
# then, it merges this membership info with mytab1 to give us mytab2
# the last step changes "gene.overlap.GeneID" to "NA" if intron and exon are 0, not sure if we want to do this if considering multiple overlaps?


# create a "target.row" column for the input region bed file (myreg) and for mymems.df (made from myreg, so same row order)
myreg$target.row <- rownames(myreg) %>% as.integer()
mymems.df$target.row <- rownames(mymems.df) %>% as.integer()


#head(myreg)
# merge myreg and mymems.df to add unqID to members info
membersmerged <- merge(myreg, mymems.df, by="target.row", all.x=TRUE)
# Remove unnecessary columns (keep only unqID and prom,exon,intron,intergenic)
membersmerged <- membersmerged[ ,c(5,6,7,8,9)]
colnames(membersmerged)[1] <- "region_id"
#head(membersmerged)

# Now, merge this table with mytab1 on "unqID"
#mytab2 <- merge(mytab1, membersmerged, by="region_id", all.x=TRUE)

## I THINK I NEED TO REVERSE THEM THEY ARENT KEEPING THE RIGHT NUMBER OF DMRS.
mytab2 <- merge(membersmerged, mytab1, by="region_id", all.x=TRUE)

# If intron and exon is "0", change "gene.overlap.GeneID" to NA (this column only represents overlaps)
mytab2$gene.overlap.GeneID[mytab2$intron == 0 & mytab2$exon == 0] <- NA


##################################################################################

# 3) add bedtools closest promoter results


# this section adds "promoter.overlap.GeneID" to mytab2

## JG 8/24/20: I will now get the gene overlap with genomation
prom.hits <- findOverlaps(myreg.gr,gene.obj$promoters)

## Getting a dataframe of region_ids and promoter geneids. 
prom.out <- data.frame(region_id = myreg.gr$region_id[queryHits(prom.hits)], 
                  "TranscriptID" = data.frame(gene.obj$promoters)$name[subjectHits(prom.hits)])



# read in the bedtools closest promoters output
#prom.out <- read.delim(opt$bcop, header=FALSE)

## keep only cols 4,8 (unqID TranscriptID, distance)
#prom.out <- prom.out[ ,c(4,8,11)]
#colnames(prom.out) <- c("region_id", "TranscriptID", "dist")


# check to see if any unqIDs overlap more than one TranscriptID
# check to see if any unqIDs overlap more than one GeneID
# merge with t2g
prom.out2 <- merge(prom.out, t2g, by="TranscriptID", all.x=TRUE)
#print(genecheck(prom.out2,2,4,3))
#print(genecheck(prom.out,1,2,3)) # returns 0 so remove duplicate unqIDs, and remove distance column

# check to see how many unqIDs overlap more than one GeneID (want it to be 0)
# zero_check <- genecheck(prom.out2,2,4,3)
# 
# noquote(sprintf("Number of unique IDs with repeated Promoter gene IDS: %i", zero_check)) 
# 
# ## JG 2/10/20 Adding automated cleaning feature
# if(zero_check == 0){
#   # If we have 0 
#   prom.out2 <- prom.out2
# } else{
#   # Results in zero, so we can remove duplicate unqIDs leaving us with one row per unqID
#   noquote("Removing Duplicate Promoter IDS")
#   prom.out2 <- prom.out2[!duplicated(prom.out[,1]), ]
# }

zero_check <- prom.out2 %>% count(region_id,GeneID) #%>% count(region_id) %>% arrange(desc(n))
noquote(sprintf("Number of unique IDs with repeated GeneIDs: %i", zero_check %>% count(region_id) %>% filter(n >1) %>% nrow)) 

#bco3 <- zero_check %>% group_by(region_id) %>% slice_max(order_by = region_id,n=1)#%>% top_n(n=1) %>% select(-n) 
prom.out3 <- zero_check %>% group_by(region_id) %>% slice_max(order_by = region_id,n=1, with_ties = F) %>% select(-n) %>% ungroup()
zero_check <- prom.out3 %>% count(region_id) %>% filter(n >1) %>% nrow
noquote(sprintf("Repeated GeneIDs after processing (should be zero): %i", zero_check))
if(zero_check != 0){
  stop("Promoter Overlap Unique DMR Annotation Failed!")
}

#prom.out2 <- prom.out

#prom.out2 <- prom.out2[ ,c(1,2)]
#colnames(prom.out2) <- c("unqID", "TranscriptID")

# merge with t2g to add GeneID
#prom.out3 <- merge(prom.out2, t2g, by="TranscriptID", all.x=TRUE)
#prom.out3 <- prom.out2[ ,c(2,4)]
colnames(prom.out3)[2] <- "promoter.GeneID"
#colnames(prom.out3) <- c("")

# merge this promoter table with mytab2
mytab3 <- merge(mytab2, prom.out3, by="region_id", all.x=TRUE)


## JG 8/24/20 THE ORDER IS OFF!!!
# if the "prom" column is "0" change promoter.GeneID to NA (this code results in a warning but generates <NA> values)
mytab3$promoter.GeneID[mytab3$prom == 0] <- NA



###########################################################

# 4) add distance to nearest TSS for regions that do not overlap genes

# refer to mytss and myreg with their target.row columns created above

# merge myreg and mytss on target.row
mergedtss <- merge(myreg, mytss, by="target.row", all.x=TRUE)
# remove unnecessary columns (keep only unqID, dist.to.feature, feature.name)
mergedtss <- mergedtss[ ,c(5,6,7)]
colnames(mergedtss) <- c("region_id", "dist.to.TSS", "TranscriptID")

# merge with t2g to add GeneID
mergedtss2 <- merge(mergedtss, t2g, by="TranscriptID", all.x=TRUE)
mergedtss2 <- mergedtss2[ ,c(2,3,4)]
colnames(mergedtss2)[3] <- "tss.GeneID"

# merge with mytab3 on "unqID" NOW region_id
mytab4 <- merge(mytab3, mergedtss2, by="region_id", all.x=TRUE)

# if intron or exon column equal 1, then change tss.GeneID to NA
mytab4$tss.GeneID[mytab4$intron == 1 | mytab4$exon == 1] <- NA
# also change "dist.to.TSS" to NA
mytab4$dist.to.tss[mytab4$intron == 1 | mytab4$exon == 1] <- NA 


# rearrange columns and export here: need to add biomart info
#write.table(mytab4, "mytab4.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

# reference is here: /u1/bd/Projects/ECP7/ECP7.2/ANNOTATE/DMC/try2

##################################################################################################

# 5) Add ensembl biomart gene info for the 3 sets of GeneIDs

# read in biomart reference
## Change to be the format
#tsv_names = c("Gene stable ID",	"Transcript stable ID",	"Source of gene name"	,"Gene type",	"Gene Name",	"Description",	"Ensembl Family Description")
#tsv_names = c("GeneID",	"TranscriptID",	"source.of.gene.name"	,"gene.type",	"gene.name",	"description",	"ensembl.family.desc")

## JG 8/24/20: Issue is that downloaded mart doesnt carry other sources (ensembl_havana, etc.) I will make a system to read in from BiomaRt. 
## Now getting all unique gene IDS to query biomart with package. 
all_gene_ids <- mytab4 %>% dplyr::select(contains("GeneID")) %>% as.list() %>% unlist() %>% as.character()
all_gene_ids <- all_gene_ids[!is.na(all_gene_ids)] %>% unique()

print(sprintf("Number of Unique gene IDs Identified: %i", length(all_gene_ids)))


## JG 8/24/20: Now initalizing biomart.
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
## CHANGE HERE IF YOU WANT OTHER ANNOTATION DATA.
the_attrib <- c("ensembl_gene_id","external_gene_name","source","gene_biotype","description","family_description")
#new_names <- 
search_by <- "ensembl_gene_id"

ref <- getBM(attributes = the_attrib,
             filters = search_by,
             values = all_gene_ids,
             mart = ensembl)

# ## There are duplicate entires when familes have more than one entry... Now merging. 
# paste_desc <- function(desc_vect){
#   print(desc_vect)
#   ## If it is just one thing return it
#   if(length(desc_vect == 1)){
#     return(desc_vect)
#   } else{
#     #print(length(desc_vect))
#     desc_vect <- desc_vect[! desc_vect == ""]
#     print(desc_vect)
#     #print(length(desc_vect))
#     if(length(desc_vect) > 1){
#       out <- paste(desc_vect, collapse = "|")
#     } else{
#       out <- desc_vect
#     }
#     return(out)
#   }## End of else there are at least 2 things.
# }
# 
# funct <- function(x) return(paste_desc(x))
#ref %>% group_by(ensembl_gene_id) %>% summarize(family_description = funct(family_description)) ## THIS DIDNT WORK

## JG 8/24/20 THERE IS ALOT OF ROOM FOR IMPROVEMENT OF THIS FUNCTION... I dont know how to make the function ignore "" and just return that.
ref <- ref %>% group_by(ensembl_gene_id) %>% mutate(family_description = paste(family_description, collapse = "|")) %>% ungroup() %>% unique()


#ref <- read.delim(opt$mart, header=TRUE, stringsAsFactors = F, na.strings = "", col.names = tsv_names)
#colnames(ref)[1] <- "GeneID"

## Removing the Transcript ID
#ref <- ref %>% dplyr::select(-TranscriptID) %>% dplyr::select(GeneID,gene.name,everything()) 

# print("BIOMART DATA")
# head(ref)
# stop()


# remove duplicate GeneID entries
#ref <- ref[!duplicated(ref$GeneID), ]

# Merge with gene.overlap.GeneID
# change colnames of ref to reflect
colnames(ref) <- c("gene.overlap.GeneID","gene.overlap.gene.name","gene.overlap.source.of.gene.name","gene.overlap.gene.type","gene.overlap.gene.description","gene.overlap.Ensembl.family.description")
# merge
mytab5 <- merge(mytab4, ref, by="gene.overlap.GeneID", all.x=TRUE)


# Merge with promoter.GeneID
# change colnames of ref to reflect
colnames(ref) <- gsub("gene.overlap", "promoter", colnames(ref))


# merge
mytab6 <- merge(mytab5, ref, by="promoter.GeneID", all.x=TRUE)


# merge with tss.GeneID
# change colnames of ref to reflect
colnames(ref) <- gsub("promoter", "tss", colnames(ref))
# merge
mytab7 <- merge(mytab6, ref, by="tss.GeneID", all.x=TRUE)




################################################################################################################################

# reorder columns
#mytab7 <- mytab7[ ,c(4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,3,2,20,1,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35)]

# reorder rows by chr and start
#colnames(mytab7)[2] <- "chr"
#mytab7 <- mytab7[order(mytab7[ ,2], mytab7[ ,3]), ]


# check for duplicate rows to remove
#  these rows to remove will have the same value in 4 columns:
#  <unqID> <tss.GeneID> <promoter.GeneID> <gene.overlap.GeneID>
#mytab8 <- mytab7[!duplicated(mytab7[c("region_id","tss.GeneID","promoter.GeneID","gene.overlap.GeneID")]), ]
#mytab8 <- mytab7

# print("MYTAB8")
# head(mytab8)
# stop()

# EXPORT FINAL TABLE
#write.table(mytab7, paste(opt$name, "mytab7.txt", sep="."), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
## JG 8/21/20: Adding New Name function
#write.table(mytab8, paste(opt$name, "mytab8.txt", sep="."), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
#out_name <- sprintf("%s/%s.mytab8.txt",out_dir,opt$name)
#write.table(mytab8, out_name, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)



## JG 8/22/20: Create and write out xslx output

## This should be mytab8
#head(mytab8)
#mytab <- read.delim(infile, header=TRUE, stringsAsFactors = F)
mytab <- mytab7
#head(mytab)
#dim(mytab) 
# change "unqID" to "DMR_ID" for proper merging
#colnames(mytab)[4] <- "region_id"

#mytab$region_id <- gsub(pattern = '_', replacement = ".", mytab$region_id)

# read in the corresponding results file
## JG 8/21/20 Open DMR file

res <- read.delim(opt$dmr, header=TRUE, stringsAsFactors = F)
## DMR FILE IS INFILE
#res <- myreg
#dim(res)


# merge on DMR_ID
#print("MERGE! 1")
foo <- merge(res, mytab, by="region_id", all.x=TRUE)

#print("MERGE WORKED!")
## Add absolute difference for excel mapping
foo <- foo %>% mutate(abs_diff = abs(adj_diff_meth)) %>% mutate(dm_type = case_when( adj_diff_meth > 0  ~ "hyper",
                                                                                     adj_diff_meth < 0 ~ "hypo"))

#head(foo)
#stop()

# reorder columns
#"DMR_ID	chr	start	end	pvalue	qvalue	meth.diff	prom	exon	intron	intergenic	gene.overlap.GeneID	promoter.GeneID	dist.to.TSS	tss.GeneID	gene.overlap.gene.name	gene.overlap.source.of.gene.name	gene.overlap.gene.type	gene.overlap.gene.description	gene.overlap.Ensembl.family.description	promoter.gene.name	promoter.source.of.gene.name	promoter.gene.type	promoter.gene.description	promoter.Ensembl.family.description	tss.gene.name	tss.source.of.gene.name	tss.gene.type	tss.gene.description	tss.Ensembl.family.description"
ord1 <- c("region_id", "pvalue", "FDRp", "dm_type","adj_diff_meth","meth_diff","abs_diff","case_mean","ctrl_mean")
ord2 <- c('prom',	'exon',	'intron',	'intergenic')
ord3 <- c('gene.overlap.GeneID',	'promoter.GeneID',	'dist.to.TSS',	'tss.GeneID')
ord4 <- c('gene.overlap.gene.name',	'gene.overlap.source.of.gene.name',	'gene.overlap.gene.type',	'gene.overlap.gene.description',	'gene.overlap.Ensembl.family.description')
ord5 <- c('promoter.gene.name',	'promoter.source.of.gene.name',	'promoter.gene.type',	'promoter.gene.description',	'promoter.Ensembl.family.description')
ord6 <- c('tss.gene.name',	'tss.source.of.gene.name',	'tss.gene.type',	'tss.gene.description',	'tss.Ensembl.family.description')
final_order <- c(ord1,ord2,ord3,ord4,ord5,ord6)
 
 
foo2 <- foo[ ,final_order]
 
 
to_parse <- foo2 %>% dplyr::select(region_id)
parsed <- to_parse %>% tidyr::separate(region_id, c('chr','start','stop'))
 
foo2 <- cbind(foo2, parsed)

foo2 <- foo2 %>% dplyr::select(region_id, chr, start, stop, everything())

# add "chr" to chr name
#foo2$chr <- paste0("chr", foo2$chr)

# sort numerically by DMR_ID
#foo3 <- foo2
#foo3$sort <- gsub('DMR_', '', foo3$region_id)
#foo4 <- foo3[order(as.numeric(as.character(foo3$sort))), ]
#foo4 <- foo4 %>% dplyr::select(-sort)
foo4 <- foo2 %>% arrange(region_id)
# 
# # export
write.table(foo4, "dmr.sig.final_table.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

# I want to make this file 'pretty' to be easily analyzed by a biologist.
# For each subtype I will keep only the relevant info

pretty_cols <- c("region_id", "chr","FDRp", "dm_type","abs_diff","meth_diff","case_mean","ctrl_mean")
#ord2 #REMINDER THAT IS IN THE BINARY MATRIX OF THE OVERLAP TYPE
prom_info <- c('promoter.GeneID','promoter.gene.name',	'promoter.source.of.gene.name',	'promoter.gene.type',	'promoter.gene.description',	'promoter.Ensembl.family.description')
exon_intron_info <- c('gene.overlap.GeneID', 'gene.overlap.gene.name',	'gene.overlap.source.of.gene.name',	'gene.overlap.gene.type',	'gene.overlap.gene.description',	'gene.overlap.Ensembl.family.description')
intergenic_info <-  c('dist.to.TSS',	'tss.GeneID', 'tss.gene.name',	'tss.source.of.gene.name',	'tss.gene.type',	'tss.gene.description',	'tss.Ensembl.family.description')
simple_genic <- c('gene.overlap.GeneID', 'gene.overlap.gene.name','gene.overlap.gene.type','gene.overlap.gene.description')

## export straight to excel file
# Must create named list.
# Removing empty comlumns to make pretty! https://stackoverflow.com/questions/15968494/how-to-delete-columns-that-contain-only-nas
foo4 <- foo4 %>% arrange(desc(abs_diff), FDRp)
prom <- foo4 %>% filter(prom == 1) %>% select_if(~!all(is.na(.))) %>% dplyr::select(pretty_cols,ord2, prom_info, simple_genic)
exon_intron <- foo4 %>% filter(exon == 1 | intron == 1) %>% dplyr::select_if(~!all(is.na(.))) %>% dplyr::select(pretty_cols, ord2,exon_intron_info)
intergenic <- foo4 %>% filter(intergenic == 1) %>% dplyr::select_if(~!all(is.na(.))) %>% dplyr::select(pretty_cols,ord2, intergenic_info)
hypo <- foo4 %>% filter(dm_type == "hypo")
hyper <- foo4 %>% filter(dm_type == "hyper")

## Make outfile
out_name <- sprintf("%s/%s.xlsx",out_dir ,opt$name)

#to_export <- list(full_dmrs = foo4, promoters = prom, exon_intron = exon_intron, itergenic = intergenic)
to_export <- list(full_dmrs = foo4, promoters = prom, exon_intron = exon_intron, itergenic = intergenic,hypo=hypo,hyper=hyper)
# print(to_export %>% names())
# #write_xlsx(to_export, 'TS_BAV_DMR_Final.xlsx')
write_xlsx(to_export, out_name)
print("Excel File Written! Now Done.")
# 
# to_export <- list(full_regions = foo4 %>% arrange(desc(abs_diff), FDRp), 
#                   promoters = prom %>% arrange(desc(abs_diff), FDRp), 
#                   exon_intron = exon_intron %>% arrange(desc(abs_diff), FDRp),
#                   itergenic = intergenic %>% arrange(desc(abs_diff), FDRp))
# 
# ## Write sorted.
# out_name <- sprintf("%s/Sorted_%s.xlsx",out_dir ,opt$name)
# #write_xlsx(to_export, 'Sorted_TS_BAV_DMR_Final.xlsx')
# write_xlsx(to_export, out_name)

