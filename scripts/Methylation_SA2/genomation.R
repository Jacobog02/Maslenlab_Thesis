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



#-------------------------------------------------------------------------------------------------------#
# HOW TO RUN:												#
#													#
# Command Line Arguments										#
#													#
# 1. -i, --infile: input regions BED file in format <chr><start><end><unqID>				#
#                   must contain canonical contigs only, no header					#
# 2. -a, --annoBED: path to ensembl annotation BED12 file   						#
# 3. -t, --t2g: path to t2g reference file: <GeneID> <TranscriptID>, no header				#
# 5. -g, --bcog: bedtools closest gene output file   		     					#
# 6. -p, --bcop: bedtools closest promoter output file							#
# 7. -b, --mart: path to bioMart gene info ref file:							#
#        <GeneID><gene.name><source.of.gene.name><gene.type><gene.description><Ensembl.family.desc>	#
# 8. -n, --name: name at the beginning of the output table						#
#-------------------------------------------------------------------------------------------------------#

library("optparse")

#----------------------#
# COMMAND LINE OPTIONS #
#----------------------#
option_list = list(
  make_option(c("-i", "--infile"), type="character", default=NULL, 
              help="path to input regions bed file: <chr> <start> <end> <unqID>", metavar="character"),
  make_option(c("-a", "--annoBED"), type="character", default=NULL, 
	      help="path to ensembl annotation BED12 file", metavar="character"),
  make_option(c("-t", "--t2g"), type="character", default=NULL,
  	      help="path to t2g reference file: <GeneID> <TranscriptID>, no header", metavar="character"),
  make_option(c("-g", "--bcog"), type="character", default=NULL,
  	      help="bedtools closest gene output file", metavar="character"),
  make_option(c("-p", "--bcop"), type="character", default=NULL,
              help="bedtools closest promoter output file", metavar="character"),
  make_option(c("-b", "--mart"), type="character", default=NULL,
  	      help="path to bioMart gene info ref file: <GeneID><gene.name><source.of.gene.name><gene.type><gene.description><Ensembl.family.desc>", metavar="character"),
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
if (is.null(opt$bcog)){
  print_help(opt_parser)
  stop("Missing bedtools closest gene output file.n", call.=FALSE)
}
if (is.null(opt$bcop)){
  print_help(opt_parser)
  stop("Missing bedtools closest promoter output file.n", call.=FALSE)
}
if (is.null(opt$mart)){
  print_help(opt_parser)
  stop("Missing bioMart gene info file.n", call.=FALSE)
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

library(genomation)
library(GenomicRanges)
library(dplyr)

# read in bed file of input regions ( -i )
myreg <- read.table(opt$infile, header=FALSE)
# convert to GRanges object
myreg.gr <- makeGRangesFromDataFrame(myreg, seqnames.field="V1", start.field="V2", end.field="V3", ignore.strand=TRUE, starts.in.df.are.0based=TRUE, keep.extra.columns=TRUE)

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

######################################################

# MERGING #

# this has been commented out for a while (even in working versions), should probably delete
###### 1) Merge results with bedtools closest gene output
###### read in results file on input regions (first column should be the unique ID)
######myres <- read.delim("/u1/bd/Projects/ECP9/ALL_ATAC/macs2/individual/p05_spmr/diffBind/new2/vitroPos_vs_vitroNeg/vitroPos_vs_vitroNeg.ALL.results.can.txt", header=TRUE)
######myres <- read.delim(opt$res, header=TRUE)
######colnames(myres)[1] <- "unqID"



# The point of this section seems to be finding the "Gene.overlap.GeneID" for each unqID
#  we end up with mytab1, which is <unqID> <gene.overlap.geneID>


# read in the bedtools closest gene overlap output file ( -g )
bco <- read.delim(opt$bcog, header=FALSE)

# clean up bco, keep columns (1,2,3,4),(5,6,7,8),10,17
# input.coords, gene.coords, strand, distance
bco <- bco[ ,c(1,2,3,4,5,6,7,8,10,17)]
colnames(bco) <- c("foo.chr","foo.start","foo.end","unqID","gene.chr","gene.start","gene.end","TranscriptID","gene.strand","Distance")

# reference the t2g table to match TranscriptID to GeneID
bco2 <- merge(bco, t2g, by="TranscriptID", all.x=TRUE)

# Function: Check to see if any input regions overlap more than one gene
# input: df: dataframe (bco2)
#       ID1: column number of unqID
#       ID2: column number of GeneID
#       dis: column number of "distance"
genecheck <- function(df, ID1, ID2, dis){
    test <- df[df[ ,dis] == "0", ]				          # grab rows that overlap a Gene
    test2 <- test[ ,c(ID1,ID2)]	 	       	      	   	          # keep only the unqID and GeneID
    test3 <- test2[!duplicated(test2), ]             	       	          # remove duplicate rows (where both unqID and GeneID is the same, keeps one)
    #check <- test3[test3[,1]==test3[,1][duplicated(test3[,1])],]         # first try to do below, worked on test but not on real data
    check <- as.data.frame(test3 %>% group_by(unqID) %>% filter(n() > 1)) # keep rows where unqID appears more than once (and must have a different GeneID per above)
    return(nrow(check))						          # return the number of rows
}

# check to see how many unqIDs overlap more than one GeneID (want it to be 0)
zero_check <- genecheck(bco2, 5, 11, 10)

noquote(sprintf("Number of unique IDs with repeated GeneIDs: %i", zero_check)) 

## JG 2/10/20 Adding automated cleaning feature
if(zero_check == 0){
  # If we have 0 
  bco3 <- bco2
} else{
  # Results in zero, so we can remove duplicate unqIDs leaving us with one row per unqID
  noquote("Removing Duplicate IDS")
  bco3 <- bco2[!duplicated(bco2$unqID), ]
}



# remove unnecessary columns, keep: <unqID><GeneID>
bco3 <- bco3[ ,c("unqID","GeneID")]
colnames(bco3) <- c("unqID","gene.overlap.GeneID")

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

# merge myreg and mymems.df to add unqID to members info
membersmerged <- merge(myreg, mymems.df, by="target.row", all.x=TRUE)
# Remove unnecessary columns (keep only unqID and prom,exon,intron,intergenic)
membersmerged <- membersmerged[ ,c(5,6,7,8,9)]
colnames(membersmerged)[1] <- "unqID"

# Now, merge this table with mytab1 on "unqID"
mytab2 <- merge(mytab1, membersmerged, by="unqID", all.x=TRUE)

# If intron and exon is "0", change "gene.overlap.GeneID" to NA (this column only represents overlaps)
mytab2$gene.overlap.GeneID[mytab2$intron == 0 & mytab2$exon == 0] <- NA

##################################################################################

# 3) add bedtools closest promoter results


# this section adds "promoter.overlap.GeneID" to mytab2


# read in the bedtools closest promoters output
prom.out <- read.delim(opt$bcop, header=FALSE)

# keep only cols 4,8 (unqID TranscriptID, distance)
prom.out <- prom.out[ ,c(4,8,11)]
colnames(prom.out) <- c("unqID", "TranscriptID", "dist")

# check to see if any unqIDs overlap more than one TranscriptID
# check to see if any unqIDs overlap more than one GeneID
# merge with t2g
prom.out2 <- merge(prom.out, t2g, by="TranscriptID", all.x=TRUE)
print(genecheck(prom.out2,2,4,3))
#print(genecheck(prom.out,1,2,3)) # returns 0 so remove duplicate unqIDs, and remove distance column

# check to see how many unqIDs overlap more than one GeneID (want it to be 0)
zero_check <- genecheck(prom.out2,2,4,3)

noquote(sprintf("Number of unique IDs with repeated Promoter gene IDS: %i", zero_check)) 

## JG 2/10/20 Adding automated cleaning feature
if(zero_check == 0){
  # If we have 0 
  #prom.out2 <- prom.out2
} else{
  # Results in zero, so we can remove duplicate unqIDs leaving us with one row per unqID
  noquote("Removing Duplicate Promoter IDS")
  prom.out2 <- prom.out2[!duplicated(prom.out[,1]), ]
}


#prom.out2 <- prom.out

#prom.out2 <- prom.out2[ ,c(1,2)]
#colnames(prom.out2) <- c("unqID", "TranscriptID")

# merge with t2g to add GeneID
#prom.out3 <- merge(prom.out2, t2g, by="TranscriptID", all.x=TRUE)
prom.out3 <- prom.out2[ ,c(2,4)]
colnames(prom.out3)[2] <- "promoter.GeneID"

# merge this promoter table with mytab2
mytab3 <- merge(mytab2, prom.out3, by="unqID", all.x=TRUE)

# if the "prom" column is "0" change promoter.GeneID to NA (this code results in a warning but generates <NA> values)
mytab3$promoter.GeneID[mytab3$prom == 0] <- NA

###########################################################

# 4) add distance to nearest TSS for regions that do not overlap genes

# refer to mytss and myreg with their target.row columns created above

# merge myreg and mytss on target.row
mergedtss <- merge(myreg, mytss, by="target.row", all.x=TRUE)
# remove unnecessary columns (keep only unqID, dist.to.feature, feature.name)
mergedtss <- mergedtss[ ,c(5,6,7)]
colnames(mergedtss) <- c("unqID", "dist.to.TSS", "TranscriptID")

# merge with t2g to add GeneID
mergedtss2 <- merge(mergedtss, t2g, by="TranscriptID", all.x=TRUE)
mergedtss2 <- mergedtss2[ ,c(2,3,4)]
colnames(mergedtss2)[3] <- "tss.GeneID"

# merge with mytab3 on "unqID"
mytab4 <- merge(mytab3, mergedtss2, by="unqID", all.x=TRUE)

# if intron or exon column equal 1, then change tss.GeneID to NA
mytab4$tss.GeneID[mytab4$intron == 1 | mytab4$exon == 1] <- NA
# also change "dist.to.TSS" to NA
mytab4$dist.to.TSS[mytab4$intron == 1 | mytab4$exon == 1] <- NA 


# rearrange columns and export here: need to add biomart info
#write.table(mytab4, "mytab4.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

# reference is here: /u1/bd/Projects/ECP7/ECP7.2/ANNOTATE/DMC/try2

##################################################################################################

# 5) Add ensembl biomart gene info for the 3 sets of GeneIDs

# read in biomart reference
ref <- read.delim(opt$mart, header=TRUE, stringsAsFactors = F, na.strings = "")
colnames(ref)[1] <- "GeneID"

# remove duplicate GeneID entries
ref <- ref[!duplicated(ref$GeneID), ]

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
mytab8 <- mytab7[!duplicated(mytab7[c("unqID","tss.GeneID","promoter.GeneID","gene.overlap.GeneID")]), ]


# EXPORT FINAL TABLE
#write.table(mytab7, paste(opt$name, "mytab7.txt", sep="."), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
## JG 8/21/20: Adding New Name function
#write.table(mytab8, paste(opt$name, "mytab8.txt", sep="."), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
out_name <- sprintf("%s/%s.mytab8.txt",out_dir,opt$name)
write.table(mytab8, out_name, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
