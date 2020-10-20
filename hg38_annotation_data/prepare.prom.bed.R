## Jacob Gutierrez gutierja@ohsu.edu
## 2/6/20
## THIS IS NOT AN ORIGINAL SCRIPT I AM NOT PRETENDING IT IS!!!!
## Copied from /home/groups/hoolock2/u0/bd/my_scripts/methylation_annot/prepare.prom.bed.R THANK YOU BRETT!

## JG 8/21/20: I changed the annotation file so I updated it.

# Purpose: use the bed12 file (made from the gtf file) to make a bed file of promoter regions

library(genomation)
library(GenomicRanges)

# read in the ensembl annotation bed12 file
#  use unique.prom = FALSE to get the promoter region for each Transcript ID with the ID name included
gene.obj <- readTranscriptFeatures("ensembl_GRCh38.98_ann.bed",
				   up.flank=3000,
				   down.flank=0,
				   unique.prom=FALSE
)

# this results in a GRangesList object holding the exons, introns, promoters, and tss's

# convert the promoters GRanges object to a data frame in bed format
prom.bed <- data.frame(seqnames=seqnames(gene.obj$promoters),
		       starts=start(gene.obj$promoters),
		       ends=end(gene.obj$promoters),
		       names=elementMetadata(gene.obj$promoters)$name,
		       scores=elementMetadata(gene.obj$promoters)$score,
		       strands=strand(gene.obj$promoters)
)

# ensure that there is no scientific notation
#prom.bed2 <- format(prom.bed, scientific=FALSE)
prom.bed2 <- prom.bed

# ## Maybe issue with character vectors
# prom.bed2$starts <- as.numeric(prom.bed2$starts)
# prom.bed2$ends <- as.numeric(prom.bed2$ends)

# convert any negative start positions to 0 NOT ZERO 1 ????
## THE as.numeric() was added 2/7/20 bc it said that everything started at 0 bc the comarpsion was a character against a numeric zero
prom.bed2$starts[prom.bed2$starts < 0] <- 0 # trying 1 

# convert any end positions that are 0 to 1
prom.bed2$ends[prom.bed2$ends == 0] <- '1'


# ensure that there is no scientific notation
# prom.bed2$starts <- format(prom.bed2$starts, scientific=FALSE)
# 
# prom.bed2$ends <- format(prom.bed2$ends, scientific=FALSE)

## Trying formatting now. 
prom.bed2 <- format(prom.bed2, scientific=FALSE, trim = TRUE)

# export the promoter bed file
write.table(prom.bed2, "ensembl_GRCh38.98_promoter.bed", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

