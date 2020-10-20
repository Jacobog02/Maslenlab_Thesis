#!/usr/bin/Rscript

#-------------------------------------------------------------------------------------------------------#
# HOW TO RUN:												#
#													#
# Command Line Arguments										#
#													#
# 1. -i, --infile: input regions BED file in format <chr><start><end><unqID>				#
#                   must contain canonical contigs only, no header					#
#-------------------------------------------------------------------------------------------------------#

library("optparse")

#----------------------#
# COMMAND LINE OPTIONS #
#----------------------#
option_list = list(
  make_option(c("-i", "--infile"), type="character", default=NULL, 
              help="path to input regions bed file: <chr> <start> <end> <unqID>", metavar="character"),
  make_option(c("-d", "--dmr"), type="character", default=NULL, 
              help="path to input regions bed file: <chr> <start> <end> <unqID>", metavar="character"),
  make_option(c("-n", "--name"), type="character", default=NULL, 
              help="Name to Use for Output", metavar="character")

); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#-----------------------------------#
# COMMAND LINE OPTION SANITY CHECKS #
#-----------------------------------#
if (is.null(opt$infile)){
  print_help(opt_parser)
  stop("Missing input annotation file.n", call.=FALSE)
}
if (is.null(opt$dmr)){
  print_help(opt_parser)
  stop("Missing input dmr file.n", call.=FALSE)
}
## Script starts here. 
# after "run_annot.sh" which runs the genomation.R script
require(dplyr)
require(writexl)

 # read in the "mytab8" output
infile <- opt$infile
out_dir <- dirname(infile)


mytab <- read.delim(infile, header=TRUE, stringsAsFactors = F)
# change "unqID" to "DMR_ID" for proper merging
colnames(mytab)[4] <- "region_id"

mytab$region_id <- gsub(pattern = '_', replacement = ".", mytab$region_id)

print(head(mytab))
# read in the corresponding results file
## JG 8/21/20 Open DMR file

res <- read.delim(opt$dmr, header=TRUE, stringsAsFactors = F)
print(head(res))
# merge on DMR_ID
print("MERGE! 1")
foo <- merge(mytab, res, by="region_id", all.x=TRUE)

print("MERGE WORKED!")
## Add absolute difference for excel mapping
foo <- foo %>% mutate(abs_diff = abs(adj_diff_meth)) %>% mutate(dm_type = case_when( adj_diff_meth > 0  ~ "hyper",
                                                                                     adj_diff_meth < 0 ~ "hypo"))

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


to_parse <- foo2 %>% select(region_id)
parsed <- to_parse %>% tidyr::separate(region_id, c('chr','start','stop'))

foo2 <- cbind(foo2, parsed)

foo2 <- foo2 %>% select(region_id, chr, start, stop, everything())

# add "chr" to chr name
foo2$chr <- paste0("chr", foo2$chr)

# sort numerically by DMR_ID
foo3 <- foo2
foo3$sort <- gsub('DMR_', '', foo3$region_id)
foo4 <- foo3[order(as.numeric(as.character(foo3$sort))), ]
foo4 <- foo4 %>% select(-sort)

# export
#write.table(foo4, "dmr.sig.final_table.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

# I want to make this file 'pretty' to be easily analyzed by a biologist.
# For each subtype I will keep only the relevant info 

pretty_cols <- c("region_id", "chr","FDRp", "dm_type","abs_diff","meth_diff","case_mean","ctrl_mean")
#ord2 #REMINDER THAT IS IN THE BINARY MATRIX OF THE OVERLAP TYPE
prom_info <- c('promoter.GeneID','promoter.gene.name',	'promoter.source.of.gene.name',	'promoter.gene.type',	'promoter.gene.description',	'promoter.Ensembl.family.description')
exon_intron_info <- c('gene.overlap.GeneID', 'gene.overlap.gene.name',	'gene.overlap.source.of.gene.name',	'gene.overlap.gene.type',	'gene.overlap.gene.description',	'gene.overlap.Ensembl.family.description')
intergenic_info <-  c('dist.to.TSS',	'tss.GeneID', 'tss.gene.name',	'tss.source.of.gene.name',	'tss.gene.type',	'tss.gene.description',	'tss.Ensembl.family.description')
simple_genic <- c('gene.overlap.GeneID', 'gene.overlap.gene.name','gene.overlap.gene.type','gene.overlap.gene.description')

## export straight to excel file 
## Must create named list. 
## Removing empty comlumns to make pretty! https://stackoverflow.com/questions/15968494/how-to-delete-columns-that-contain-only-nas
prom <- foo4 %>% filter(prom == 1) %>% select_if(~!all(is.na(.))) %>% select(pretty_cols,ord2, prom_info, simple_genic)
exon_intron <- foo4 %>% filter(exon == 1 | intron == 1) %>% select_if(~!all(is.na(.))) %>% select(pretty_cols, ord2,exon_intron_info)
intergenic <- foo4 %>% filter(intergenic == 1) %>% select_if(~!all(is.na(.))) %>% select(pretty_cols,ord2, intergenic_info)


out_name <- sprintf("%s/%s.xlsx",out_dir ,opt$name)

to_export <- list(full_regions = foo4, promoters = prom, exon_intron = exon_intron, itergenic = intergenic)
print(to_export %>% names())
#write_xlsx(to_export, 'TS_BAV_DMR_Final.xlsx')
write_xlsx(to_export, out_name)

to_export <- list(full_regions = foo4 %>% arrange(desc(abs_diff), FDRp), 
                  promoters = prom %>% arrange(desc(abs_diff), FDRp), 
                  exon_intron = exon_intron %>% arrange(desc(abs_diff), FDRp),
                  itergenic = intergenic %>% arrange(desc(abs_diff), FDRp))

## Write sorted.
out_name <- sprintf("%s/Sorted_%s.xlsx",out_dir ,opt$name)
#write_xlsx(to_export, 'Sorted_TS_BAV_DMR_Final.xlsx')
write_xlsx(to_export, out_name)
# How many Sig DMRs?
#nrow(foo2[foo2$Sidak_Pval < 0.1 & abs(foo2$perc_meth_change) > 10, ])

