#!/usr/bin/Rscript
#
# Load salmon produced transcript counts 
# Expects one file:
#
#  samples.txt, sample information,
#    col1 will be used as colnames ... colLAST specifies path to salmon quant.sf files 
# e.g.
#  
#    SRR4048970 control rep1 ./quants/SRR4048970/quant.sf 
#    SRR4048971 control rep2 ./quants/SRR4048971/quant.sf 
#    SRR4048972 control rep3 ./quants/SRR4048972/quant.sf 
#    SRR4048973 15d-PGJ2 rep1 ./quants/SRR4048973/quant.sf 
#    SRR4048974 15d-PGJ2 rep2 ./quants/SRR4048974/quant.sf 
#    SRR4048975 15d-PGJ2 rep3 ./quants/SRR4048975/quant.sf
#
time <- format(Sys.time(), "%Y%m%d-%H%M%S") 
library(tximport) 
library(methods) 
library(readr) 
library(DESeq2) 
library("BiocParallel") 

args <- commandArgs(trailingOnly = TRUE) 
samplefile <- args[1] 
mapfile <- args[2] 
design <- args[3]
outfile1 <- args[4] 
outfile2 <- args[5] 
threads <- as.integer(args[6]) 
register(MulticoreParam(threads))


design_matrix<- read.csv(file.path(design))
design_matrix$cond<-as.factor(design_matrix$cond)
design_matrix$pair<-as.factor(design_matrix$pair)
design_matrix$cond<-relevel(design_matrix$cond, ref="pre")

# samples.txt, e.g.
samples <- read.delim(file.path(samplefile), header = FALSE)


## tx_gene_map.txt e.g.
t2g <- read.delim(file.path(mapfile), header = FALSE)


# load quant files specified in sample_list.txt col4 are always in last column
files <- file.path(samples[,ncol(samples)])


# estimate adjusted counts from TPMs 
# Based on tximport: 
# https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html
# I ammended this from $coutns to $abundance as I think this is how we get the TPMs https://support.bioconductor.org/p/84883/
txi.tx <- tximport(files = files, type="salmon", txOut = TRUE) 
colnames(txi.tx$abundance) <- samples[,1]
y<-as.data.frame(txi.tx$abundance)



# use counts with DESeqDataSetFromMatrix but DESeqDataSetFromTximport takes care of it and does not do anything to the counts if countsFromAbundance is either "scaledTPM" or lengthScaledTPM also look at 
# http://rdrr.io/bioc/DESeq2/src/R/AllClasses.R => DESeqDataSetFromTximport
## tx ==========================================================================
#dds.tx <- DESeqDataSetFromTximport(txi.tx, colData = samples, design = ~1)
#dds.tx <- DESeq(dds.tx, parallel=T)
#dds.tx <- estimateSizeFactors(dds.tx) 
#y <- counts(dds.tx, normalized=T)
#dds.tx.counts <- rlog(dds.tx, blind=T) 
#y <- assay(dds.tx.counts)


# problem with "X" being prepended in column variables when sample names start with number ... solution print header first data after
gz1 <- gzfile(outfile1, "w") 
#header <- t(c("tx", as.character(colnames(y)))) 
#write.table(header, file=gz1, row.names=FALSE, col.names=FALSE, append=FALSE, quote=FALSE, sep="\t") 
write.table(y, file=gz1, append=FALSE, quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE) 
close(gz1)



## genes =======================================================================
## I ammended this to extract the tpms instead of the counts in order to run the rule get_tx_stats

txi <- summarizeToGene(txi.tx, tx2gene=t2g, ignoreTxVersion = FALSE) 
colnames(txi$abundance) <- samples[,1]
y <- as.data.frame(txi$abundance)

#colnames(txi$counts) <- samples[,1] 

#dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~1)


#dds <- DESeq(dds, parallel=T)
#dds <- estimateSizeFactors(dds) 
#y <- counts(dds, normalized=T)
#dds.counts <- rlog(dds, blind=T) 
#y <- assay(dds.counts)

gz2 <- gzfile(outfile2, "w") 
#header <- t(c("gene", as.character(colnames(y)))) 
#write.table(header, file=gz2, row.names=FALSE, col.names=FALSE, append=FALSE, quote=FALSE, sep="\t") 
write.table(y, file=gz2, append=FALSE, quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE) 
close(gz2)
