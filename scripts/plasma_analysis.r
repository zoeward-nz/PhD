setwd("Z:/HiSeq_OG4071")
#https://stackoverflow.com/questions/24129269/compile-r-script-to-markdown
### load the libraries
library(BiocParallel)
library(tximport)
library(readr)
library("DESeq2")
library("ggplot2")
library(TissueEnrich)
library('dplyr')

###########################################################
##### Import some files for DESEQ analysis
###########################################################

class_code_transcriptlist <- read.delim("analysis/results/04_novel/class_code_list.txt", header=F, stringsAsFactors=FALSE)
colnames(class_code_transcriptlist)<- c("class_code","gene_name","transcript_id", "number_of_exons")

## import gencode info 
## this data is taken from the gencode.gtf (the co-ordinates, ensemble gene names, gene type and gene name)
gencode_lncrna<-read.delim("STAR_GRCh38/gencode.v33.long_noncoding_RNAs.txt", header=F, sep="",stringsAsFactors=FALSE)
gencode_protein_coding<-read.delim("STAR_GRCh38/gencode.v33.protein.coding.txt",sep="", header=F, stringsAsFactors=FALSE)
gencode_both<-rbind(gencode_protein_coding,gencode_lncrna)

######################
#### design matrix
######################
## design matirx has a sample ID column and condition column
design_matrix<-read.csv("design_matrix.csv")
# I am also going to add a column for the first contrast where HF- and HF+ are one group to compare with the HVOLs
# the design matrix 2 now has 3 columns 'sample_id','condition' and 'condition_contrast1'
# The groups are A = HVOL, B = CDCS HF- , C = CDCS HF+
condition_contrast1 = c(rep("B",time=58), rep("A",time=30), rep("B", time = 1))

###########################################################
##### plot a PCA first
##### rlog transform the data to plot PCA (vst is similar - much faster)
##### 'Working with rlog-transformed data' 'https://bioc.ism.ac.jp/packages/2.14/bioc/vignettes/DESeq2/inst/doc/beginner.pdf'
##########################################################

vst<-varianceStabilizingTransformation(dds)
for_pca<-as.data.frame(assay(vst))
pca_data_100<-prcomp(t(for_pca))
pca_data_perc_100=round(100*pca_data_100$sdev^2/sum(pca_data_100$sdev^2),1)
df_pca_data_100=data.frame(PC1 = pca_data_100$x[,1], PC2 = pca_data_100$x[,2], 
                           sample = colnames(for_pca), condition=design_matrix$condition)

png("ggplot_PCA_plasma_genes_PC1_PC2.png")
ggplot(df_pca_data_100,  aes(PC1,PC2, color = condition)) +
  theme(legend.position="bottom")+
  geom_point(size=2)+
  labs(x=paste0("PC1 (",pca_data_perc_100[1],")"),
       y=paste0("PC2 (",pca_data_perc_100[2],")"))+
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(aes_string(label=design_matrix$sample_id),size=3.5, nudge_x=-0.4, nudge_y=-30,  show.legend = FALSE)
dev.off()


pca_data_perc_100=round(100*pca_data_100$sdev^2/sum(pca_data_100$sdev^2),1)
df_pca_data_100=data.frame(PC2 = pca_data_100$x[,2], PC3 = pca_data_100$x[,3], 
                           sample = colnames(for_pca), condition=design_matrix$condition)

png("ggplot_PCA_plasma_genes_PC2_PC3.png")
ggplot(df_pca_data_100,  aes(PC2,PC3, color = condition)) +
  theme(legend.position="bottom")+
  geom_point(size=2)+
  labs(x=paste0("PC2 (",pca_data_perc_100[2],")"),
       y=paste0("PC3 (",pca_data_perc_100[3],")"))+
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(aes_string(label=design_matrix$sample_id),size=3.5, nudge_x=-0.4, nudge_y=-30,  show.legend = FALSE)
dev.off()

#########################################################################
##### Discard outliers (from the PCA analysis)
#########################################################################

## I'm getting rid of these samples
to_discard<-c("S25","S43","S76")
design_matrix2<-subset(design_matrix, !design_matrix$sample_id %in% to_discard)
design_matrix2<-cbind(design_matrix2, condition_contrast1)

###########################################################
##### First I do an ANOVA on sequencing depth
#### I exported the # reads aligned from multiqc into a .csv file
##########################################################

read_depth<-read.csv("read_depth.csv")
depth_ANO<-merge(design_matrix2,read_depth, by="sample_id")
depth_ANO$group<-ifelse(depth_ANO$condition=="A", "HVOL", ifelse(depth_ANO$condition=="B","CDCS-","CDCS+"))
## attach an object so we don't have to use read_depth$condition to access column - we just use condition
attach(depth_ANO)
par(mar=c(5,3,3,2))
boxplot(M.Aligned~group)
ANOVA<-aov(M.Aligned~group)
summary(ANOVA)
# this gives the 95% CI of the 3 groups
plot(TukeyHSD(ANOVA), las=1)

png("Read depths of the three groups.png")
ggplot(depth_ANO,aes(x=group, y=M.Aligned, fill=group)) + geom_boxplot()
dev.off()


###########################################################
##### I import salmon counts with tximport
##### Import the trancripts for novels
##########################################################

# I need to create a new DESeqDataSet using this design matrix to do this one off contrast
files.tx <- file.path("analysis/results/05b_salmon_counts", design_matrix2$sample_id, "quant.sf")
names(files.tx)<-design_matrix2$sample_id
all(file.exists(files.tx))

##import the quant.sf files for transcript level
txi.tx<-tximport(files.tx,type="salmon", txOut=TRUE, countsFromAbundance = 'no')
dds.tx <- DESeqDataSetFromTximport(txi.tx, design_matrix2, ~condition_contrast1)
idx<- (counts(dds.tx)[apply(counts(dds.tx),1,function(x) sum(x>0)/length(x))>0.5,])
dds.tx <- dds.tx[rownames(idx),]
dds.tx<- DESeq(dds.tx, parallel=T)
############################################################
###### novels 
############################################################
# get the TPMS for the txs
tx_table_TPMs<-as.data.frame(txi.tx$abundance)
tx_table_TPMs <- tx_table_TPMs[rownames(idx),]
##read the lst of potential novels from my pipeline in
novel <-  read.delim("analysis/results/04_novel/novel.cpat.transcript.list", header=F, stringsAsFactors=FALSE)
colnames(novel)<-c("Chromosome","Start","Stop","Strand", "MSTRG_ID","class_code","gene_name")
novel<-merge(novel[,-7], tx_table_TPMs, by.x="MSTRG_ID",by.y='row.names')
##merge novel with class_code_transcriptlist so we can get the exon number info
novel<-merge(novel[,c(1:5)],class_code_transcriptlist, 
             by.x="MSTRG_ID",by.y="transcript_id", all.x=TRUE)
novel<-novel[which(novel$number_of_exons>1),]
novel<-novel[which(novel$class_code=="u"),]
novel_all_TPMs<-merge(novel,tx_table_TPMs,by.x="MSTRG_ID",by.y='row.names')
#work out some stats on the novels
library(matrixStats)
novel$no.SamplesNull<-rowSums(novel_all_TPMs[,-c(1:8)] == 0)
novel$no.SamplesNotnull<-rowSums(novel_all_TPMs[,-c(1:8)]>0)
novel$no.SamplesGE1<-rowSums(novel_all_TPMs[,-c(1:8)]>1)
novel$max_TPM<-round(apply(novel_all_TPMs[,-c(1:8)],1,max))
novel$min_TPM<-round(apply(novel_all_TPMs[,-c(1:8)],1,min))
novel$mean_TPM<-round(apply(novel_all_TPMs[,-c(1:8)],1,mean))
novel$std_TPM<-round(rowSds(as.matrix(novel_all_TPMs[,-c(1:8)])),digits=2)
quants<-rowQuantiles(as.matrix(novel_all_TPMs[,-c(1:8)], probs = 0.25,0.5,0.75))
novel$quant25<-round(quants[,2],digits=2)
novel$quant50<-round(quants[,3],digits=2)
novel$quant75<-round(quants[,4],digits=2)

write.csv( as.data.frame(novel), file="/media/chi/novel.csv", row.names=TRUE)


###########################################################
###########################################################
##### gene level on  gencode 
##########################################################
###########################################################

# I created this file by taking the ensembl gene and ensembl tx id from GENCODE
my_tx2gene<-read.delim("STAR_GRCh38/gencode.v33.tx2gene.txt", sep="",header=F, stringsAsFactors=FALSE)
files <- file.path("analysis/results/05b_salmon_gene_counts", design_matrix2$sample_id, "quant.sf")
names(files)<-design_matrix2$sample_id
all(file.exists(files))

##import the quant.sf files for transcript level
my.txi<-tximport(files, type="salmon", tx2gene= my_tx2gene)

my.dds <- DESeqDataSetFromTximport(my.txi, design_matrix2,  ~condition_contrast1)
# I want to filter on TPMs before we run the DESEQ anlaysis (so we're not penalised for multiple testing for really low abundant genes - 
# e.g. genes that could be noise. From the test (coded below) we chose a cut-off of 90% sample having TPM >1)
gene_table_TPMs<-as.data.frame(my.txi$abundance)
idx<-gene_table_TPMs[apply(gene_table_TPMs , 1, function(x) sum(x>1)/(length(x))>=0.9),]
gene_table_1TPM <- gene_table_TPMs[rownames(idx),]
my.dds_1TPM<-my.dds[rownames(idx),]
my.dds_1TPM<- DESeq(my.dds_1TPM, parallel=T)

############################################
##### B & C vs A (= CDCS vs HVOL)
############################################
my.res_con_BC_A<-results(my.dds_1TPM, contrast=c("condition_contrast1","B","A"),alpha = 0.05)
my.res_con_BC_A_df<-as.data.frame(my.res_con_BC_A)
my.res_con_BC_A_df<-merge(gencode_both[,c(4:6)], my.res_con_BC_A_df, by.x='V4',by.y='row.names')
names(my.res_con_BC_A_df)[1:3]<-c("ensembl_id","gene_name","gene_type")
# lets take the ones where logfc increases with reards to HVOL (A)
my.res_con_BC_A_df<-my.res_con_BC_A_df[which(my.res_con_BC_A_df$log2FoldChange>0),]
my.res_con_BC_A_df_005<-my.res_con_BC_A_df[which(my.res_con_BC_A_df$padj<0.05 ),]
my.res_con_BC_A_df_095<-my.res_con_BC_A_df[which(my.res_con_BC_A_df$padj>=0.05 ),]
#this orders with regards to padj (lowest padj first)
my.res_con_BC_A_df_005<-my.res_con_BC_A_df_005[order(my.res_con_BC_A_df_005$padj),]
ens_id<-my.res_con_BC_A_df_005$ensembl_id
gene_id<-my.res_con_BC_A_df_005$gene_name
# The groups are A = HVOL, B = CDCS HF- , C = CDCS HF+
cond = as.factor(c(rep("HF-",time=29), rep("HF+",time=29),rep("HVOL",time=30), rep("HF-", time = 1)))

## draw the plots for all 134 padj < 0.05

for(i in 1:length(ens_id)){
  ind = match(ens_id[i],row.names(gene_table_TPMs))
  test<-as.data.frame(gene_table_TPMs[ind,])
  test=as.data.frame(t(test))
  test$cond=factor(cond, levels = c("HVOL", "HF-", "HF+"))
  p=ggplot(test, aes_string(x="cond", y=ens_id[i], fill="cond")) +  xlab("cond") + ylab( gene_id[i])+ 
    geom_boxplot()
  ggsave(p, file=paste0("/media/chi/Figures/HVOLvsCDCS/padj_005",gene_id[i],".png"), width = 14, height = 10, units = "cm")
}

# Subset the ones with low padj and label the outliying samples
### https://stackoverflow.com/questions/62906687/labeling-outliers-of-boxplots-in-a-loop-r/62907443#62907443
my.res_con_BC_A_df_0001<-my.res_con_BC_A_df_005[which(my.res_con_BC_A_df_005$padj<0.001),]
ens_id<-my.res_con_BC_A_df_0001$ensembl_id
gene_id<-my.res_con_BC_A_df_0001$gene_name

library(tibble)
is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}

for(i in 1:length(ens_id)){
  ind = match(ens_id[i],row.names(gene_table_TPMs))
  test<-as.data.frame(gene_table_TPMs[ind,])
  test=as.data.frame(t(test))
  test$cond=factor(cond, levels = c("HVOL", "HF-", "HF+"))
  dat <- test %>% tibble::rownames_to_column(var="outlier") %>% group_by(cond) %>%  mutate_at(vars(ens_id[i]), list(is_outlier = ~replace(., !is_outlier(.), NA))) 
  #mutate(across(ens_id[i], ~ replace(., !is_outlier(.), NA)))
  dat$outlier[which(is.na(dat$is_outlier))] <- as.numeric(NA)
  #print(head(dat))}
  p=ggplot(dat, aes_string(y=ens_id[i], x="cond",fill="cond")) + geom_boxplot()  + ylab(gene_id[i])+ geom_text(aes(label=outlier),na.rm=TRUE,nudge_x=0.15)
  ggsave(p, file=paste0("/media/chi/Figures/HVOLvsCDCS/padj_0001/",gene_id[i],".png"), width = 14, height = 10, units = "cm")
}


#####################################################################
##### caluclate library size from linear to use for the circs 
##### then circ analysis
#####################################################################
# we can mormalise the circs by using the library size from the genes
# figure 2 https://www.nature.com/articles/s41467-019-13840-9


## if tximport is used then we can't use sizeFactors(). We have to use normalizationFactors() instead
## This gives a factor for each gene in each sample. I then took a median of each column 
## I compared this with the output from estimateSizeFactorsForMatrix using counts(my.dds) as the input matrix
## and they were very similar (after comparing library_size_counts_dds and library_size_for_circs)
library_size_counts_dds = estimateSizeFactorsForMatrix(counts(my.dds_1TPM))

test3<-as.data.frame(normalizationFactors(my.dds_1TPM))
library_size_for_circs<-apply(test3,2,FUN=median);rm(test3)

# first import the circs
sample_names=as.character(design_matrix2$sample_id)
setwd("/HANGER/zoe/zoe_HiSeq_OG4071/analysis/results/07_circs/")
## initialise the list to store the tables in:
list.data<-list()

## I had to put the tryCatch in as for the MiSeq data some files were empty and it was crashing
for(i in sample_names){
  #print(i)
  #print(read.table(paste0(i,'/circularRNA_known.bed')))
  list.data[[i]] <- tryCatch(read.table(paste0(i,'/circularRNA_known.bed'), header=FALSE, sep="\t",stringsAsFactors=FALSE),error=function(e) NULL)
  list.data[[i]] <-list.data[[i]][,c(1:3,13)]
}


## import the file with all of the circs ids from the samples combined
all_circs <- read.table("all_circular_sorted.bed", header=FALSE, sep="\t",stringsAsFactors=FALSE)
#take the chrme, strand, gene name gene id
all_circs<- all_circs[,c(1:3,6,15,16)]
# start the merge of this table with the first sample circs
final_circs <- left_join(all_circs, list.data[[1]], by =c("V1", "V2","V3"))
## loop through the rest

for(i in 2:length(list.data)){
  #print(sample_names[i])
  final_circs <- left_join(final_circs, list.data[[i]], by =c("V1", "V2","V3"))
}
## replace nas with 0
final_circs[is.na(final_circs)] <- 0
colnames(final_circs) <- c("chr", "start", "stop","strand","gene_symbol","Transcript_ID",names(list.data))

final_circs<- final_circs[apply(final_circs[,c(7:ncol(final_circs))],1,function(x) sum(x>0)/length(x))>0.9,]

# prepare df for DESEQ I will keep the rownames and ids to combine later
meta_circs<-final_circs[,c(1:6)]
final_circs_DESEQ<-final_circs[,-c(1:6)]

# To get normalised counts of the circs
# now we have to apply the library size from the genes calculated above to the df 
# going to use mapply https://stackoverflow.com/questions/18382883/what-is-the-right-way-to-multiply-data-frame-by-vector

final_circs2<-data.frame(mapply(`/`,final_circs[,c(7:81)],library_size_for_circs))
# with mapply we've lost the data that isn't involved in the multiplication eg tx ids etc
# I'm going to do a dirty cbind but then check that the tx ids are in the same order as the original
final_circs3<-cbind(final_circs[,c(1:6)], final_circs2)
table(final_circs3$Transcript_ID==final_circs$Transcript_ID)
#final_circs<-final_circs3; rm(final_circs2,final_circs3)
final_circs2.names<-names(final_circs2)[(names(final_circs2) %in% filtered_names)]
final_circs2<-final_circs2[,final_circs2.names]

# log the data I checked that the data was normalised by eyeballing some histograms
# ## first replace any zeros with 1s (we can't log zero)
final_circs2[final_circs2 == 0]<-1
#circtools.comb.5M.normalised<-merge(circtools.comb_meta,circtools.comb.5M.normalised,by='row.names')
# # natural log the df
ln.final_circs2<-as.data.frame(apply(final_circs2, 2, function(x) log(x)))


#now get the column numbers of the HVOLS and CDCS groups to run the T=Test
HVOL.5M<-c(50:74)
CDCS.5M<-c(1:49,75)

pvals <- rep(NA, nrow(ln.final_circs2))
for(i in 1:nrow(ln.final_circs2)) pvals[i] <- t.test(ln.final_circs2[i,c(HVOL.5M)],ln.final_circs2[i,c(CDCS.5M)])$p.value
ln.final_circs2<-cbind(pvals,ln.final_circs2)
row.names(ln.final_circs2)<-row.names(final_circs)
ln.final_circs2<-merge(meta_circs,ln.final_circs2,by='row.names')
row.names(ln.final_circs2)<-ln.final_circs2$Row.names;ln.final_circs2$Row.names<-NULL


padj<-p.adjust(ln.final_circs2$pvals,method="BH")
ln.final_circs2<-cbind(ln.final_circs2[,c(1:7)],padj,ln.final_circs2[,c(8:82)])
ln.final_circs2.0001<-ln.final_circs2[which(ln.final_circs2$padj<0.001),]

# overlap between circexplorer and circtools (Circtools seemed to be one off in the start coordinate (it's a bed))
ln.circtools.comb.5M.normalised.0001_annot$Start<-ln.circtools.comb.5M.normalised.0001_annot$Start-1
overlap<-merge(ln.final_circs2.0001[,c(1:8)],ln.circtools.comb.5M.normalised.0001_annot[,c(1:10)],by.x=c("chr","start","stop","strand"),by.y=c("Chr","Start","End","Strand"))

cn=row.names(ln.final_circs2.0001)
for(i in 1:length(cn)){
  ind = match(cn[i],row.names(ln.final_circs2.0001))
  test_circ2<-as.data.frame(ln.final_circs2.0001[ind,c(9:ncol(ln.final_circs2.0001))])
  test_circ2=as.data.frame(t(test_circ2))
  test_circ2$cond=as.factor(c(rep("CDCS",49),rep("HVOL", 25),rep("CDCS",1)))
  #print(tail(test_circ2))}
  test_circ2$cond=factor(cond, levels = c("HVOL", "CDCS"))
  #    print(ln.final_circs2$gene_symbol[i])}
  p=ggplot(test_circ2, aes_string(x="cond", y=as.name(cn[i]), fill="cond")) +  xlab("cond") + ylab( ln.final_circs2$gene_symbol[i])+
    geom_boxplot()
  ggsave(p, file=paste0("/media/chi/Figures/circexplorer.5M/",ln.final_circs2$gene_symbol[i],".png"), width = 14, height = 10, units = "cm")
  
}

################################################################################
################ test of varying tpms as input to filter  ######################
################################################################################
################################################################################

gene_table_TPMs<-as.data.frame(txi$abundance)
## I'm going to only take the genes that are not 'MSTRG'
idx<- gene_table_TPMs[!grepl("MSTRG.",row.names(gene_table_TPMs)),]
gene_table_TPMs <- gene_table_TPMs[rownames(idx),]


tpm<-c(0,0.025,0.05,0.075,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2.5,5,7.5,10)
no.genes<-numeric()
for (i in 1:length(tpm)){
  idx<-gene_table_TPMs[apply(gene_table_TPMs , 1, function(x) sum(x>tpm[i])/(length(x))>=0.9),]
  no.genes[i] = nrow(idx)
}

# https://stackoverflow.com/questions/4877357/how-to-plot-all-the-columns-of-a-data-frame-in-r
require(ggplot2)
require(reshape2)

df<-data.frame(tpm,no.genes)
df <- melt(df ,  id.vars = 'tpm', variable.name = 'number_of_genes')


png("Varying 0-10 TPM filter.png")
ggplot(df, aes(tpm,value)) + geom_line(aes())  +   
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9,10))

dev.off()

###########################################################
##### Plot TissueEnrich
##########################################################
## I want to get a normalised table of the raw counts for TissueEnrich
dds<- estimateSizeFactors(dds)
norm_data_gene<- as.data.frame(round(counts(dds, normalized=TRUE)))
# I remove the 'MSTRG' genes - we only want ENSG
norm_data_gene<-norm_data_gene[!grepl("MSTRG.",row.names(norm_data_gene)),]
tiss_enr<-merge(gencode_both[,c(4:5)], norm_data_gene, by.x='V4',by.y='row.names')
colnames(tiss_enr)[1:2]<-c("ensembl_id","gene_name")
sample_names<-colnames(norm_data_gene)

#sample_names = colnames(norm_data_gene)[52]
#sample_names_top_100<-tiss_enr[,names(tiss_enr) %in% sample_names[i],drop=F]
for(i in 1:length(sample_names)){
  sample_names_top_100_1<-tiss_enr[,1:2]
  sample_names_top_100<-tiss_enr[,names(tiss_enr) %in% sample_names[i],drop=F]
  sample_names_top_100<- cbind(sample_names_top_100_1,sample_names_top_100);
  #sample_names_top_100<-sample_names_top_100[!grepl("MT-",sample_names_top_100$gene_name),]
  sample_names_top_100<-sample_names_top_100[order(-sample_names_top_100[3]),]
  sample_names_top_100<-sample_names_top_100[1:100,]
  
  # for outlier
  gs<-GeneSet(geneIds=sample_names_top_100$gene_name,organism="Homo Sapiens",geneIdType=SymbolIdentifier())
  output<-teEnrichment(inputGenes = gs)
  
  seEnrichmentOutput<-output[[1]]
  enrichmentOutput<-setNames(data.frame(assay(seEnrichmentOutput),row.names = rowData(seEnrichmentOutput)[,1]), colData(seEnrichmentOutput)[,1])
  enrichmentOutput$Tissue<-row.names(enrichmentOutput)
  #head(enrichmentOutput)
  
  # png(paste("Figures/tissue_enrichment_",sample_names[i],"_removed.png"))
  p = ggplot(enrichmentOutput,aes(x=reorder(Tissue,-Log10PValue),y=Log10PValue,label = Tissue.Specific.Genes,fill = Tissue))+
    geom_bar(stat = 'identity')+
    labs(title=sample_names[i],x='', y = '-LOG10(P-Adjusted)')+
    theme_bw()+
    theme(legend.position="none")+
    theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),panel.grid.major= element_blank(),panel.grid.minor = element_blank())
  
  ggsave(p, file=paste0("Figures/tissue_enrichment_",sample_names[i],".png"), width = 14, height = 10, units = "cm")
}
