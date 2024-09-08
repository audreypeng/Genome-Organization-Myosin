library("DESeq2")
library("AnnotationDbi")
library("org.Sc.sgd.db")
library(biomartr)
library(EnhancedVolcano)
library(qqman)

directory <- "/Users/audreypeng/Desktop/RNAseq_analysis/deseq2/counts"
sampleFiles <- grep("AA",list.files(directory),value=TRUE)
sampleName <- c('DMSO_1', 'DMSO_2', 'AA2_1', 'AA2_2', 'AA5_1', 'AA5_2')
condition <- c('DMSO','DMSO','2h','2h','5h', '5h')

sampleTable <- data.frame(sampleName = sampleName,
                          fileName = sampleFiles,
                          condition = condition)
sampleTable
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~ condition)
ddsHTSeq$condition <- factor(ddsHTSeq$condition, levels = c("DMSO","2h", "5h"))

## Pre-filtering
nrow(ddsHTSeq)
keep <- rowSums(counts(ddsHTSeq)) > 50
dds <- ddsHTSeq[keep,]

## estimate dispersions and fit to NBM with LRT test
dds <- DESeq(dds, test = "LRT", reduced = ~1)

res2 = results(dds, contrast=c("condition", "2h","DMSO"))
res5 = results(dds, contrast=c("condition", "5h","DMSO"))
res2v5 = results(dds, contrast=c("condition", "5h","2h"))


### getting gene names
ens.str <- substr(row.names(res2), 1, 15)
res2$genename <- mapIds(org.Sc.sgd.db,
                     keys=ens.str,
                     column="GENENAME",
                     keytype="ENSEMBL",
                     multiVals="first")

ens.str <- substr(row.names(res5), 1, 15)
res5$genename <- mapIds(org.Sc.sgd.db,
                        keys=ens.str,
                        column="GENENAME",
                        keytype="ENSEMBL",
                        multiVals="first")

### Put chr positions in for res2
# Finding the ensembl
biomartr::getMarts()
# Finding datasets
print(biomartr::getDatasets(mart = "fungi_mart"), n=70)
# What attributes it has
head( biomartr::getAttributes(mart    = "fungi_mart", 
                              dataset = "scerevisiae_eg_gene"), 10 )

gene_set = row.names(res2)
result_BM <- biomartr::biomart( genes      = gene_set, # genes were retrieved using biomartr::getGenome()
                                mart       = "fungi_mart", # marts were selected with biomartr::getMarts()
                                dataset    = "scerevisiae_eg_gene", # datasets were selected with biomartr::getDatasets()
                                attributes = c("chromosome_name", "start_position", "end_position", "strand"), 
                                filters    = "ensembl_gene_id") # specify what ID type was stored in the fasta file retrieved with biomartr::getGenome()
row.names(result_BM) = result_BM[,1]
result_BM = result_BM[,2:5]

res2df <- as.data.frame(res2)
res2df = merge(res2df, result_BM, by="row.names", all=T)
########

### now for res5
gene_set = row.names(res5)
result_BM <- biomartr::biomart( genes      = gene_set, # genes were retrieved using biomartr::getGenome()
                                mart       = "fungi_mart", # marts were selected with biomartr::getMarts()
                                dataset    = "scerevisiae_eg_gene", # datasets were selected with biomartr::getDatasets()
                                attributes = c("chromosome_name", "start_position", "end_position", "strand"), 
                                filters    = "ensembl_gene_id") # specify what ID type was stored in the fasta file retrieved with biomartr::getGenome()
row.names(result_BM) = result_BM[,1]
result_BM = result_BM[,2:5]

res5df <- as.data.frame(res5)
res5df = merge(res5df, result_BM, by="row.names", all=T)
####

##save deseq output
write.table(res2df, file = "/Users/audreypeng/Desktop/RNAseq_analysis/deseq2/deseq_2h.txt", sep="\t", quote=F, row.names=F)
write.table(res5df, file = "/Users/audreypeng/Desktop/RNAseq_analysis/deseq2/deseq_5h.txt", sep="\t", quote=F, row.names=F)
###

## Volcano plot
tiff("/Users/audreypeng/Desktop/RNAseq_analysis/deseq2/deseq2_volcanoplot_2h.tiff", res = 300, width=2300, height=2000)
EnhancedVolcano(res2df,
                lab = res2df$genename,
                x = 'log2FoldChange',
                y = 'padj')
dev.off()

tiff("/Users/audreypeng/Desktop/RNAseq_analysis/deseq2/deseq2_volcanoplot_5h.tiff", res = 300, width=2300, height=2000)
EnhancedVolcano(res5df,
                lab = res5df$genename,
                x = 'log2FoldChange',
                y = 'padj')
dev.off()



### set cutoff with padj (BH adjusted) <0.05 and L2FC to at least 2 fold changes (1 and -1)

sig_up_2 	= subset(res2df, res2df$padj < .01 & res2df$log2FoldChange > 1)
sig_down_2= subset(res2df, res2df$padj < .01 & res2df$log2FoldChange < -1)

#sig_nochange_2 = subset(res2df, res2df$padj < .01 & res2df$log2FoldChange > -0.1 & res2df$log2FoldChange < 0.1)
#sig_nochange_2 = sig_nochange_2[!is.na(sig_nochange_2$genename),]

sig_up_2$chromosome_name = as.numeric(as.roman(sig_up_2$chromosome_name))
sig_down_2$chromosome_name = as.numeric(as.roman(sig_down_2$chromosome_name))
#sig_nochange_2$chromosome_name = as.numeric(as.roman(sig_nochange_2$chromosome_name))

##checkcheck how many up and down?
sig_up_2[order(sig_up_2$padj),]
nrow(sig_up_2)
nrow(sig_down_2)
##

sig_up_5 	= subset(res5df, res5df$padj < .01 & res5df$log2FoldChange > 1)
sig_down_5= subset(res5df, res5df$padj < .01 & res5df$log2FoldChange < -1)

#sig_nochange_5 = subset(res5df, res5df$padj < .01 & res5df$log2FoldChange > -1 & res5df$log2FoldChange < 1)
#sig_nochange_5 = sig_nochange_5[!is.na(sig_nochange_5$genename),]

sig_up_5$chromosome_name = as.numeric(as.roman(sig_up_5$chromosome_name))
sig_down_5$chromosome_name = as.numeric(as.roman(sig_down_5$chromosome_name))

## to change chr names from roman numerals to numeric - not necessary but makes it easier for me
sig_up_5$chromosome_name = as.numeric(as.roman(sig_up_5$chromosome_name))
sig_down_5$chromosome_name = as.numeric(as.roman(sig_down_5$chromosome_name))

## tried to manhattan plot
manhattan(sig_up_2, chr="chromosome_name", bp="start_position", snp="genename", p="padj" )
## doesnt show centromeres and not clear enough either
## https://stackoverflow.com/questions/33727432/how-to-plot-positions-along-a-chromosome-graphic
## refer to chromosomeplot.R - modified version!


#save sig up and down files as txt
write.table(sig_up_2, file = "/Users/audreypeng/Desktop/RNAseq_analysis/deseq2/upreg_2h.txt", sep="\t", quote=F, row.names=F)
write.table(sig_down_2, file = "/Users/audreypeng/Desktop/RNAseq_analysis/deseq2/downreg_2h.txt", sep="\t", quote=F, row.names=F)
write.table(sig_up_5, file = "/Users/audreypeng/Desktop/RNAseq_analysis/deseq2/upreg_5h.txt", sep="\t", quote=F, row.names=F)
write.table(sig_down_5, file = "/Users/audreypeng/Desktop/RNAseq_analysis/deseq2/downreg_5h.txt", sep="\t", quote=F, row.names=F)


#save just gene names for GO
write.table(sig_up_2$Row.names, file = "/Users/audreypeng/Desktop/RNAseq_analysis/deseq2/upreg_2h_ensembl.txt", sep="\t", quote=F, row.names=F)
write.table(sig_down_2$Row.names, file = "/Users/audreypeng/Desktop/RNAseq_analysis/deseq2/downreg_2h_ensembl.txt", sep="\t", quote=F, row.names=F)
write.table(sig_up_5$Row.names, file = "/Users/audreypeng/Desktop/RNAseq_analysis/deseq2/upreg_5h_ensembl.txt", sep="\t", quote=F, row.names=F)
write.table(sig_down_5$Row.names, file = "/Users/audreypeng/Desktop/RNAseq_analysis/deseq2/downreg_5h_ensembl.txt", sep="\t", quote=F, row.names=F)
##############################################


sig_up_2[order(sig_up_2$chromosome_name),]

