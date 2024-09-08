library(dplyr)

#This code is to obtain index numbers to be used for Hi-C analysis. 
#This is the code for genes downregulated after 2h of rapamycin treatment which are predicted to be Fkh1 bound.
#sig_down_2 and TF files are changed accordingly.

# Rap 2h Downregulated 
sig_down_2 = read.table("/Users/audreypeng/Desktop/RNAseq_analysis/deseq2/downreg_2h.txt", sep="\t", header = TRUE)

# list of genes from YEASTract+
TF = read.table("down2_fkh1.txt", header = FALSE, sep="")

# Adding index numbers to be read into Hi-C analysis files
TF_df = as.data.frame(t(TF))
colnames(TF_df) = "genename"
TF_fulldf = inner_join(TF_df,sig_down_2)
TF_fulldf = TF_fulldf[order(TF_fulldf$log2FoldChange),]
TF_fulldf = TF_fulldf[,c(1,2,4,8:11)]
TF_fulldf$mid = ave(TF_fulldf$start_position, TF_fulldf$end_position)
TF_fulldf$idx = trunc(c(TF_fulldf$mid/5000))
TF_idx = TF_fulldf[,c(1,2,3,5,9)]
colnames(TF_idx)[2] ="ensembl"
write.table(TF_idx, file = "/Users/audreypeng/Desktop/RNAseq_analysis/deseq2/down2_fkh1_idx.txt", sep="\t", quote=F, row.names=F)
