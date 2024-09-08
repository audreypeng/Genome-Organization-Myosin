setwd("/Users/audreypeng/Desktop/RNAseq_analysis/deseq2")
ds2h = read.delim('deseq_2h.txt', header = T, sep = '\t')
ds5h = read.delim('deseq_5h.txt', header = T, sep = '\t')

ds2h = ds2h[,c(1,3,7:11)]
ds2h = ds2h[order(ds2h$chromosome_name),]

ds2h_nomito = subset(ds2h, !(ds2h$chromosome_name == 'Mito'))
ds2h_nomito$chromosome_name = as.numeric(as.roman(ds2h_nomito$chromosome_name))

ds2h_nomito$mid = ave(ds2h_nomito$start_position, ds2h_nomito$end_position)
ds2h_nomito$idx = trunc(c(ds2h_nomito$mid/5000))
ds2h_nomito = ds2h_nomito[,c(1:5,9)]
colnames(ds2h_nomito)[1] ="genesymbol"

write.table(ds2h_nomito, file = "/Users/audreypeng/Desktop/RNAseq_analysis/deseq2/ds2h_idx.txt", sep="\t", quote=F, row.names=F)
write.table(ds2h_nomito, file = "/Users/audreypeng/Desktop/HiC_analysis/ds2h_idx.txt", sep="\t", quote=F, row.names=F)


ds5h = ds5h[,c(1,3,7:11)]
ds5h = ds5h[order(ds5h$chromosome_name),]

ds5h_nomito = subset(ds5h, !(ds5h$chromosome_name == 'Mito'))
ds5h_nomito$chromosome_name = as.numeric(as.roman(ds5h_nomito$chromosome_name))

ds5h_nomito$mid = ave(ds5h_nomito$start_position, ds5h_nomito$end_position)
ds5h_nomito$idx = trunc(c(ds5h_nomito$mid/5000))
ds5h_nomito = ds5h_nomito[,c(1:5,9)]
colnames(ds5h_nomito)[1] ="genesymbol"

write.table(ds5h_nomito, file = "/Users/audreypeng/Desktop/RNAseq_analysis/deseq2/ds5h_idx.txt", sep="\t", quote=F, row.names=F)
write.table(ds5h_nomito, file = "/Users/audreypeng/Desktop/HiC_analysis/ds5h_idx.txt", sep="\t", quote=F, row.names=F)

