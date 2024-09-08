library("ggplot2") # for the plot
library("ggrepel") # for spreading text labels on the plot, you can replace with `geom_text` if you want
library("scales") # for axis labels notation

# chromosome sizes and centromere locations retrieved from YeastMine version updated on 28 Feb 2023.

# data with sig up or down genes based on L2FC < 1 or < -1, padj<0.01 using Deseq2
sample_down <- structure(list(gene = sig_down_2$genename, 
                             chromosome = as.character(sig_down_2$chromosome_name), 
                             start = sig_down_2$start_position,
                             end = sig_down_2$end_position))
sample_down = as.data.frame(sample_down)                        

sample_up <- structure(list(gene = sig_up_2$genename, 
                             chromosome = as.character(sig_up_2$chromosome_name), 
                             start = sig_up_2$start_position,
                             end = sig_up_2$end_position))
sample_up = as.data.frame(sample_up)                         

# R64-1-1 chromosome sizes
chrom_sizes <- structure(list(chromosome = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16"),
                              size = c(230218,813184,316620,1531933,576874,270161,1090940,562643,439888,745751,666816,1078177,924431,784333,1091291,948066)))
chrom_sizes = as.data.frame(chrom_sizes)

#R64-1-1 centromere locations
centromeres <- structure(list(chromosome = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16"),
                              start = c(151465, 238207, 114385, 449711, 151987, 148510, 496920, 105586, 355629, 436307, 440129, 150828, 268031, 628758, 326584, 555957),
                              end = c(151582, 238323, 114501, 449821,152104, 148627,497038, 105703,355745,436425,440246, 150947, 268149, 628875,326702, 556073)))
centromeres = as.data.frame(centromeres)


## DOWN
pdf("/Users/audreypeng/Desktop/RNAseq_analysis/deseq2/down_2h_chrplot.pdf")
ggplot(data = chrom_sizes) + 
  # base rectangles for the chroms, with numeric value for each chrom on the x-axis
  geom_rect(aes(xmin = as.numeric(chromosome) - 0.2, 
                xmax = as.numeric(chromosome) + 0.2, 
                ymax = size, ymin = 0), 
            colour="black", fill = "white") + 
  # rotate the plot 90 degrees
  coord_flip() +
  # black & white color theme 
  theme(axis.text.x = element_text(colour = "black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank()) + 
  # give the appearance of a discrete axis with chrom labels
  scale_x_discrete(name = "chromosome", limits = names(chrom_sizes$chromosome)) +
  # add bands for centromeres
  geom_rect(data = centromeres, aes(xmin = as.numeric(chromosome) - 0.2, 
                                    xmax = as.numeric(chromosome) + 0.2, 
                                    ymax = end, ymin = start), colour="yellow") +
  # add bands for CNA value
  geom_rect(data = sample_down, aes(xmin = as.numeric(chromosome) - 0.2, 
                                  xmax = as.numeric(chromosome) + 0.2, 
                                   ymax = end, ymin = start), colour="blue") + 

  #geom_rect(data = sample_up, aes(xmin = as.numeric(chromosome) - 0.2, 
   #                                xmax = as.numeric(chromosome) + 0.2, 
    #                               ymax = end, ymin = start), colour="orange") +
  # supress scientific notation on the y-axis
  scale_y_continuous(labels = comma) +
  ylab("region (bp)")
  dev.off()                       


### UP
pdf("/Users/audreypeng/Desktop/RNAseq_analysis/deseq2/up_2h_chrplot.pdf")
ggplot(data = chrom_sizes) + 
  # base rectangles for the chroms, with numeric value for each chrom on the x-axis
  geom_rect(aes(xmin = as.numeric(chromosome) - 0.2, 
                xmax = as.numeric(chromosome) + 0.2, 
                ymax = size, ymin = 0), 
            colour="black", fill = "white") + 
  # rotate the plot 90 degrees
  coord_flip() +
  # black & white color theme 
  #theme(axis.text.x = element_text(colour = "black"), 
        #panel.grid.major = element_blank(), 
        #panel.grid.minor = element_blank(), 
        #panel.background = element_blank()) + 
  # give the appearance of a discrete axis with chrom labels
  scale_x_discrete(name = "chromosome", limits = names(chrom_sizes$chromosome)) +
  # add bands for centromeres
  geom_rect(data = centromeres, aes(xmin = as.numeric(chromosome) - 0.2, 
                                    xmax = as.numeric(chromosome) + 0.2, 
                                    ymax = end, ymin = start), colour="yellow") +
  
  geom_rect(data = sample_up, aes(xmin = as.numeric(chromosome) - 0.2, 
                                  xmax = as.numeric(chromosome) + 0.2, 
                                 ymax = end, ymin = start), colour="red") +
  # supress scientific notation on the y-axis
  
  scale_y_continuous(labels = comma) +
  ylab("region (bp)")
dev.off()

###UP and DOWN
pdf("/Users/audreypeng/Desktop/RNAseq_analysis/deseq2/both_2h_chrplot.pdf")
ggplot(data = chrom_sizes) + 
  # base rectangles for the chroms, with numeric value for each chrom on the x-axis
  geom_rect(aes(xmin = as.numeric(chromosome) - 0.2, 
                xmax = as.numeric(chromosome) + 0.2, 
                ymax = size, ymin = 0), 
            colour="black", fill = "white") + 
  # rotate the plot 90 degrees
  coord_flip() +
  # black & white color theme 
  #theme(axis.text.x = element_text(colour = "black"), 
  #panel.grid.major = element_blank(), 
  #panel.grid.minor = element_blank(), 
  #panel.background = element_blank()) + 
  # give the appearance of a discrete axis with chrom labels
  scale_x_discrete(name = "chromosome", limits = names(chrom_sizes$chromosome)) +
  # add bands for centromeres
  geom_rect(data = centromeres, aes(xmin = as.numeric(chromosome) - 0.2, 
                                    xmax = as.numeric(chromosome) + 0.2, 
                                    ymax = end, ymin = start), colour="yellow") +
  # add bands for CNA value
  geom_rect(data = sample_down, aes(xmin = as.numeric(chromosome) - 0.2, 
                                    xmax = as.numeric(chromosome) + 0.2, 
                                    ymax = end, ymin = start), colour="blue") + 
  
  geom_rect(data = sample_up, aes(xmin = as.numeric(chromosome) - 0.2, 
                                  xmax = as.numeric(chromosome) + 0.2, 
                                  ymax = end, ymin = start), colour="red") +
  # supress scientific notation on the y-axis
  
  scale_y_continuous(labels = comma) +
  ylab("region (bp)")
dev.off()
