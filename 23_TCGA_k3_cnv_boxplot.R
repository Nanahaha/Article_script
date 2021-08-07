## ------ boxplot of copy number ------ ##

library(stringr)
library(reshape2)
library(ggplot2)
library(ggpubr)  
library(RColorBrewer)  

path <- getwd()
cnv <- read.table(file = str_c(path,"TCGA-OV.masked_cnv.tsv"),
                  header = T, sep = "\t", fill = T, check.names = F)
cnv[1:4,]

cnv$cnv_stat <- NA
cnv$cnv_stat[which(cnv$value > 0.3)] <- "Amplification"
cnv$cnv_stat[which(cnv$value < -0.3)] <- "Deletion"
cnv$cnv_stat[which(cnv$value>=-0.3 & cnv$value<=0.3)] <- "non_cnv"

cnv_df <- dcast(data = cnv[,c(1,6)], formula = sample~cnv_stat)
colnames(cnv_df)[1] <- "Patients"

cnv_df2 <- cnv_df[-which(cnv_df$Patients == "TCGA-23-1023-01R"),]
cnv_df2$Patients <- sapply(cnv_df2$Patients, function(x) str_sub(x, end = -2))
which(duplicated(cnv_df2$Patients))

subclass <- read.table(file = str_c(path,"k3\\TCGA_consensus_k3_hclust_subclass2.txt"),
                       header = T, fill = T, check.names = F)
subclass$Class <- sapply(subclass$Class,function(x) str_c("C",x))
subclass[1:4,]

cnv_df2_class <- merge(subclass, cnv_df2, by = "Patients", all = F)  
cnv_df2_class[1:5,]


## amplification boxplot
my_comparisons <- list(c("C1","C2"),c("C1","C3"),c("C2","C3"))
pdf(file = str_c(path,"k3\\TCGA_k3_cnv_amplification_boxplot.pdf"),
    width = 3, height = 3.5)
p <- ggplot(cnv_df2_class, aes(x = Class, y = Amplification, fill = Class)) +
  geom_boxplot() +
  theme_classic() + 
  scale_fill_manual(values = brewer.pal(3,"Set2")) +
  labs(y = "Number of Amplification") +
  theme(axis.title.x = element_blank(),
        legend.position="none") +  
  stat_compare_means(comparisons = my_comparisons,
                     tip.length = 0.01,
                     label = "p.signif")
print(p)
dev.off()

## deletion boxplot
my_comparisons <- list(c("C1","C2"),c("C1","C3"),c("C2","C3"))
pdf(file = str_c(path,"k3\\TCGA_k3_cnv_deletion_boxplot.pdf"),
    width = 3, height = 3.5)
p <- ggplot(cnv_df2_class, aes(x = Class, y = Deletion, fill = Class)) +
  geom_boxplot() +
  theme_classic() +  
  scale_fill_manual(values = brewer.pal(3,"Set2")) +
  labs(y = "Number of Deletion") +
  theme(axis.title.x = element_blank(),
        legend.position="none") +  
  stat_compare_means(comparisons = my_comparisons,
                     tip.length = 0.01,
                     label = "p.signif")
print(p)
dev.off()

save.image(str_c(path,"23_TCGA_k3_cnv_boxplot.RData"))

