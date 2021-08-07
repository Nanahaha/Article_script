## ------ ssGSEA evaluate immune cell infiltration ------ ##

library(tidyverse)
library(GSVA)
library(pheatmap)
library(RColorBrewer) 
library(ggpubr)
library(reshape2) 
library(limma)  

path <- getwd()

data <- data.frame(fread(file = str_c(path,"TCGA-OV_HTSeq_TPM_01survnoNA_rowsamp_log2.txt"),
                         header = T, fill = T, check.names = F), check.names = F) 
data[1:4,1:10]

gene_set <- read.csv(file = str_c(path,"The gene set of a specific subpopulation of immune cells.csv"),
                     header = T, fill = T, check.names = F)[,1:2]
head(gene_set)
list <- split(as.matrix(gene_set)[,1], gene_set[,2])

rownames(data) <- data$Patients
expr_matrix <- t(data[,-c(1:8)])
expr_matrix[1:4,1:4]

gsva_matrix <- gsva(expr_matrix, 
                    list,
                    method = "ssgsea",
                    kcdf = "Gaussian",
                    abs.ranking = TRUE)
#Warning message:
#  In .filterFeatures(expr, method) :
#  251 genes with constant expression values throuhgout the samples.

gsva_matrix[1:4,1:4]
pheatmap(gsva_matrix)

gsva_df <- data.frame("cell_type" = rownames(gsva_matrix), gsva_matrix, check.names = F)
gsva_df[1:4,1:4]

write.table(gsva_df, file = str_c(path,"k3\\TCGA_k3_ssGSEA_immu_scores.txt"),
            sep = "\t", quote = F, row.names = F)  # row.names一定要为T，因为是细胞名


load(str_c(path,"18_TCGA_k3_ssGSEA_immu_score_plot.RData"))


##################################################################
## boxplot
plotDensities(gsva_matrix, legend = F)

data_rowsample <- t(gsva_matrix)
data_rowsample[1:4,1:4]
colnames(data_rowsample)

data_rowsample_df <- data.frame("Patients" = rownames(data_rowsample), 
                                data_rowsample, check.names = F)

subclass <- read.table(file = str_c(path,"k3\\TCGA_consensus_k3_hclust_subclass2.txt"),
                       header = T, fill = T)
subclass$Class <- sapply(subclass$Class,function(x) str_c("C",x))

scores_class <- merge(subclass, data_rowsample_df, by = "Patients", all = F)
scores_class[1:5,1:5]

write.table(scores_class, file = str_c(path,"k3\\TCGA_k3_ssGSEA_immu_scores_class.txt"),
            sep = "\t", quote = F, row.names = F)


mydata <- melt(scores_class[,-1], id = "Class")
mydata[1:4,]

min(mydata$value)
max(mydata$value)

dev.size()
pdf(file = str_c(path,"k3\\TCGA_k3_ssGSEA_immu_score_boxplot.pdf"),
    width = 8, height = 4)
p <- ggplot(mydata, aes(x = variable, y = value, fill = Class))+
  geom_boxplot()+
  theme_classic()+  
  scale_fill_manual(values = brewer.pal(5,"Set2"))+
  labs(y = "Immunity Signature")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
print(p)
dev.off()

## statistical analysis
scores_class2 <- scores_class
colnames(scores_class2)[-c(1:2)] <- sapply(1:(ncol(scores_class2)-2), function(x) str_c("cell_",x))
colnames(scores_class2)

my_comparisons<-list(c("C1","C2"),c("C1","C3"),c("C2","C3"))

for (i in 3:ncol(scores_class2)) {
  formu <- as.formula(str_c(colnames(scores_class2)[i],"~ Class"))
  compa <- compare_means(formu, data = scores_class2)
  write.table(compa,file = str_c(path,"k3\\ssGSEA_compare_results\\","`",colnames(scores_class)[i],"`","_compare.txt"),
              sep = '\t',quote = F,row.names = F)
}

save.image(str_c(path,"19_TCGA_k3_ssGSEA_immu_score_plot.RData"))
