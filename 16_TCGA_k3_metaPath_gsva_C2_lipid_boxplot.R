## ------ drawing boxplot for lipid metabolism pathway ------ ##

library(stringr)
library(dplyr)
library(RColorBrewer) 
library(ggpubr) 
library(reshape2) 
library(limma)  

path <- getwd()

data <- read.table(file = str_c(path,"k3\\TCGA_k3_metaPath_gsva_scores.txt"),
                   header = T, fill = T, sep = "\t", check.names = F)
data[1:4,1:3]
subclass <- read.table(file = str_c(path,"k3\\TCGA_consensus_k3_hclust_subclass2.txt"),
                       header = T, fill = T)
pathway_class <- read.table(file = str_c(path,"metabolism114_signatures_genes_KEGGclass.txt"),
                            header = T, fill = T, sep = "\t", check.names =F)

subclass$Class <- sapply(subclass$Class,function(x) str_c("C",x))
table(subclass$Class)
#  C1  C2  C3 
# 106 119 148  

CLASS <- "C2"
c_spe_pathways <- read.table(file = str_c(path,"k3\\metaPath_gsva_diff_analysis\\","TCGA_",CLASS,"_fc0fdr05_gsva.txt"),
                             sep = "\t", header = T, fill = T)

c_pathways_scores <- data[which(rownames(data) %in% rownames(c_spe_pathways)),]
c_pathways_scores [1:4,1:3]

path_class_data <- data.frame("Pathway_signatures" = rownames(c_pathways_scores),
                              c_pathways_scores, check.names = F)
path_class_data[1:4,1:3]
path_class_data2 <- merge(pathway_class[,c(1,3)], path_class_data,
                          by = "Pathway_signatures", all=F)
path_class_data2_order  <- path_class_data2[order(path_class_data2$Pathway_TCGAclass2),]
path_class_data2_order[1:4,1:3]
rownames(path_class_data2_order) <- path_class_data2_order$Pathway_signatures
path_class_data2_order[1:4,1:3]
table(path_class_data2_order$Pathway_TCGAclass2)
#Amino acid metabolism              Carbohydrate metabolism 
#4                                   11 
#Drug metabolism                    Energy metabolism 
#3                                    1 
#Glycan biosynthesis and metabolism                     Lipid metabolism 
#5                                   17 
#Metabolism of cofactors and vitamins                     Other metabolism 
#8                                   10                                  10


## --- pick out lipid, Output first, adjust the order manually, then input--- ##
onepath <- path_class_data2_order[which(path_class_data2_order$Pathway_TCGAclass2 == "Lipid metabolism"),]
dim(onepath)
onepath[1:4,1:5]
write.table(onepath, file = str_c(path,"k3\\TCGA_k3_metaPath_gsva_C2_lipid_scores.txt"),
            row.names = F, sep = "\t", quote = F)

plotDensities(onepath[,-c(1:2)],legend = F)

## Draw the unsorted first
temp1 <- t(onepath[,-c(1:2)])
temp2 <- data.frame("Patients" = row.names(temp1), temp1, check.names = F)
temp_class <- merge(subclass,temp2, by = "Patients", all = F)
temp_class_order <- temp_class[order(temp_class$Class),]
temp_class_order[1:4,1:3]
rownames(temp_class_order) <- temp_class_order$Patients
temp_class_order[1:5,1:3]

mydata <- melt(temp_class_order[,-1], id = "Class")

dev.size()
pdf(file = str_c(path,"k3\\TCGA_k3_metaPath_gsva_C2_lipid_boxplot.pdf"),
    width = 8, height = 4)
p <- ggplot(mydata, 
            aes(x = variable, y = value, fill = Class))+
  geom_boxplot()+
  theme_classic()+  
  scale_fill_manual(values = brewer.pal(5,"Set2"))+
  labs(y = "GSVA Score")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 6, angle = 45, hjust = 1))
print(p)
dev.off()

## analysis of difference
scores_class2 <- temp_class_order
colnames(scores_class2)[-c(1:2)] <- str_c("pathway",seq(1:(ncol(temp_class_order)-2)))

my_comparisons<-list(c("C1","C2"),c("C1","C3"),c("C2","C3"))

for (i in 3:ncol(scores_class2)) {
  formu <- as.formula(str_c(colnames(scores_class2)[i],"~ Class"))
  compa <- compare_means(formu, data = scores_class2)
  write.table(compa,file = str_c(path,"k3\\metaPath_gsva_C2_lipid_compare_results\\","`",colnames(temp_class_order)[i],"`","_compare.txt"),
              sep = '\t',quote = F,row.names = F)
}


## ------  Draw the sorted second ------ ##
onepath2 <- read.table(file = str_c(path,"k3\\TCGA_k3_metaPath_gsva_C2_lipid_scores_ranged.txt"),
                       header = T, fill = T, sep = "\t", check.names = F)
onepath2[1:4,1:5]
rownames(onepath2) <- onepath2$Pathway_signatures
temp12 <- t(onepath2[,-c(1:2)])
temp22 <- data.frame("Patients" = row.names(temp12), temp12, check.names = F)
temp_class2 <- merge(subclass,temp22, by = "Patients", all = F)
temp_class2_order <- temp_class2[order(temp_class2$Class),]
temp_class2_order[1:4,1:3]
rownames(temp_class2_order) <- temp_class2_order$Patients
temp_class2_order[1:5,1:4]


mydata2 <- melt(temp_class2_order[,-1], id = "Class")

dev.size()
pdf(file = str_c(path,"k3\\TCGA_k3_metaPath_gsva_C2_lipid_boxplot_ranged.pdf"),
    width = 8, height = 4)
p <- ggplot(mydata2, 
            aes(x = variable, y = value, fill = Class))+
  geom_boxplot()+
  theme_classic()+  
  scale_fill_manual(values = brewer.pal(5,"Set2"))+
  labs(y = "GSVA Score")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 6, angle = 45, hjust = 1))
print(p)
dev.off()


save.image(str_c(path,"16_TCGA_k3_metaPath_gsva_C2_lipid_boxplot.RData"))
