## ------- immune checkpoint genes in GEO dataset----- ##

library(tidyverse)
library(data.table)
library(limma)  
library(car)  
library(RColorBrewer) 
library(reshape2) 
library(ggpubr) 

SET <- "GSE32062"
PATH <- getwd()
path <- str_c(PATH,SETS,"\\")

data <- data.frame(fread(file = str_c(path,SET,"_rma_quatile_prob_symbol_exprSet.txt"),
                         sep = "\t", header = T, fill = T, check.names = F), check.names = F)
data[1:4,1:4]

rownames(data)  <- data$GENE_SYMBOL
data[1:4,1:4]

plotDensities(data[,-c(1:2)], legend = F)

data_rowsample <- t(data[,-c(1:2)])
data_rowsample[1:4,1:4]

checkpoint <- c("PDCD1","CD274","PDCD1LG2","CCL2","CTLA4","CXCR4","IL1A",
                "IL6","LAG3","TGFB1","TNFRSF4","TNFRSF9","TNFSF4","CD4","CD86",
                "CD80","CD276","VTCN1")
if(length(which(!checkpoint %in% colnames(data_rowsample))) > 0){
  data_rs_cp <- data_rowsample[,checkpoint[-which(!checkpoint %in% colnames(data_rowsample))]]
}else{
  data_rs_cp <- data_rowsample[,checkpoint]
}

data_rs_cp[1:4,1:4]

data_rs_cp_df <- data.frame("Patients" = rownames(data_rs_cp), 
                            data_rs_cp, check.names = F)
data_rs_cp_df[1:4,1:4]

subclass <- read.table(file = str_c(path,"k3\\",SET,"_k3_NTP_45classifier_subclass.txt"),
                       header = T, fill = T)
subclass[1:4,]

data_rs_cp_df_class <- merge(subclass, data_rs_cp_df, by = "Patients", all = F)
data_rs_cp_df_class[1:5,1:5]
# 将class变成因子
data_rs_cp_df_class$Class <- factor(data_rs_cp_df_class$Class,
                                    levels = c("C1","C2","C3"),
                                    ordered = T)

## ------ boxplot图 ------ ##
mydata <- reshape2::melt(data_rs_cp_df_class[,-1], id = "Class")
mydata[1:4,]

min(mydata$value)
max(mydata$value)

dev.size()
pdf(file = str_c(path,"k3\\",SET,"_k3_immu_checkpoint_genes_boxplot.pdf"),
    width = 12, height = 5)
p <- ggplot(mydata, aes(x = variable, y = value, fill = Class))+
  geom_boxplot()+
  theme_classic()+
  scale_fill_manual(values = brewer.pal(3,"Set2"))+
  labs(y = "Relative expression (log2)")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
print(p)
dev.off()

data_rs_cp_df_class[1:4,]

colnames(data_rs_cp_df_class)

my_comparisons<-list(c("C1","C2"),c("C1","C3"),c("C2","C3"))

for (i in 3:ncol(data_rs_cp_df_class)) {
  formu <- as.formula(str_c(colnames(data_rs_cp_df_class)[i],"~ Class"))
  compa <- compare_means(formu, data = data_rs_cp_df_class)
  write.table(compa,file = str_c(path,"k3\\immuCheckpointGenes_compare_results\\",colnames(data_rs_cp_df_class)[i],"_compare.txt"),
              sep = '\t',quote = F,row.names = F)
}


save.image(str_c(path,"28_",SET,"_k3_immu_checkpoint_genes_plot.RData"))
