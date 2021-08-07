## ESTIMATE algorithm estimate immune score and stromal score: https://www.sohu.com/a/396935566_777125

library(estimate)
library(stringr)
library(RColorBrewer)
library(ggplot2)
library(ggpubr) 
library(plyr)
library(data.table)

path <- getwd()
setwd(path)

data <- data.frame(fread(file = str_c(path,"TCGA-OV_HTSeq_TPM_01survnoNA_rowsamp_log2.txt"),
                         header = T, fill = T, check.names = F), check.names = F) 
data[1:4,1:9]
rownames(data) <- data$Patients
data_colsamp <- t(data[,-c(1:8)])
data_colsamp[1:4,1:4]
write.table(data_colsamp, file = str_c(path,"TCGA-OV_HTSeq_TPM_01survnoNA_log2_colsamp_forestimate.txt"),
            row.names = T, sep = "\t", quote = F)

data_colsamp_df <- data.frame("gene_symbol" = rownames(data_colsamp), data_colsamp, check.names = F)
data_colsamp_df[1:4,1:4]
write.table(data_colsamp_df, file = str_c(path,"TCGA-OV_HTSeq_TPM_01survnoNA_log2_colsamp.txt"),
            row.names = F, sep = "\t", quote = F)

filterCommonGenes(input.f = str_c(path,"TCGA-OV_HTSeq_TPM_01survnoNA_log2_colsamp_forestimate.txt"),
                  output.f = str_c(path,"TCGA-OV_HTSeq_TPM_01survnoNA_log2_colsamp_forestimate.gct"),
                  id = "GeneSymbol")

estimateScore(input.ds = str_c(path,"TCGA-OV_HTSeq_TPM_01survnoNA_log2_colsamp_forestimate.gct"),
              output.ds = str_c(path,"k3\\TCGA_k3_estimate_scores.gct"),
              platform = "illumina")

scores <- read.table(file = str_c(path,"k3\\TCGA_k3_estimate_scores.gct"),
                     skip = 2, header = T, check.names = F) 
scores[1:4,1:4]
rownames(scores) <- scores[,1]
scores_t <- t(scores[,-c(1:2)])
scores_t[1:4,]
newname <- sapply(rownames(scores_t),function(x){
  a <- strsplit(x, split = ".", fixed = T)
  a2 <- str_c(a[[1]][1],"-",a[[1]][2],"-",a[[1]][3],"-",a[[1]][4])
})
rownames(scores_t) <- newname
scores_t[1:4,]
scores_t_df <- data.frame("Patients" = rownames(scores_t), scores_t, check.names = F)
scores_t_df[1:4,1:4]
write.table(scores_t_df, file = str_c(path,"k3\\TCGA_k3_estimate_scores_rowsamp.txt"),
            sep = "\t", quote = F, row.names = F)


# --- boxplot --- #
subclass <- read.table(file = str_c(path,"k3\\TCGA_consensus_k3_hclust_subclass2.txt"),
                       header = T, fill = T)
subclass$Class <- sapply(subclass$Class,function(x) str_c("C",x))

scores_class <- merge(subclass, scores_t_df, by = "Patients", all = F)

scores_class$Class <- factor(scores_class$Class,levels = c("C1","C2","C3"),
                             ordered = T)
scores_class[1:4,]

write.table(scores_class, file = str_c(path,"k3\\TCGA_k3_estimate_scores_rowsamp_class.txt"),
            sep = "\t", quote = F, row.names = F)

for (i in 3:4) {
  formu <- as.formula(str_c(colnames(scores_class)[i],"~ Class"))
  compa <- compare_means(formu, data = scores_class)
  write.table(compa,file = str_c(path,"k3\\TCGA_K3_",colnames(scores_class)[i],"_compare.txt"),
              sep = '\t',quote = F,row.names = F)
  my_comparisons_stro<-list(c("C1","C2"),c("C1","C3"),c("C2","C3"))
  my_comparisons_immo<-list(c("C1","C2"),c("C1","C3"),c("C2","C3"))
  
  pdf(file = str_c(path,"k3\\TCGA_k3_estimate_",colnames(scores_class)[i],"_boxplot.pdf"),
      width = 4, height = 4)
  
  if(i == 3){
    p <- ggplot(scores_class,aes(x = Class,y = StromalScore,fill = Class))+
      geom_boxplot()+
      theme_classic()+ 
      scale_fill_manual(values = brewer.pal(5,"Set2"))+
      labs(y = "Stromal Score")+
      theme(axis.title.x = element_blank(),
            axis.title.y = element_text(size = 15),
            axis.text.x = element_text(size = 15),
            axis.text.y = element_text(size = 15))+
      stat_compare_means(comparisons = my_comparisons_stro,
                         tip.length = 0.01,
                         label = "p.signif")
  }else{
    p <- ggplot(scores_class,aes(x = Class,y = ImmuneScore,fill = Class))+
      geom_boxplot()+
      theme_classic()+ 
      scale_fill_manual(values = brewer.pal(5,"Set2"))+
      labs(y = "Immune Score")+
      theme(axis.title.x = element_blank(),
            axis.title.y = element_text(size = 15),
            axis.text.x = element_text(size = 15),
            axis.text.y = element_text(size = 15))+
      stat_compare_means(comparisons = my_comparisons_immo,
                         tip.length = 0.01,
                         label = "p.signif")
  }
  
  print(p)
  
  dev.off()
}

save.image(str_c(path,"17_TCGA_k3_immune_estimate.RData"))

