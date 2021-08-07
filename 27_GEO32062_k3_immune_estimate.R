## ESTIMATE algorithm for GEO dataset

library(estimate)
library(stringr)
library(RColorBrewer)
library(ggplot2)
library(ggpubr) 
library(plyr)
library(data.table)


SET <- "GSE32062"
PATH <- getwd()
path <- str_c(PATH,SETS,"\\")
setwd(path)

data <- data.frame(fread(file = str_c(path,SET,"_rma_quatile_prob_symbol_exprSet.txt"),
                         header = T, fill = T, check.names = F), check.names = F) 
data[1:4,1:9]
rownames(data) <- data$GENE_SYMBOL
write.table(data[,-c(1:2)], file = str_c(path,SET,"_rma_quatile_prob_symbol_exprSet_forestimate.txt"),
            row.names = T, sep = "\t", quote = F)

filterCommonGenes(input.f = str_c(path,SET,"_rma_quatile_prob_symbol_exprSet_forestimate.txt"),
                  output.f = str_c(path,SET,"_rma_quatile_prob_symbol_exprSet_forestimate.gct"),
                  id = "GeneSymbol")

estimateScore(input.ds = str_c(path,SET,"_rma_quatile_prob_symbol_exprSet_forestimate.gct"),
              output.ds = str_c(path,"k3\\",SET,"_k3_estimate_scores.gct"),
              platform = "agilent")

scores <- read.table(file = str_c(path,"k3\\",SET,"_k3_estimate_scores.gct"),
                     skip = 2, header = T, check.names = F) 
scores[1:3,1:4]
rownames(scores) <- scores[,1]
scores_t <- t(scores[,-c(1:2)])
scores_t[1:4,]
scores_t2 <- data.frame("Patients" = rownames(scores_t), scores_t, check.names = F)
scores_t2[1:4,]
write.table(scores_t2, file = str_c(path,"k3\\",SET,"_k3_estimate_scores_rowsample.txt"),
            sep = "\t", quote = F, row.names = F)


# --- statistical analysis, boxplot --- #
subclass <- read.table(file = str_c(path,"k3\\",SET,"_k3_NTP_45classifier_subclass.txt"),
                       header = T, fill = T)


subclass[1:4,]

scores_class <- merge(subclass, scores_t2, by = "Patients", all = F)

scores_class$Class <- factor(scores_class$Class,levels = c("C1","C2","C3"),
                             ordered = T)
scores_class[1:4,]


for (i in 3:4) {
  formu <- as.formula(str_c(colnames(scores_class)[i],"~ Class"))
  compa <- compare_means(formu, data = scores_class)
  write.table(compa,file = str_c(path,"k3\\",SET,"_k3_",colnames(scores_class)[i],"_compare.txt"),
              sep = '\t',quote = F,row.names = F)
  
  my_comparisons_stro<-list(c("C1","C2"),c("C1","C3"),c("C2","C3"))
  my_comparisons_immo<-list(c("C1","C2"),c("C1","C3"),c("C2","C3"))
  
  pdf(file = str_c(path,"k3\\",SET,"_k3_estimate_",colnames(scores_class)[i],"_boxplot.pdf"),
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

save.image(str_c(path,"27_",SET,"_k3_immune_estimate.RData"))

