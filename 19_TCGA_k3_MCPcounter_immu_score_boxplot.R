## ------ Merging and statistical analysis of the data after Timer ------ ##

library(stringr)
library(dplyr)
library(ggplot2)
library(RColorBrewer) 
library(ggpubr) 
library(reshape2) 
library(limma)  

path <- getwd()

filelist <- list.files(path = str_c(path,"k3\\TIMER2_results\\"),
                       pattern = ".csv",
                       full.names = T)
mydata <- lapply(filelist, function(x){
  read.csv(x, header=T, fill = T, check.names = F)})  

timer_data <- mydata[[1]]
for (i in 2:length(mydata)) {
  timer_data <- merge(timer_data, mydata[[i]], by = "cell_type", all = F)
}

timer_data2 <- timer_data %>% mutate(Algorithm = NA,
                                     .after = "cell_type")
colnames(timer_data)

timer_data2$Algorithm <- sapply(timer_data2$cell_type, function(x) str_split(x, pattern = "_")[[1]][2])

timer_data2$cell_type <- sapply(timer_data2$cell_type, function(x) str_split(x, pattern = "_")[[1]][1])

timer_data2 <- timer_data2[order(timer_data2$Algorithm),]
timer_data2[1:5,1:5]
write.table(timer_data2, file = str_c(path,"k3\\TCGA_k3_Timer_results.txt"),
            sep = "\t", quote = F, row.names = F)

table(timer_data2$Algorithm)
# pick up MCPCounter
mcprow <- which(timer_data2$Algorithm=="MCPCOUNTER")

mcp_re <- timer_data2[mcprow,]
mcp_re[1:4,1:4]
rownames(mcp_re) <- mcp_re$cell_type

plotDensities(mcp_re[,-c(1:2)],legend = F)

mcp_re_log <- as.matrix(log2(mcp_re[,-c(1:2)]))
mcp_re_log[1:6,1:4]
mcp_re_log[which(mcp_re_log=="-Inf")] <- 0
plotDensities(mcp_re_log,legend = F)

mcp_re_samplerow <- t(mcp_re_log)
mcp_re_samplerow[1:5,1:5]

mcp_re_samplerow3 <- data.frame("Patients" = rownames(mcp_re_samplerow),
                                mcp_re_samplerow,
                                check.names = F)
## 导入分组信息
subclass <- read.table(file = str_c(path,"k3\\TCGA_consensus_k3_hclust_subclass2.txt"),
                       header = T, fill = T)
subclass$Class <- sapply(subclass$Class,function(x) str_c("C",x))

scores_class <- merge(subclass, mcp_re_samplerow3, by = "Patients", all = F)
scores_class[1:4,1:4]
min(scores_class[,-c(1:2)])
max(scores_class[,-c(1:2)])
which(scores_class[,-c(1:2)]==0)

write.table(scores_class, file = str_c(path,"k3\\TCGA_k3_MCPcounter_immu_scoreslog_class.txt"),
            sep = "\t", quote = F, row.names = F)

mydata <- melt(scores_class[,-1], id = "Class")
mydata[1:4,]

min(mydata$value)
max(mydata$value)

#Cancer associated fibroblast，Endothelial cell和Myeloid dendritic cell，,"T cell")
cells <- c("T cell","Cancer associated fibroblast","Endothelial cell","Myeloid dendritic cell")

dev.size()
pdf(file = str_c(path,"k3\\TCGA_k3_MCP_immu_score_boxplot.pdf"),
    width = 3.5, height = 3.5)
p <- ggplot(mydata[which(mydata$variable %in% cells),], 
            aes(x = variable, y = value, fill = Class))+
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
colnames(scores_class2)[-c(1:2)] <- c("B","fibroblast","cytotoxicity","endothelial","mycophage","monocyte","myeloid",
                                      "neutrophil","NK","CD8","T")
my_comparisons<-list(c("C1","C2"),c("C1","C3"),c("C2","C3"))

for (i in 3:ncol(scores_class2)) {
  formu <- as.formula(str_c(colnames(scores_class2)[i],"~ Class"))
  compa <- compare_means(formu, data = scores_class2)
  write.table(compa,file = str_c(path,"k3\\TIMER2_compare_results\\",colnames(scores_class2)[i],"_compare.txt"),
              sep = '\t',quote = F,row.names = F)
}

save.image(str_c(path,"19_TCGA_k3_MCPcounter_immu_score_plot.RData"))
