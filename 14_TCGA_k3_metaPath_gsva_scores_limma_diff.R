## ------ Analysis of gsva score difference, p.value = 0.05, Not considering logfc ------ ##致

library(stringr)
library(limma)

path <- getwd()

data <- read.table(file = str_c(path,"k3\\TCGA_k3_metaPath_gsva_scores.txt"),
                   header = T, fill = T, sep = "\t", check.names = F)
data[1:4,1:4]
plotDensities(data, legend = F)

subclass <- read.table(file = str_c(path,"k3\\TCGA_consensus_k3_hclust_subclass2.txt"),
                       header = T, fill = T)
str(data)
data[1:4,1:4]
subclass[1:4,]
C1 <- which(subclass[,2]=="1")
C2 <- which(subclass[,2]=="2")
C3 <- which(subclass[,2]=="3")

mean(as.matrix(data[4,C1]))/mean(as.matrix(data[4,-C1]))

## Comparing single groups with other samples
class_nrow <- list(C1,C2,C3)
k3_limmmafit <- list()
for (i in 1:3) {
  a1 <- subclass
  
  now_class <- class_nrow[[i]]
  row_n <- 1:nrow(a1)
  rest_class <- row_n[-which(row_n %in% class_nrow[[i]])] 
  
  a1[now_class,2] <- "c1"
  a1[rest_class,2] <- "rest"
  table(a1$Class)
  
  a1$Class <- factor(a1$Class,
                     levels = c("c1","rest"),
                     ordered = T)
  design <- model.matrix(~0+a1$Class)
  colnames(design) <- levels(a1$Class)
  rownames(design) <- a1$Sample
  
  contrast_matrix <- makeContrasts("c1-rest",levels = design)
  
  fit <- lmFit(data,design)
  fit2 <- contrasts.fit(fit,contrast_matrix)
  fit2 <- eBayes(fit2)
  k3_limmmafit[[i]] <- fit2
  
  allGenes <- topTable(fit2, coef = "c1-rest", number = Inf, 
                       adjust.method = "fdr", sort.by = "logFC")
  write.table(allGenes ,file = str_c(path,"k3\\metaPath_gsva_diff_analysis\\","TCGA_C",i,"_diffall_gsva.txt"),
              sep = '\t',quote = F,row.names = T)
  
  diffGenes <- topTable(fit2, coef = "c1-rest", number = Inf, 
                        adjust.method = "fdr", sort.by = "logFC", 
                        p.value = 0.05, lfc = 0)
  write.table(diffGenes,file = str_c(path,"k3\\metaPath_gsva_diff_analysis\\","TCGA_C",i,"_fc0fdr05_gsva.txt"),
              sep = "\t", quote = F, row.names = T)
}

save.image(str_c(path,"14_TCGA_k3_metaPath_gsva_scores_limma_diff.RData"))


