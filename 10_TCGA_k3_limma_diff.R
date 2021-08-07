## ------ Analysis of differences between subgroups ------ ##

library(stringr)
library(limma)
library(edgeR)
library(data.table)

path <- getwd()

data2 <- data.frame(fread(file = str_c(path,"TCGA-OV_HTSeq_Counts_idsymbol_noNAdup_samch.txt"),
                          header = T, fill = T, check.names = F), 
                    check.names = F) 
data2[1:4,1:4]
subclass <- read.table(file = str_c(path,"k3\\TCGA_consensus_k3_hclust_subclass2.txt"),
                       header = T, fill = T)

## row gene, column sample
data2[1:4,1:4]
rownames(data2) <- data2$gene_symbol
data <- data2[,-c(1:2)]
data[1:4,1:4]

## delete last letter of sample name, check weather exist duplicated sample
newcolname <- sapply(colnames(data), function(x){
  str_sub(x, end = -2)
})
which(duplicated(newcolname))
# integer(0)
colnames(data) <- newcolname
data[1:4,1:4]
data3 <- data[,which(colnames(data) %in% subclass$Patients)]
data3[1:4,1:4]

## Transpose to sample as row
data5 <- t(data3)
data5[1:4,1:4]
data5_df <- data.frame("Patients"=rownames(data5), data5, check.names = F)

## Reorder the dataset according to the sample order of the subclass
data4 <- data3[,subclass$Patients]
data4[1:4,1:4]
data4_df <- data.frame("gene_symbol"=rownames(data4), data4, check.names = F)
data4_df[1:4,1:4]
write.table(data4_df, file = str_c(path,"TCGA-OV_HTSeq_Counts_symbol_orderclass_373patients_colsamp.txt"),
            sep = "\t", quote = F, row.names = F)

subclass[1:4,]
C1 <- which(subclass[,2]=="1")
C2 <- which(subclass[,2]=="2")
C3 <- which(subclass[,2]=="3")

#load(str_c(path,"TCGA_ICGC_k3_limma_diff.RData"))

mean(as.matrix(data4[which(rownames(data4)=="COL1A1"),C2]))
mean(as.matrix(data4[which(rownames(data4)=="COL1A1"),-C2]))

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
  rownames(design) <- a1$Patients
  
  DGElist <- DGEList(counts = data4, group = a1$Class)
  keep_gene <- rowSums(cpm(DGElist) > 10) >= 2
  table(keep_gene)
  DGElist <- DGElist[keep_gene, ,keep.lib.sizes =FALSE]
  
  DGElist <- calcNormFactors( DGElist )
  v <- voom(DGElist, design, plot = TRUE, normalize = "quantile")
  fit <- lmFit(v, design, method = "ls")
  contrast <- makeContrasts("c1-rest", levels = design)
  contrast
  
  fit2 <- contrasts.fit(fit, contrast)
  fit2 <- eBayes(fit2)
  k3_limmmafit[[i]] <- fit2
  
  allGenes <- topTable(fit2, coef = "c1-rest", number = Inf, 
                       adjust.method = "fdr", sort.by = "logFC")
  write.table(allGenes ,file = str_c(path,"k3\\diff_analysis\\","TCGA_C",i,"_diffall_genes.txt"),
              sep = "\t", quote = F, row.names = T)
  
  diffGenes <- topTable(fit2, coef = "c1-rest", number = Inf, 
                        adjust.method = "fdr", sort.by = "logFC", 
                        p.value = 0.01, lfc = 1.5)
  write.table(diffGenes,file = str_c(path,"k3\\diff_analysis\\","TCGA_C",i,"_fc15fdr01_genes.txt"),
              sep = "\t", quote = F, row.names = T)
  
  diffGenes_logFC0 <- diffGenes[which(diffGenes$logFC > 0),]
  write.table(diffGenes_logFC0,file = str_c(path,"k3\\diff_analysis\\","TCGA_C",i,"_fc15fdr01_increasegenes.txt"),
              sep = "\t", quote = F, row.names = T)
  
}

save.image(str_c(path,"10_TCGA_k3_limma_diff.RData"))
