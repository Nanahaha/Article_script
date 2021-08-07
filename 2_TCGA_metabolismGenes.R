## ------ 2751 metabolism-related genes based on entrez  ------ ##

library(sva) 
library(tidyverse)
library(ggsci) 
library(data.table) 

path <- getwd()

## read in entrez number of 38292 genes from g:Profiler website
entrez <- read.table(file = str_c(path,"gProfiler_hsapiens_2021-5-10下午2-22-15-38292ENTREZGENEacc.txt"),
                     header = T, fill = T, check.names = F, sep = ",")
entrez[1:4,]
colnames(entrez)[c(1,2)] <- c("gene_id","entrez_id")

## ------ extract 2751 metabolism-related genes from 38292 genes ------ ##
meta_genes <- fread(file = str_c(path,"metabolism-relevant_genes.txt"),
                    header = T, fill = T, sep = "\t", quote = "", data.table = T)
meta_genes[1:4,]
dim(meta_genes)
# 2752    4
colnames(meta_genes)[2:3] <- c("gene_symbol","entrez_id")
meta_genes2 <- meta_genes[-which(duplicated(meta_genes$gene_symbol)),]
dim(meta_genes2)
# 2751    4
meta_genes2[1:5,1:3]

mergeMetaGenes <- entrez$name[which(entrez$entrez_id %in% meta_genes2$entrez_id)]

tcga_exp_all <- data.frame(fread(file = str_c(path,"TCGA-OV_HTSeq_TPM_01survnoNA_rowsamp_log2.txt"),
                                 sep = "\t", header =T, stringsAsFactors = F, data.table = T, check.names = F), 
                           check.names = F)

tcga_exp_all[1:4,1:10]
table(tcga_exp_all[,3])
# Alive  Dead 
# 143   230 
colnames(tcga_exp_all)[3:4] <- c("Survival_Status","OS_Time")
which(is.na(tcga_exp_all$OS_Time))


## process expression matrix
rownames(tcga_exp_all) <- tcga_exp_all$Patients
tcga_exp <- t(tcga_exp_all[,-c(1:8)])
tcga_exp[1:4,1:4]
dim(tcga_exp)
# 38292   373 
tcga_exp2 <- data.frame("gene_symbol" = rownames(tcga_exp), tcga_exp, check.names = F)
tcga_exp2[1:4,1:4]

exp_meta_df <- tcga_exp2[which(tcga_exp2$gene_symbol %in% mergeMetaGenes),]
dim(exp_meta_df)
# 2718  374

write.table(exp_meta_df, file = str_c(path,"TCGA_metabolism_GeneMatrix.txt"),
            quote = F, sep = "\t", row.names = F)


save.image(str_c(path,"2_TCGA_metabolismGenes.RData"))




