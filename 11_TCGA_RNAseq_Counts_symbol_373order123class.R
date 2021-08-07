## ------- process matrix for GSEA analysis, Sort the classes by 123 and output ------ ##

library(stringr)
library(data.table)

path <- getwd()

count <- data.frame(fread(file = str_c(path,"TCGA-OV_HTSeq_Counts_symbol_orderclass_373patients_colsamp.txt"),
                          header = T, fill = T, check.names = F), check.names = F)

count[1:4,1:4]
rownames(count) <- count$gene_symbol

pheno <- read.table(file = str_c(path,"k3\\TCGA_consensus_k3_hclust_subclass2.txt"),
                    header = T, fill = T, check.names = F)

pheno[1:4,]

count_rowsamp <- t(count[,-1])
count_rowsamp[1:4,1:4]

count_rowsamp_df <- data.frame("Patients" = rownames(count_rowsamp),
                               count_rowsamp)

pheno_count <- merge(pheno, count_rowsamp_df, by = "Patients", all = F)
dim(pheno_count)

pheno_count2 <- pheno_count[order(pheno_count$Class),]
pheno_count2[1:10,1:3]

write.table(pheno_count2[,1:2], file = str_c(path,"k3\\TCGA_consensus_k3_hclust_subclass2_ordered.txt"),
            sep = "\t", quote = F, row.names = F)

rownames(pheno_count2) <- pheno_count2$Patients
pheno_count2_colsamp <- t(pheno_count2[,-c(1:2)])
pheno_count2_colsamp[1:4,1:4]

pheno_count2_colsamp_df <- data.frame("gene_symbol" = rownames(pheno_count2_colsamp),
                                      pheno_count2_colsamp, check.names = F)
pheno_count2_colsamp_df[1:4,1:4]
write.table(pheno_count2_colsamp_df, file = str_c(path,"TCGA-OV_HTSeq_Counts_symbol_order123class_373patients_colsamp.txt"),
            sep = "\t", quote = F, row.names = F)

save.image(str_c(path,"11_TCGA_RNAseq_Counts_symbol_373order123class.RData"))
