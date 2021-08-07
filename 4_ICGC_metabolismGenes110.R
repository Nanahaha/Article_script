
library(dplyr)
library(stringr)
library(data.table) 

##  read in ICGC matrix 
path <- getwd()
icgc_exp_all <- fread(file = str_c(path,"ICGC_OV_AU_Psolid_survnoNA_TPM_log2.txt"),
                      sep = "\t", header =T, stringsAsFactors = F, data.table = T, check.names = F)

icgc_exp_all[1:4,1:10]
table(icgc_exp_all[,3])
# alive deceased 
# 15       66  
colnames(icgc_exp_all)[2:3] <- c("OS_Time","Survival_Status")
icgc_exp_all$Survival_Status[which(icgc_exp_all$Survival_Status=="alive")] <- "Alive"
icgc_exp_all$Survival_Status[which(icgc_exp_all$Survival_Status=="deceased")] <- "Dead"
which(is.na(icgc_exp_all$OS_Time))
table(icgc_exp_all[,3])
#Alive  Dead 
#15    66 

## extract survival data
surv_df <- data.frame("Patients"=icgc_exp_all$exp_sample, icgc_exp_all[,2:3], check.names = F)
write.table(surv_df, file = str_c(path,"ICGC_OV_AU_surv81.txt"),
            sep = "\t", quote = F, row.names = F)

icgc_exp <- t(icgc_exp_all[,-c(1:7)])
icgc_exp[1:4,1:4]
colnames(icgc_exp) <- icgc_exp_all$exp_sample
icgc_exp[1:4,1:4]
dim(icgc_exp)
# 56206   81
icgc_exp2 <- data.frame("gene_id" = rownames(icgc_exp), icgc_exp, check.names = F)
icgc_exp2[1:4,1:4]


################################################################################################
## read in transformed data from g:profiler website, merge
path2 <- "F:\\Laboratory\\OvarianCancer_MelacularClassification\\OvarianCancer_TCGA\\"
gProfiler <- read.table(file = str_c(path2,"gProfiler_hsapiens_2021-5-9下午6-51-01.txt"),
                        sep = ",", header = T, fill = T, check.names = F)
gProfiler[1:4,1:4]
colnames(gProfiler)[c(1,3)] <- c("gene_id","gene_symbol")

icgc_exp2_symbol <- merge(gProfiler[,c(1,3)], icgc_exp2, by = "gene_id", all = F)
dim(icgc_exp2_symbol)
# 49986  83
icgc_exp2_symbol[1:5,1:3]
icgc_exp2_symbol2 <- icgc_exp2_symbol[-which(icgc_exp2_symbol$gene_symbol=="nan"),]
dim(icgc_exp2_symbol2)
# 36405   83
icgc_exp2_symbol2[1:4,1:3]
icgc_exp2_symbol2$gene_symbol[which(duplicated(icgc_exp2_symbol2$gene_symbol))]
# delete duplicated gene symbol, retain matix value
if(sum(duplicated(icgc_exp2_symbol2$gene_symbol)) > 0){
  icgc_exp2_symbol3 <- aggregate(x = icgc_exp2_symbol2, by = list(icgc_exp2_symbol2$gene_symbol), FUN = max)
}else{
  icgc_exp2_symbol3 <- data.frame("Group.1" = icgc_exp2_symbol2$gene_symbol, icgc_exp2_symbol2, check.gene_symbols = F)
}

which(duplicated(icgc_exp2_symbol3$gene_symbol))
rownames(icgc_exp2_symbol3) <- icgc_exp2_symbol3$gene_symbol
icgc_exp2_symbol3[1:5,1:4]

write.table(icgc_exp2_symbol3[,-1],file = str_c(path,"ICGC_OV_AU_TPM_log2_81samp_idsymbol.txt"),
            quote = F, sep = "\t", row.names = F)

###############################################################################################
## extract 110 metaboolism-related genes 
meta_final_genesname <- read.table(file = str_c(path2,"TCGA_metabolism_MADcox_110genes.txt"),
                                         sep = "\t", header = F, check.names = F)[1:2,]
meta_final_genesname[,1:4]
meta_final_genesname2 <- t(meta_final_genesname[,-1])
meta_final_genesname2[1:4,]
colnames(meta_final_genesname2)[1] <- "gene_symbol"
GEO_final_gene <- merge(meta_final_genesname2, icgc_exp2_symbol3[,-c(1:2)], by = "gene_symbol", all.x = T)[,-2]
dim(GEO_final_gene)
# 110  82
GEO_final_gene[1:5,1:5]
GEO_final_gene$gene_symbol[is.na(GEO_final_gene[,2])]
# no "SRXN1"
GEO_final_gene2 <- na.omit(GEO_final_gene)
dim(GEO_final_gene2)
# 109 82  109 genes
which(is.na(GEO_final_gene2))
GEO_final_gene2[1:4,1:4]
write.table(GEO_final_gene2, file = str_c(path,"ICGC_metabolism_colsamp_109genes.txt"),
            quote = F, sep = "\t", row.names = F)

save.image(str_c(path,"4_ICGC_metabolismGenes110.RData"))
