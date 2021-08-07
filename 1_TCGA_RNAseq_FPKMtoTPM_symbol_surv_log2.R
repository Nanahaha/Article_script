## ------  TCGAbiolinks download FPKM data，symbol transform  ------ ##

library(TCGAbiolinks) 
library(SummarizedExperiment)
library(stringr)

path <- getwd()
setwd(path)

PRO <- "TCGA-OV"

####################################################################################
query <- GDCquery(project = PRO, 
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq - FPKM",
                  legacy = F) 
GDCdownload(query,
            directory = str_c(path,"GDCdata\\"), 
            method="api") 

expdat <- GDCprepare(query) 

expdat_frame <- assay(expdat) 
expdat_frame[1:5,1:5]
expdat_df <- data.frame("gene_id"=rownames(expdat_frame), expdat_frame, check.names = F)
expdat_df[1:4,1:4]
write.table(expdat_df, file = str_c(PRO,"_HTSeq_FPKM.txt"),
            sep = "\t", quote = F, row.names = F) 


## --- symbol transform ------ ##
## read in transformed data from g:profiler website, merge
gProfiler <- read.table(file = str_c(path,"gProfiler_hsapiens_2021-5-9下午6-51-01.txt"),
                        sep = ",", header = T, fill = T, check.names = F)
gProfiler[1:4,1:4]
colnames(gProfiler)[c(1,3)] <- c("gene_id","gene_symbol")

fpkm_symbol <- merge(gProfiler[,c(1,3)], expdat_df, by = "gene_id", all = F)
fpkm_symbol2 <- fpkm_symbol[-which(fpkm_symbol$gene_symbol=="nan"),]

write.table(fpkm_symbol2,file = str_c(PRO,"_HTSeq_FPKM_idsymbol.txt"),
            quote = F, sep = "\t", row.names = F)

length(which(duplicated(fpkm_symbol2$gene_symbol)))
## delete duplicated gene symbol, retain max value
if(length(which(duplicated(fpkm_symbol2$gene_symbol))) > 0){
  fpkm_symbol3 <- aggregate(x = fpkm_symbol2, by = list(fpkm_symbol2$gene_symbol), FUN = max)
}else{
  fpkm_symbol3 <- data.frame("Group.1" = fpkm_symbol2$gene_symbol, fpkm_symbol2, check.names = F)
}
which(duplicated(fpkm_symbol3$gene_symbol))
rownames(fpkm_symbol3) <- fpkm_symbol3$gene_symbol

write.table(fpkm_symbol3[,-1],file = str_c(PRO,"_HTSeq_FPKM_idsymbol_noNAdup.txt"),
            quote = F, sep = "\t", row.names = F)


## transform FPKM to TPM：https://www.jianshu.com/p/9dfb65e405e8
fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
expdat_tpm <- apply(fpkm_symbol3[,-c(1:3)], 2, fpkmToTpm)
expdat_tpm[1:5,1:4]
expdat_tpm_df <- data.frame(fpkm_symbol3[,c(2:3)], expdat_tpm, check.names = F)
expdat_tpm_df[1:4,1:4]

write.table(expdat_tpm_df, file=str_c(PRO,"_HTSeq_TPM_idsymbol_noNAdup.txt"),
            sep = "\t", quote = F, row.names = F) 

## process sample name
expdat_tpm2 <- expdat_tpm
sample <- colnames(expdat_tpm2)
sample2 <- sapply(sample,function(x) str_c(unlist(str_split(x,pattern = '-'))[1:4],
                                           collapse = "-"))
colnames(expdat_tpm2) <- sample2
expdat_tpm2[1:5,1:3]
expdat_tpm2_df <- data.frame(expdat_tpm_df[,1:2], expdat_tpm2, check.names = F)
expdat_tpm2_df[1:4,1:4]

write.table(expdat_tpm2_df, file=str_c(PRO,"_HTSeq_TPM_idsymbol_noNAdup_samch.txt"),
            sep = "\t", quote = F, row.names = F) 


############################################
## download survival data
clinical <- GDCquery_clinic(project = PRO,
                            type = "clinical")

write.table(clinical,file=str_c(PRO,"_clinical.txt"),sep = "\t", quote=F)

## extract survival data
surv <- clinical[,c("submitter_id","days_to_last_follow_up",
                    "vital_status","days_to_death","figo_stage")]
surv[1:4,1:4]

surv$days_to_death[is.na(surv$days_to_death)] <- surv$days_to_last_follow_up[is.na(surv$days_to_death)]
surv$Days <- surv$days_to_death
surv$Status <- surv$vital_status
colnames(surv)[1] <- "Patients_cli"

write.table(surv, file=str_c(PRO,"_surv.txt"),sep = "\t", quote=F)



## ------ retain 01 sample, 据merge expression and survival data ------ ##
expdat_tpm2[1:4,1:4]
COLNAME <- sapply(colnames(expdat_tpm2), function(x) str_sub(x,end = -2))
colnames(expdat_tpm2) <- COLNAME
prim_col <- which(sapply(COLNAME, function(x) str_split(x,pattern = "-")[[1]][4]) == "01")
expdat_tpm2_prim <- expdat_tpm2[,prim_col]
expdat_tpm2_prim[1:5,1:3]
expdat_tpm2_prim_rowsam <- t(expdat_tpm2_prim)
expdat_tpm2_prim_rowsam[1:5,1:4]
expdat_tpm2_prim_rowsam2 <- data.frame("Patients_cli" = sapply(rownames(expdat_tpm2_prim_rowsam),
                                                               function(x) str_sub(x,end = -4)),
                                       "Patients" = rownames(expdat_tpm2_prim_rowsam), 
                                       expdat_tpm2_prim_rowsam, 
                                       check.names = F)
expdat_tpm2_prim_rowsam2[1:5,1:3]

expdat_tpm2_prim_rowsam2_surv <- merge(surv, expdat_tpm2_prim_rowsam2, 
                                       by = "Patients_cli", all = F)
expdat_tpm2_prim_rowsam2_surv[1:4,1:10]
dim(expdat_tpm2_prim_rowsam2_surv)

write.table(expdat_tpm2_prim_rowsam2_surv, file = str_c(PRO,"_HTSeq_TPM_01allsurv_rowsamp.txt"),
            quote = F, sep = "\t", row.names =F)


## delete sample without survival data
expdat_tpm2_prim_rowsam2_surv[1:3,1:3]
if(sum(is.na(expdat_tpm2_prim_rowsam2_surv$days_to_death)) > 0){
  expdat_tpm2_prim_rowsam2_surv_noNA <- expdat_tpm2_prim_rowsam2_surv[-which(is.na(expdat_tpm2_prim_rowsam2_surv$days_to_death)),]
}else{
  expdat_tpm2_prim_rowsam2_surv_noNA <- expdat_tpm2_prim_rowsam2_surv
}
expdat_tpm2_prim_rowsam2_surv_noNA [1:4,1:4]
dim(expdat_tpm2_prim_rowsam2_surv_noNA)

write.table(expdat_tpm2_prim_rowsam2_surv_noNA, 
            file = str_c(path,PRO,"_HTSeq_TPM_01survnoNA_rowsamp.txt"),
            quote = F, sep = "\t", row.names =F)

## separating out survival data
write.table(expdat_tpm2_prim_rowsam2_surv_noNA[,1:8], 
            file = str_c(path,PRO,"_surv373.txt"),
            quote = F, sep = "\t", row.names =F)


## log(TPM+1) transform
expdat_tpm_surv_log <- expdat_tpm2_prim_rowsam2_surv_noNA
expdat_tpm_surv_log[1:4,1:10]
expdat_tpm_surv_log[,-c(1:8)] <- log2(expdat_tpm_surv_log[,-c(1:8)]+1)
expdat_tpm_surv_log[1:4,c(9:12)]
which(expdat_tpm_surv_log<0)
# integer(0)
write.table(expdat_tpm_surv_log, file = str_c(path,PRO,"_HTSeq_TPM_01survnoNA_rowsamp_log2.txt"),
            quote = F, sep = "\t", row.names =F)

save.image(str_c(path,"1_TCGA_RNAseq_FPKMtoTPM_symbol_surv_log2.RData"))
