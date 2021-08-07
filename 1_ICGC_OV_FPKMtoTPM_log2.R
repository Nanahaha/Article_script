## ------ transform ICGC FPKM to TPM ------ ##

library(tidyverse)

PATH <- getwd()
setwd(PATH)

## extract sample information 
specimen_file <- read.delim(file = "specimen.OV-AU.tsv",
                            sep = "\t",header = T,fill = T) 
sepcimen_simplify_file <- specimen_file %>% mutate(specimen_type=ifelse(specimen_type=="Primary tumour - solid tissue","C","N")) %>%
  select(1,5,7)

## read in expression matrix (Manually extracted exp_seq.OV-AU.tsv in linux)
exp_file <- read_delim(file = "exp_seq_OV-AU_simplify.txt",
                       delim = "\t",col_names =  T) 
exp_file[1:5,]

## merge
merge_df <- merge(exp_file,sepcimen_simplify_file,by="icgc_specimen_id",all.x=TRUE)
merge_df[1:5,]

## gather format
tmp <- merge_df[,c(1,5,6,3,4)] %>% tibble::as_tibble() %>% unite(sample_id,icgc_specimen_id,icgc_donor_id.y,specimen_type,sep="-") 
tmp[1:5,]
class(tmp)
length(unique(exp_file$icgc_donor_id))    # [1] 93
length(unique(exp_file$icgc_specimen_id)) # [1] 111

exp_matrix <- tmp %>% group_by(sample_id) %>% mutate(id=1:n()) %>% 
  spread(sample_id,normalized_read_count) %>% select(-"id")

exp_matrix[1:5,1:3]
#rownames(exp_matrix)<-exp_matrix[,1]
write.table(exp_matrix,file="ICGC_OV_AU_exp_seq_FPKM.txt",
            sep = '\t',quote = F,row.names = F) 

## transform FPKM to TPM：https://www.jianshu.com/p/9dfb65e405e8
fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
exp_tpm <- apply(exp_matrix[,-1],2,fpkmToTpm)
exp_tpm[1:3,1:3]
exp_tpm <- data.frame("gene_id" = exp_matrix$gene_id, exp_tpm, check.names = F)

write.table(exp_tpm,file = "ICGC_OV_AU_exp_seq_TPM.txt",
            sep = '\t',quote = F, row.names = F) 

## delet normal sample
exp_tpm <- read.delim(file = "ICGC_OV_AU_exp_seq_TPM.txt",
                      header = T,fill = F,check.names = F)

cancer_col <- which(sapply(colnames(exp_tpm),function(x) str_split(x,pattern = '-')[[1]][3]) == "C")
exp_tpm_cancer <- exp_tpm[,c(1,cancer_col)]

write.table(exp_tpm_cancer,file = "ICGC_OV_AU_cancer_exp_seq_TPM.txt",
            sep = '\t',quote = F, row.names = F) 


##############################################################
## orginazing survival data
## read in donor.OV-AU 
donor_file <- read.delim(file = "donor.OV-AU.tsv",
                            sep = "\t",header = T,fill = T)
donor_surv <- donor_file %>% select(1,17,6) 

Psolid_donor_file <- read.delim(file = "ICGC_OV_AU_primary_solid_info.txt",
                                sep = "\t",header = T,fill = T) 
Psolid_donor_info <- Psolid_donor_file %>% select(4,6,9,11) 

Psolid_donor_surv <- merge(donor_surv, Psolid_donor_info, by = "icgc_donor_id", all = F)
if(sum(is.na(Psolid_donor_surv$donor_survival_time)) > 0){
  Psolid_donor_surv <- Psolid_donor_surv[-which(is.na(Psolid_donor_surv$donor_survival_time)),]
}

write.table(Psolid_donor_surv,file="ICGC_OV_AU_Psolid_surv.txt",sep = '\t', 
            quote=F,row.names = F)

## merge expression and survival data
Psolid_donor_surv <- read.table(file = "ICGC_OV_AU_Psolid_surv.txt",
                                header = T,fill = T,sep = "\t", check.names = F)
exp_tpm_cancer <- read.table(file = "ICGC_OV_AU_cancer_exp_seq_TPM.txt",
                             header = T,fill = T,check.names = F)
exp_tpm_cancer[1:5,1:3]

exp_tpm_cancer_t <- t(exp_tpm_cancer[,-1])
colnames(exp_tpm_cancer_t) <- exp_tpm_cancer$gene_id
exp_tpm_cancer_t[1:5,1:5]
cancer_donor_id <- sapply(rownames(exp_tpm_cancer_t),function(x) str_split(x,pattern = '-')[[1]][2])
exp_tpm_cancer_t_frame <- data.frame("icgc_donor_id" = cancer_donor_id, 
                                     "exp_sample" = rownames(exp_tpm_cancer_t),
                                     exp_tpm_cancer_t, check.names = F)
exp_tpm_cancer_t_frame[1:5,1:4]

exp_tpm_cancer_t_frame_surv <- merge(Psolid_donor_surv, exp_tpm_cancer_t_frame,
                                     by = "icgc_donor_id", all = F)
exp_tpm_cancer_t_frame_surv[1:4,1:10]

write.table(exp_tpm_cancer_t_frame_surv, file = "ICGC_OV_AU_Psolid_surv_TPM_all.txt",
            quote = F, sep = "\t", row.names =F)

## delete sample without survival data
path_icgc <- "F:\\Laboratory\\OvarianCancer_MelacularClassification\\OvarianCancer_ICGC\\OV_AU\\"
exp_tpm_cancer_t_frame_surv <- fread(file = str_c(path_icgc,"ICGC_OV_AU_Psolid_survall_TPM.txt"),
                                 sep = "\t", header =T, stringsAsFactors = F, data.table = T, check.names = F)
exp_tpm_cancer_t_frame_surv[1:3,1:3]
if(sum(is.na(exp_tpm_cancer_t_frame_surv$donor_survival_time)) > 0){
  exp_tpm_cancer_t_frame_surv_noNA <- exp_tpm_cancer_t_frame_surv[-which(is.na(exp_tpm_cancer_t_frame_surv$donor_survival_time)),]
}else{
  exp_tpm_cancer_t_frame_surv_noNA <- exp_tpm_cancer_t_frame_surv
}

write.table(exp_tpm_cancer_t_frame_surv_noNA, file = str_c(path_icgc,"ICGC_OV_AU_Psolid_survnoNA_TPM.txt"),
            quote = F, sep = "\t", row.names =F)

## log(TPM+1) transform
exp_tpm_cancer_t_frame_surv_noNA <- read.table(file = str_c(path_icgc,"ICGC_OV_AU_Psolid_survnoNA_TPM.txt"),
                                           header = T, fill = T, sep = "\t", check.names = F)

exp_tpm_cancer_t_frame_surv_noNA[1:4,1:10]
exp_tpm_cancer_t_frame_surv_noNA_log <- exp_tpm_cancer_t_frame_surv_noNA
exp_tpm_cancer_t_frame_surv_noNA_log[,-c(1:7)] <- log2(exp_tpm_cancer_t_frame_surv_noNA_log[,-c(1:7)]+1)
exp_tpm_cancer_t_frame_surv_noNA_log[1:4,c(8:12)]
which(exp_tpm_cancer_t_frame_surv_noNA_log<0)
# integer(0)
write.table(exp_tpm_cancer_t_frame_surv_noNA_log, file = str_c(path_icgc,"ICGC_OV_AU_Psolid_survnoNA_TPM_log2.txt"),
            quote = F, sep = "\t", row.names =F)

save.image(str_c(path,"1_ICGC_OV_FPKMtoTPM_log2.R"))



