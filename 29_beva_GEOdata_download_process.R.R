## ------- GEOquery download bevacizumab expression profile ------ ##

library(GEOquery)
library(stringr)
library(dplyr)
library(data.table)


SET <- "GSE140082"
PATH <- getwd()
path <- str_c(PATH,"bevacizumab_exp_submap\\")
gse <- getGEO(SET, GSEMatrix = T, destdir = path)

metdata <- pData(gse[[1]])
write.table(metdata, file = str_c(path,SET,"_phenodata.txt"),
            sep = "\t", quote = F, row.names = F)

expma <- exprs(gse[[1]])
expma[1:4,1:4]
expma2 <- data.frame("ID"=rownames(expma), expma, check.names = F)
expma2[1:4,1:4]
write.table(expma2, file = str_c(path,SET,"_exp_matrix.txt"),
            sep = "\t", quote = F, row.names = F)

GPL <- "GPL14951"
gpl <- getGEO(GPL, destdir = path) %>%
  Table() 
save(gpl, file = str_c(path,GPL,"_annot.RData"))
gpl[1:4,]
colnames(gpl)
gpl$Symbol

probe <- gpl %>% 
  select("ID","Symbol","Entrez_Gene_ID")

expma2_sym <- merge(probe[,1:2], expma2, by = "ID", all = F)
expma2_sym[1:4,1:4]
expma2_sym$Symbol

if(length(which(duplicated(expma2_sym$Symbol)))>0){
  expma2_sym2 <- aggregate(x = expma2_sym, by = list(expma2_sym$Symbol), FUN = max)
}else{
  expma2_sym2 <- data.frame("Group.1" = expma2_sym$Symbol, expma2_sym, check.names = F)
}
expma2_sym2[1:4,1:4]
which(duplicated(expma2_sym2$Symbol))

write.table(expma2_sym2[,-1], file = str_c(path,SET,"_exp_matrix_idsym_nodup.txt"),
            sep = "\t", quote = F, row.names = F)

###########################################################################

metdata_beva <- metdata %>% filter(characteristics_ch1.2=="treatment: bevacizumab") %>%
  select("geo_accession","characteristics_ch1.2","characteristics_ch1.4")
table(metdata_beva$characteristics_ch1.4)

metdata_beva2 <- metdata_beva[-which(metdata_beva$characteristics_ch1.4=="debulking_status: Inoperable"),]

metdata_beva2$response <- NA
metdata_beva2$response[which(metdata_beva2$characteristics_ch1.4=="debulking_status: OPTIMAL")] <- "Optimal"
metdata_beva2$response[which(metdata_beva2$characteristics_ch1.4=="debulking_status: SUB-OPTIMAL")] <- "Suboptimal"
metdata_beva2[1:4,]

rownames(expma2_sym2) <- expma2_sym2$Symbol
expma2_sym2[1:4,1:4]
expma2_sym2_t <- t(expma2_sym2[,-c(1:3)])
expma2_sym2_t[1:4,1:4]
expma2_sym2_t_df <- data.frame("geo_accession"=rownames(expma2_sym2_t), expma2_sym2_t,
                               check.names = F)
expma2_sym2_t_df[1:4,1:4]
# merge
expma2_sym2_t_df_beva <- merge(metdata_beva2, expma2_sym2_t_df, 
                               by = "geo_accession", all = F)
expma2_sym2_t_df_beva[1:4,1:7]
dim(expma2_sym2_t_df_beva)
# 198 20823

write.table(expma2_sym2_t_df_beva[,1:4], file = str_c(path,SET,"_phenodata_198beva.txt"),
            sep = "\t", quote = F, row.names = F)

forcls0 <- expma2_sym2_t_df_beva[,c(1,4)]
forcls0$response_num <- NA
forcls0$response_num[which(forcls0$response=="Suboptimal")] <- 1
forcls0$response_num[which(forcls0$response=="Optimal")] <- 2
forcls0[1:4,]
forcls <- t(forcls0[,c(1,3)])
forcls[,1:4]
write.table(forcls, file = str_c(path,SET,"_phenodata_198beva_space.txt"),
            sep = " ", quote = F, row.names = F, col.names = F)

expma2_sym2_t_df_beva[1:4,1:6]
rownames(expma2_sym2_t_df_beva) <- expma2_sym2_t_df_beva$geo_accession
exp_beva_rowgene <- t(expma2_sym2_t_df_beva[,-c(1:4)])
exp_beva_rowgene[1:4,1:4]
exp_beva_rowgene_df <- data.frame("gene_symbol"=rownames(exp_beva_rowgene), exp_beva_rowgene,
                                  check.names = F)
exp_beva_rowgene_df[1:4,1:4]
write.table(exp_beva_rowgene_df, file = str_c(path,SET,"_exp_matrix_sym_nodup_colsamp_beva.txt"),
            sep = "\t", quote = F, row.names = F)


## --- Find genes shared with TCGA ---  ##
path2 <- getwd()
tcga_exp <- data.frame(fread(file = str_c(path2,"TCGA-OV_HTSeq_TPM_01survnoNA_log2_colsamp.txt"),
                             header = T, fill = T, check.names = F), check.names = F)
tcga_exp[1:4,1:4]

exp_beva_rowgene_df_tcga <- exp_beva_rowgene_df[which(exp_beva_rowgene_df$gene_symbol %in% tcga_exp$gene_symbol),]
dim(exp_beva_rowgene_df_tcga)
# 16951   199
write.table(exp_beva_rowgene_df_tcga, file = str_c(path,SET,"_exp_matrix_sym_nodup_colsamp_beva_tcga.txt"),
            sep = "\t", quote = F, row.names = F)

tcga_exp_geobeva <- tcga_exp[which(tcga_exp$gene_symbol %in% exp_beva_rowgene_df$gene_symbol),]
dim(tcga_exp_geobeva)
# 16951   374
write.table(tcga_exp_geobeva, file = str_c(path2,"TCGA-OV_HTSeq_TPM_01survnoNA_log2_colsamp_beva.txt"),
            sep = "\t", quote = F, row.names = F)


save.image(file = str_c(path,"29_beva_GEOdata_download_process.RData"))
