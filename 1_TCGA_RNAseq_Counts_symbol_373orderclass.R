## ------ TCGAbiolinks download Counts data，symbol transform ------ ##

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
                  workflow.type = "HTSeq - Counts",
                  legacy = F) 
GDCdownload(query,
            directory = str_c(path,"GDCdata\\"), 
            method="api") 

expdat <- GDCprepare(query) 

expdat_frame <- assay(expdat) 
expdat_frame[1:5,1:5]
expdat_df <- data.frame("gene_id"=rownames(expdat_frame), expdat_frame, check.names = F)
expdat_df[1:4,1:4]
write.table(expdat_df, file = str_c(PRO,"_HTSeq_Counts.txt"),
            sep = "\t", quote = F, row.names = F) 


## --- symbol transform ------ ##
## read in transformed data from g:profiler website, merge
gProfiler <- read.table(file = str_c(path,"gProfiler_hsapiens_2021-5-9下午6-51-01.txt"),
                        sep = ",", header = T, fill = T, check.names = F)
gProfiler[1:4,1:4]
colnames(gProfiler)[c(1,3)] <- c("gene_id","gene_symbol")

counts_symbol <- merge(gProfiler[,c(1,3)], expdat_df, by = "gene_id", all = F)
counts_symbol2 <- counts_symbol[-which(counts_symbol$gene_symbol=="nan"),]

write.table(counts_symbol2,file = str_c(PRO,"_HTSeq_Counts_idsymbol.txt"),
            quote = F, sep = "\t", row.names = F)

length(which(duplicated(counts_symbol2$gene_symbol)))
## delete duplicated gene symbol, retain max value
if(length(which(duplicated(counts_symbol2$gene_symbol))) > 0){
  counts_symbol3 <- aggregate(x = counts_symbol2, by = list(counts_symbol2$gene_symbol), FUN = max)
}else{
  counts_symbol3 <- data.frame("Group.1" = counts_symbol2$gene_symbol, counts_symbol2, check.names = F)
}
which(duplicated(counts_symbol3$gene_symbol))
rownames(counts_symbol3) <- counts_symbol3$gene_symbol

write.table(counts_symbol3[,-1],file = str_c(PRO,"_HTSeq_Counts_idsymbol_noNAdup.txt"),
            quote = F, sep = "\t", row.names = F)

## process sample name
expdat_counts <- counts_symbol3[,-c(1:3)]
expdat_counts[1:4,1:4]
sample <- colnames(expdat_counts)
sample2 <- sapply(sample,function(x) str_c(unlist(str_split(x,pattern = '-'))[1:4],
                                                 collapse = "-")) 
sample2[1:4]
colnames(expdat_counts) <- sample2
expdat_counts[1:5,1:3]
expdat_counts_df <- data.frame(counts_symbol3[,2:3], expdat_counts, check.names = F)
expdat_counts_df[1:4,1:4]

write.table(expdat_counts_df, file=str_c(PRO,"_HTSeq_Counts_idsymbol_noNAdup_samch.txt"),
            sep = "\t", quote = F, row.names = F) 

save.image(str_c(path,"1_TCGA_RNAseq_Counts_symbol_373orderclass.RData"))
