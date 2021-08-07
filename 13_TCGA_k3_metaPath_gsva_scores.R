## ------ score metabolism related pathway using GSVA ------ ##

library(stringr)
library(GSVA)
library(pheatmap)
library(data.table)

path <- getwd()
data <- data.frame(fread(file = str_c(path,"TCGA-OV_HTSeq_TPM_01survnoNA_rowsamp_log2.txt"),
                         header = T, fill = T, check.names = F), check.names = F) # 带生存数据的
data[1:4,1:10]

rownames(data) <- data$Patients

# creat gene set list
signatures <- file(str_c(path,"metabolism114_signatures_genes.txt"),open = "r")
signatures_list <- list()
n <- 1
while (TRUE) {
  line <- readLines(signatures,n = 1)
  if(length(line) == 0){
    break
  }
  del2 <- unlist(str_split(line,"\t"))
  del3 <- del2[-which(del2=="")]
  signatures_list[[del3[1]]] <- del3[-1]
  n <- n+1
}
close(signatures)

expr_matrix <- t(data[,-c(1:8)])
expr_matrix[1:4,1:4]
es.max <- gsva(expr_matrix, signatures_list, min.sz = 1) 
pheatmap(es.max)
es.max[1:4,1:4]

write.table(es.max, file = str_c(path,"k3\\TCGA_k3_metaPath_gsva_scores.txt"),
            sep = "\t", quote = F, row.names = T)

save.image(str_c(path,"13_TCGA_k3_metaPath_gsva_scores.RData"))



