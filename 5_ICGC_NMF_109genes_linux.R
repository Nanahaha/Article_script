## linux
cd /data/users/liuxiaona/Laboratory/OvarianCancer_MelacularClassification/OvarianCancer_ICGC/OV_AU
R

library(NMF)
library(stringr) 

path <- getwd()
data <- read.table(file = str_c(path,"/ICGC_metabolism_colsamp_109genes.txt"),
                   sep = "\t", header = T, fill = T, check.names = F)
data[1:5,1:4]
rownames(data) <- data$gene_symbol
dim(data)
which(is.na(data))
which(data<0)

estim.r <- nmf(data[,-1], 2:10, nrun = 1000, seed = 123456)

save.image(str_c(path,"/5_ICGC_NMF_109genes_linux.RData"))

pdf(file = str_c(path,"/ICGC_NMF_cophenetic_109genes.pdf"),
    width = 10,height = 8)
plot(estim.r)
dev.off()

## heatmap
pdf(file = str_c(path,"/ICGC_NMF_heatmap_109genes.pdf"),
    width = 18,height = 15)
consensusmap(estim.r, annCol=NA,labCol=NA, labRow=NA)
dev.off()

