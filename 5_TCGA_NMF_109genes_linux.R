## linux
cd /data/users/liuxiaona/Laboratory/OvarianCancer_MelacularClassification/OvarianCancer_TCGA
R

library(NMF)
library(stringr) 

path <- getwd()
data <- read.table(file = str_c(path,"/TCGA_metabolism_MADcox_109genes.txt"),
                   sep = "\t", header = T, fill = T, check.names = F)
data[1:5,1:4]
rownames(data) <- data$Patients
data_t <- t(data[,-1])
data_t[1:5,1:3]
which(data_t<0)
dim(data_t)
data_t[1:4,1:4]

estim.r <- nmf(data_t, 2:10, nrun = 1000, seed = 123456)

save.image(str_c(path,"/5_TCGA_NMF_109genes_linux.RData"))

pdf(file = str_c(path,"/TCGA_NMF_cophenetic_109genes.pdf"),
    width = 10,height = 8)
plot(estim.r)
dev.off()

## heatmap
pdf(file = str_c(path,"/TCGA_NMF_heatmap_109genes.pdf"),
    width = 18,height = 15)
consensusmap(estim.r, annCol=NA,labCol=NA, labRow=NA)
dev.off()

