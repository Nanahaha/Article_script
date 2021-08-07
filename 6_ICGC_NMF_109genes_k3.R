## linux
cd /data/users/liuxiaona/Laboratory/OvarianCancer_MelacularClassification/OvarianCancer_ICGC
R

library(NMF)
library(tidyverse) 
library(cluster)
library(stringr)

path <- getwd()
data <- read.table(file = str_c(path,"/ICGC_metabolism_colsamp_109genes.txt"),
                   sep = "\t", header = T, fill = T, check.names = F)
data[1:5,1:4]
rownames(data) <- data$gene_symbol
dim(data)
which(is.na(data))
which(data<0)
data[1:5,1:4]

###########################################################################
meth <- nmfAlgorithm(version = "R")
meth <- c(names(meth), meth)
length(meth)
res.multi.method <- nmf(data[,-1], 3, nrun = 200, meth, 
                        seed = 123456, .options="t")
png(filename = str_c(path,"k3_nmf_method_selection.png"))
plot(res.multi.method)  
dev.off()
## ".R#offset" method is best

res <- nmf(data[,-1], 3, nrun = 1000, method = ".R#offset", 
           seed = 123456, .options = "t")

save.image(str_c(path,"6_ICGC_NMF_109genes_k3.RData"))

pdf(file = str_c(path,"/k3/ICGC_NMF_109genes_k3_consensusmap.pdf"),
    width = 10, height = 10)
consensusmap(res)
dev.off()

## hclust extract subtype
conse <- consensus(res)
result <- dist(conse, method = "euclidean")
result_hc <- hclust(d = result)
class_result <- cutree(result_hc, k = 3) # https://www.r-bloggers.com/lang/chinese/619
sample_class <- data.frame("Patients" =  colnames(data[,-1]), "Class" = class_result)
table(sample_class$Class)
# 1  2  3
#34 27 20

write.table(sample_class, file = str_c(path,"/k3/ICGC_consensus_k3_hclust_subclass.txt"),
            quote = F, sep = "\t", row.names = F)






