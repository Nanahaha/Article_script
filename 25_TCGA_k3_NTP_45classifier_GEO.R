## ------ nearest template prediction for GEO dataset ------ ##


library(CMScaller)
library(stringr)
library(data.table)

SETS <- "GSE32062"
PATH <- getwd()
path <- str_c(PATH,SETS,"\\")
data <- as.data.frame(fread(file = str_c(path,SETS,"_rma_quatile_prob_symbol_exprSet.txt"),
                            sep = "\t", header = T, fill = T, check.names = F), check.names = F)
data[1:4,1:4]  
rownames(data) <- data$GENE_SYMBOL
data2 <- t(scale(t(data[,-c(1:2)]), center = T, scale = T))  
data2[1:4,1:4]

path2 <- str_c(PATH,"k3\\")
classifier <- read.table(file = str_c(path2,"Classifier_45genes.txt"), header = T, 
                         fill = T, sep = "\t", check.names = F)
classifier[1:4,]

## nearest template prediction
res <- ntp(emat = data2, templates = classifier, doPlot = T)
class_pre <- data.frame("Patients" = rownames(res), "Class" = res$prediction)
table(class_pre$Class)
#C1  C2  C3 
#62  86 112  
write.table(class_pre, file = str_c(path,"k3\\GSE32062_k3_NTP_45classifier_subclass.txt"),
            sep = "\t", quote = F, row.names = F)

save.image(str_c(PATH,"25_TCGA_k3_NTP_51classifier_GEO.RData"))
