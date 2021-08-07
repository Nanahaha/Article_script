## ------ the mutation rate of BRCA1 and BRCA2 in each subtype ------ ##

library(stringr)
library(reshape2)  
library(ggplot2)
library(RColorBrewer)  
library(tableone)   

path <- getwd()
snv_df <- read.table(file = str_c(path,"TCGA-OV.mutect2_snv_num_screen_class_df.txt"),
                      header = T, fill = T, check.names = F)
snv_df[1:4,1:4]
min(snv_df[,-c(1:2)])
max(snv_df[,-c(1:2)])
which(colnames(snv_df) == "BRCA1")

snv_df2 <- as.matrix(snv_df[,-c(1:2)])
snv_df2[which(snv_df2 > 0)] <- 1

BRCA_snv <- cbind(snv_df[,1:2], snv_df2[,c("BRCA1","BRCA2")])
BRCA_snv$Class <- factor(BRCA_snv$Class, 
                         levels = c("C1","C2","C3"),
                         ordered = T)
BRCA_snv$BRCA1 <- factor(BRCA_snv$BRCA1, 
                         levels = c(0,1),
                         ordered = T)
BRCA_snv$BRCA2 <- factor(BRCA_snv$BRCA2, 
                         levels = c(0,1),
                         ordered = T)
BRCA_snv[1:4,1:4]
table(BRCA_snv[,c(2,3)])

tab1 <- CreateTableOne(vars = colnames(BRCA_snv)[c(3:4)], data = BRCA_snv, 
                       factorVars = colnames(BRCA_snv)[c(3:4)],
                       strata="Class")
tab11 <- print(tab1, printToggle = FALSE, noSpaces = TRUE)
write.table(tab11, file = str_c(path, "k3\\BRCA_class_pvalue.txt"),
            quote = F, sep = "\t")

BRCA_C1C2 <- BRCA_snv[-which(BRCA_snv$Class=="C3"),]
BRCA_C1C2$Class <- factor(BRCA_C1C2$Class, levels = c("C1","C2"), ordered = T)
tabC1C2 <- CreateTableOne(vars = colnames(BRCA_C1C2)[c(3:4)], data = BRCA_C1C2, 
                       factorVars = colnames(BRCA_C1C2)[c(3:4)],
                       strata="Class")
tabC1C211 <- print(tabC1C2, printToggle = FALSE, noSpaces = TRUE)
write.table(tabC1C211, file = str_c(path, "k3\\BRCA_C1C2_pvalue.txt"),
            quote = F, sep = "\t")
# C1和C3
BRCA_C1C3 <- BRCA_snv[-which(BRCA_snv$Class=="C2"),]
BRCA_C1C3$Class <- factor(BRCA_C1C3$Class, levels = c("C1","C3"), ordered = T)
tabC1C3 <- CreateTableOne(vars = colnames(BRCA_C1C3)[c(3:4)], data = BRCA_C1C3, 
                          factorVars = colnames(BRCA_C1C3)[c(3:4)],
                          strata="Class")
tabC1C311 <- print(tabC1C3, printToggle = FALSE, noSpaces = TRUE)
write.table(tabC1C311, file = str_c(path, "k3\\BRCA_C1C3_pvalue.txt"),
            quote = F, sep = "\t")

BRCA_C2C3 <- BRCA_snv[-which(BRCA_snv$Class=="C1"),]
BRCA_C2C3$Class <- factor(BRCA_C2C3$Class, levels = c("C2","C3"), ordered = T)
tabC2C3 <- CreateTableOne(vars = colnames(BRCA_C2C3)[c(3:4)], data = BRCA_C2C3, 
                          factorVars = colnames(BRCA_C2C3)[c(3:4)],
                          strata="Class")
tabC2C311 <- print(tabC2C3, printToggle = FALSE, noSpaces = TRUE)
write.table(tabC2C311, file = str_c(path, "k3\\BRCA_C2C3_pvalue.txt"),
            quote = F, sep = "\t")


for (i in 3:4) {
  mytable <- table(BRCA_snv[,c(2,i)])
  mytable_fre <- prop.table(mytable,1)  
  mytable_fre2 <- data.frame(mytable_fre)
  mytable_fre2$gene_sta <- NA
  mytable_fre2$gene_sta[which(mytable_fre2[,2] == 0)] <- "WT"
  mytable_fre2$gene_sta[which(mytable_fre2[,2] == 1)] <- "Mut"
  mytable_fre2$gene_sta <- factor(mytable_fre2$gene_sta, 
                                   levels = c("WT","Mut"),
                                   ordered = T)
  dev.size()
  pdf(file = str_c(path,"k3\\",colnames(mytable_fre2)[2],"_class_freq_bar.pdf"),
      width = 2.7, height = 3)
  p <- ggplot(mytable_fre2, aes(x = Class, y = Freq, fill = gene_sta)) +
    geom_bar(stat = "identity") + 
    scale_fill_manual(values = brewer.pal(8,"Set2")[c(8,2)]) +
    labs(y = "Ratio") + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.title.x = element_blank(),
          legend.title = element_blank(),
          legend.position = c(0.15,0.87),
          legend.text = element_text(size = 8),  
          legend.key.width = unit(3.5, "mm") ,
          legend.key.height = unit(3.5, "mm")  
    )
  print(p)
  dev.off()
}

save.image(str_c(path,"22_TCGA_k3_BRCA_freq_p_bar.RData"))

