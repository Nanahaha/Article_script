library(tidyverse)
library(survival)
library(survminer)
library(stringr)

path <- getwd()
allsample_os <- read.table(file = str_c(path,"ICGC_OV_AU_surv81.txt"),
                      sep = "\t", header = T)
class <- read.table(file = str_c(path,"k3\\ICGC_consensus_k3_hclust_subclass_reclass.txt"),
                    sep = "\t", header = T)

class_os <- merge(class, allsample_os, by = "Patients", all = F)

table(class_os$Survival_Status)
# Alive  Dead 
#15    66
table(class_os$Class)
#  1  2  3 
# 27 34 20  

class_os$Survival_Status[which(class_os$Survival_Status=="Dead")] <- 1
class_os$Survival_Status[class_os$Survival_Status!=1] <- 0
table(class_os$Survival_Status)
# 0   1 
# 15 66
class_os[1:4,]

############### 
## survival plot
fit <- survfit(Surv(class_os$OS_Time, as.numeric(class_os$Survival_Status)) ~ class_os$Class, 
               data = class_os)

pdf(file = str_c(path,"k3\\ICGC_k3_reclass_survival.pdf"),
    width = 4, height = 4)
p1 <- ggsurvplot(fit,
                 pval = TRUE, pval.size = 4.5,
                 pval.coord = c(0.15,0.15), 
                 surv.median.line = "hv", 
                 legend.title = " ",
                 legend.labs = c("C1", "C2", "C3"),
                 legend = c(0.9,0.75), 
                 ylab = "Overall Survival",
                 xlab = "Time (Days)",
                 censor.shape = 124,censor.size = 2, 
                 conf.int = FALSE,
                 ggtheme = theme_bw(), 
                 risk.table = F,
                 palette = "npg",
                 font.y = c(14)
                 )
p1$plot
dev.off()


save.image(str_c(path,"9_ICGC_109genes_k3_reclass_survival.RData"))
