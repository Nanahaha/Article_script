## ------ survival of predicted subtypes in GEO dataset  ------- ##

library(tidyverse)
library(survival)
library(survminer)

SETS <- "GSE32062" 
PATH <- getwd()
PATH1 <- str_c(PATH,SETS,"\\")

geo_os <- read.table(file = str_c(PATH1,SETS,"_surv.txt"), 
                      sep = '\t', header = T)
geo_class <- read.table(file = str_c(PATH1,"k3\\",SETS,"_k3_NTP_45classifier_subclass.txt"),
                         sep = '\t', header = T)
colnames(geo_os)[c(1,4:5)] <- c(colnames(geo_class)[1],"OS_Time","Survival_Status")

geo_class_os <- merge(geo_class, geo_os[,c(1,4:5)], by = "Patients", all = F)
geo_class_os2 <- na.omit(geo_class_os)
table(geo_class_os2$Survival_Status)
#  0   1 
# 139 121 
table(geo_class_os2$Class)
# C1  C2  C3 
# 62  86 112 

geo_class_os2$Class <- factor(geo_class_os2$Class,levels = c("C1", "C2", "C3"))
levels(geo_class_os2$Class)

geo_class_os2$Survival_Status[which(geo_class_os2$Survival_Status=="Dead")] <- 1
geo_class_os2$Survival_Status[geo_class_os2$Survival_Status!=1] <- 0


fit <- survfit(Surv(geo_class_os2$OS_Time, as.numeric(geo_class_os2$Survival_Status)) ~ geo_class_os2$Class, 
               data = geo_class_os2)

pdf(file = str_c(PATH1,"k3\\",SETS,"_NTP_k3_45classifier_survival.pdf"),  
    width = 4, height = 4)
p1 <- ggsurvplot(fit,
                 pval = TRUE, pval.size = 4.5,
                 pval.coord = c(0.15,0.15), 
                 surv.median.line = "hv", 
                 legend.title = " ",
                 legend.labs = c("C1", "C2", "C3"),  
                 legend = c(0.9,0.75), 
                 ylab = "Overall Survival",
                 xlab = "Time (Months)",
                 censor.shape = 124,censor.size = 2, 
                 conf.int = FALSE,
                 ggtheme = theme_bw(), 
                 risk.table = F,
                 palette = "npg",
                 font.y = c(14)
                 )
p1$plot
dev.off()

save.image(str_c(PATH,"26_",SETS,"_NTP_k3_45classifier_survival.RData"))  

