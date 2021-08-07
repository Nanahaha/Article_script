## ------ filter metabolism-related expression matrix: (MAD) >0.5 & single factor cox p-value<0.05 ----- ##

library(tidyverse) 
library(stats) 
library(survival)
library(data.table) 

path <- getwd()

data <- read.table(file = str_c(path,"TCGA_metabolism_GeneMatrix.txt"),
                   header = T, fill = T, check.names = F)
dim(data)
# 2718  374
data[1:5,1:3]
rownames(data) <- data$gene_symbol
data_rowsamp <- t(data[,-1])
data_rowsamp[1:4,1:4]

###################################### median absolute deviation (MAD) >0.5
filter1 <- c()
for (i in 1:ncol(data_rowsamp)) {
  if(mad(data_rowsamp[,i])<=0.5){
    filter1 <- c(filter1,i)
  }
}
data3 <- data.frame("Patients" = rownames(data_rowsamp), data_rowsamp[,-filter1], check.names = F)

################################################################## single factor cox
tcga_os_all <- read.table(file = str_c(path,"TCGA-OV_surv373.txt"),
                               sep = "\t", header = T, fill = T, check.names = F)
data_os <- merge(tcga_os_all, data3, by = "Patients", all =F)
data_os[1:4,1:9]
dim(data_os)

colnames(data_os)[c(4:5)] <- c("Survival_Status","OS_Time")
table(data_os$Survival_Status)
# Alive  Dead 
# 143   230

data_os$Survival_Status[which(data_os$Survival_Status=="Alive")] <- 0
data_os$Survival_Status[which(data_os$Survival_Status=="Dead")] <- 1
data_os$Survival_Status <- as.numeric(data_os$Survival_Status)
table(data_os$Survival_Status)
#  0   1 
# 143 230 
str(data_os)
data_os[1:4,1:9]

covariates <- colnames(data_os)[-c(1:8)]  
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste("Surv(OS_Time, Survival_Status)~", x))) 

univ_models <- lapply( univ_formulas, function(x){coxph(x, data = data_os)})
# Extract data 
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2); #coeficient beta
                         HR <-signif(x$coef[2], digits=2); #exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, wald.test, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                       "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res <- t(as.data.frame(univ_results, check.names = FALSE))
res2 <- as.data.frame(res)
res2[1:3,]

filter2 <- rownames(res2)[which(res2$p.value>=0.05)]
length(filter2)
# 1821

data4 <- data3 %>% select(-filter2)



write.table(data4, file = str_c(path,"TCGA_metabolism_MADcox_118genes.txt"),
            sep = "\t", quote = F, row.names = F)
# delete "ARSH"    "GALNT15" "GPAT4"   "KYAT1"   "NT5C3A"  "PLPP5"   "SDHAF3"  "SLC51A" gene manually, 
# because this eight gene doesn't exist in GEO dataset.
# delete "SRXN1" gene manually, because this gene doesn't exist in ICGC dataset.
# retain 109 metabolism-related genes finally 

save.image(str_c(path,"3_TCGA_metabolismGenes_filter.RData"))

