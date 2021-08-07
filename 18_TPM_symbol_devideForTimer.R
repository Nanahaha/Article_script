##  spliting TPM data for Timer

library(stringr)
library(data.table)

path <- getwd()
TPM_data <- as.data.frame(fread(file = str_c(path,"TCGA-OV_HTSeq_TPM_01survnoNA_rowsamp.txt"),
                                header = T, fill = T, check.names = F),
                          check.names = F) 
TPM_data[1:4,1:10]
rownames(TPM_data) <- TPM_data$Patients
TPM_data_t <- t(TPM_data[,-c(1:8)]) 
TPM_data_t[1:4,1:4]
TPM_data_t2 <- data.frame("gene_symbol" = rownames(TPM_data_t), TPM_data_t, check.names = F)
TPM_data_t2[1:4,1:4]


#########################################################################################
## spliting data for Timer, genes as row names and samples as column names
dim(TPM_data_t)
#[1] 38292   373

i <- 0
while(i < 301){
  write.table(TPM_data_t[,(i+1):(i+50)],file = str_c(path,"k3\\TIMER2_inputFile\\TIMER_matrix",i+1,"_",i+50,".txt"),
              sep = "\t", quote = F, row.names = T)
  i <- i+50
  
}
write.table(TPM_data_t[,351:ncol(TPM_data_t)],file = str_c(path,"k3\\TIMER2_inputFile\\TIMER_matrix351_",ncol(TPM_data_t),".txt"),
            sep = "\t", quote = F, row.names = T)

save.image(str_c(path,"18_TPM_symbol_devideForTimer.RData"))
