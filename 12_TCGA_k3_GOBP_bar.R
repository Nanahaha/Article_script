## ------ Mapping GO BP------ ## 

library(ggplot2)
library(stringr)
library(RColorBrewer)  # 颜色搭配条

path <- getwd()

filelist <- list.files(path = str_c(path,"k3\\diff_analysis\\DAVID_fc15fdr01\\"),
                       pattern = ".txt",
                       full.names = T)
myfile <- lapply(filelist[1:3], function(x) read.table(x, header = T, fill = T, 
                                                  sep = "\t", check.names = F))

## --- Plotting GO enrichment analysis bar graph --- ##
class <- c("C1", "C2", "C3")

for (i in 1:length(myfile)) {

  a <- myfile[[i]]
  
  a$GO_name <- sapply(a$Term, function(x) str_split(x, pattern = "~")[[1]][2])
  
  a2 <- a
  

  ## print have y axis text
  pdf(file = str_c(path,"k3\\TCGA_",class[i],"_BP_bar.pdf"),
      width = 7, height = 3.5)
  plot <- ggplot(a2[order(a2$PValue, decreasing = F),][1:10,],  
                 aes(x = reorder(GO_name,Count), y = Count, fill = -1*log10(FDR))) + 
    geom_bar(stat = "identity", width = 0.8) +
    coord_flip() + 
    scale_fill_gradient(expression(-log["10"](FDR)), low = "#0072B5", high = "#BC3C28")+  
    labs(y = "Gene Count",
         x = "") + 
    theme_bw() + 
    theme(panel.grid = element_blank())
  print(plot)
  dev.off()
  
  ## print no y axis text
  pdf(file = str_c(path,"k3\\TCGA_",class[i],"_BP_barNoYtext.pdf"),
      width = 3.5, height = 3.5)
  plot <- ggplot(a2[order(a2$PValue, decreasing = F),][1:10,],  
                 aes(x = reorder(GO_name,Count), y = Count, fill = -1*log10(FDR))) + 
    geom_bar(stat = "identity", width = 0.8) +
    coord_flip() +  
    scale_fill_gradient(expression(-log["10"](FDR)), low = "#0072B5", high = "#BC3C28")+  
    labs(y = "Gene Count",
         x = "") + 
    theme_bw() + 
    theme(panel.grid = element_blank(),
          axis.text.y = element_blank())
  print(plot)
  dev.off()

}




