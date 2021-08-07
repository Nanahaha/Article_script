library(NMF)
library(tidyverse) 
library(Rtsne)
library(RColorBrewer)  

path <- getwd()

data <- read.table(file = str_c(path,"TCGA_metabolism_MADcox_109genes.txt"),
                   sep = "\t", header = T, fill = T, check.names = F)
data[1:5,1:5]

tcga_class <- read.table(file = str_c(path,"k3\\TCGA_consensus_k3_hclust_subclass2.txt"),
                       sep = "\t", header = T, fill = T, check.names = F)
tcga_class[1:5,]

data_class <- merge(tcga_class, data, by = "Patients", all = F)
data_class[1:5,1:5]
str(data_class)
data_class$Class <- factor(data_class$Class, 
                           levels = c("1","2","3"),
                           ordered = T)
table(data_class$Class)
write.table(data_class, file = str_c(path,"TCGA_metabolism_MADcox_109genes_class.txt"),
            sep = "\t", quote = F, row.names = F)

##
set.seed(520)
tsne_out <- Rtsne(as.matrix(data_class[,-c(1:2)]),
                  perplexity=120, verbose=TRUE)

#######################################################################
## drawing
tsne_point <- data.frame(data_class$Class, tsne_out$Y)
colnames(tsne_point) <- c("Class","Coordinate 1", "Coordinate 2")

pdf(file = str_c(path,"k3\\TCGA_k3_tSNE.pdf"),
    width = 3.5, height = 3.5)
p <- ggplot(tsne_point, aes(x = `Coordinate 1`, y = `Coordinate 2`)) +
  geom_point(aes(colour = Class), 
             size = 2, alpha = 0.9) + 
  scale_colour_manual(breaks = c("1", "2", "3"),
                      labels = c("C1", "C2", "C3"),
                      values = brewer.pal(3,"Dark2")[1:3]) +
  theme(panel.grid.major = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_text(hjust = 0.5),
        legend.title=element_blank(),
        legend.position = c(.01, .98),
        legend.justification = c("left", "top"),
        legend.box.just = "left",
        legend.margin = margin(6, 6, 6, 6)) 
print(p)
dev.off()  

save.image(str_c(path,"7_TCGA_k3_tSNE.RData"))
