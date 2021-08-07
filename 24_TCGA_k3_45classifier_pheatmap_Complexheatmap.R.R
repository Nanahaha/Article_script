## ------ heatmap of 51-genes classifier in TCGA dataset ------ ##

library(stringr)
library(data.table)
library(ggplot2)
library(RColorBrewer)  
library(pheatmap)
library(circlize)  
library(CMScaller)

path <- getwd()
data <- as.data.frame(fread(file = str_c(path,"TCGA-OV_HTSeq_TPM_01survnoNA_log2_colsamp.txt"),
                            sep = "\t", header = T, fill = T, check.names = F), check.names = F) 
data[1:4,1:4]  
rownames(data) <- data$gene_symbol

subclass <- read.table(file = str_c(path,"k3\\TCGA_consensus_k3_hclust_subclass2.txt"),
                       header = T, fill = T, check.names = F)
subclass$Class <- sapply(subclass$Class,function(x) str_c("C",x))
subclass[1:4,]

data_t <- t(data[,-1])
data_t_df <- data.frame("Patients" = rownames(data_t), data_t, check.names = F)
data_t_df[1:4,1:4]
data_class <- merge(subclass, data_t_df, by = "Patients", all = T)
table(data_class$Class)
# C1  C2  C3 
# 106 119 148 
data_class_order <- data_class[order(data_class$Class),]
data_class_order$Class
data_class_order[1:4,1:4]

classifier <- read.table(file = str_c(path,"k3\\Classifier_45genes.txt"), header = T, 
                         fill = T, sep = "\t", check.names = F)
classifier[1:4,]

classifier$class2 <- NA
classifier$class2[which(classifier$class=="C1")] <- "Specific genes of C1"
classifier$class2[which(classifier$class=="C2")] <- "Specific genes of C2"
classifier$class2[which(classifier$class=="C3")] <- "Specific genes of C3"
classifier[1:4,]

genes45_exp <- data_class_order[,c("Patients","Class",classifier$probe)]
genes45_exp[1:4,1:7]
dim(genes45_exp)
# 373  47 
rownames(genes45_exp) <- genes45_exp$Patients
genes45_exp_colsamp <- t(genes45_exp[,-c(1:2)])
genes45_exp_colsamp[1:4,1:4]
min(genes45_exp_colsamp)

genes45_exp_colsamp_scale <- t(scale(t(genes45_exp_colsamp)))
genes45_exp_colsamp_scale[1:4,1:4]
min(genes45_exp_colsamp_scale)
max(genes45_exp_colsamp_scale)

## ------ heatmap ------ ##
COLOR <- brewer.pal(7,"RdBu")

col_CLASS_COLOR <- brewer.pal(3,"Set2")

unique(annotation_row$Pathway_class)

col_colors <- list(
  Class = c(C1 = col_CLASS_COLOR[1],
            C2 = col_CLASS_COLOR[2],
            C3 = col_CLASS_COLOR[3])
)

row_colors <- list(
  class2 = c("Specific genes of C1" = col_CLASS_COLOR[1],
             "Specific genes of C2" = col_CLASS_COLOR[2],
             "Specific genes of C3" = col_CLASS_COLOR[3])
)

col_anno <- HeatmapAnnotation(Class = genes45_exp$Class,
                              col = col_colors,
                              which = "column",
                              show_legend = F,  # 不显示分类图例
                              show_annotation_name = F  #不显示class
)

row_anno <- HeatmapAnnotation(class2 = classifier$class2,
                              col = row_colors,
                              which = "row",
                              show_legend = F,  # 不显示分类图例
                              show_annotation_name = F  #不显示class
)

class_num <- as.vector(table(genes45_exp$Class))

pdf(file = str_c(path,"k3\\Classifier_45genes_TCGA_complexheatmap.pdf"),
    width = 8, height = 6)
heat <- Heatmap(genes45_exp_colsamp_scale,  
                col = rev(COLOR), 
                cluster_rows = F,
                show_row_dend = F,
                show_row_names = F,
                row_names_side = "left",
                row_names_gp = gpar(fontsize = 7),
                row_title_side = "left",
                row_title = unique(classifier$class2),
                row_split = classifier$class2,
                row_gap = unit(0, "mm"),
                row_title_gp = gpar(col = col_CLASS_COLOR),  
                row_title_ro = 90,  
                left_annotation = row_anno,  
                cluster_columns = F,  
                show_column_dend = F,
                show_column_names = F,
                column_split = genes45_exp$Class,
                column_title = c(str_c("C1 (n = ",class_num[1],")"),
                                 str_c("C2 (n = ",class_num[2],")"),
                                 str_c("C3 (n = ",class_num[3],")")), 
                column_title_side = "top",  
                column_title_gp = gpar(col = col_CLASS_COLOR),
                column_gap = unit(0, "mm"),
                top_annotation = col_anno 
)
print(heat)
dev.off()


annotation_col = data.frame(Class = genes45_exp$Class)
rownames(annotation_col) = genes45_exp$Patients

col_CLASS_COLOR <- brewer.pal(3,"Set2")

ann_colors = list(
  Class = c(C1 = col_CLASS_COLOR[1],
            C2 = col_CLASS_COLOR[2],
            C3 = col_CLASS_COLOR[3])
)

pdf(file = str_c(path,"k3\\Classifier_45genes_TCGA_pheatmap.pdf"),
    width = 8, height = 6)
pheatmap(genes45_exp_colsamp,  
         scale = "row", 
         cluster_cols = F,
         cluster_rows = F,
         color = rev(COLOR), 
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         show_colnames = F,
         show_rownames = F
)
dev.off()

save.image(str_c(path,"25_TCGA_k3_45classifier_pheatmap_Complexheatmap.RData"))
