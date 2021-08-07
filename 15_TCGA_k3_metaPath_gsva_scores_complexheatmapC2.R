## ------ Drawing heatmap for GSVA score of C2 subtype ------ ##

library(ComplexHeatmap)
library(circlize)  
library(gplots)  
library(stringr)
library(ggplot2)
library(RColorBrewer)  

path <- getwd()

data <- read.table(file = str_c(path,"k3\\TCGA_k3_metaPath_gsva_scores.txt"),
                   header = T, fill = T, sep = "\t", check.names = F)
data[1:4,1:3]
subclass <- read.table(file = str_c(path,"k3\\TCGA_consensus_k3_hclust_subclass2.txt"),
                       header = T, fill = T)
pathway_class <- read.table(file = str_c(path,"metabolism114_signatures_genes_KEGGclass.txt"),
                            header = T, fill = T, sep = "\t", check.names =F)

subclass$Class <- sapply(subclass$Class,function(x) str_c("C",x))
table(subclass$Class)
#  C1  C2  C3 
# 106 119 148 

CLASS <- "C2"
c_spe_pathways <- read.table(file = str_c(path,"k3\\metaPath_gsva_diff_analysis\\","TCGA_",CLASS,"_fc0fdr05_gsva.txt"),
                             sep = "\t", header = T, fill = T)

c_pathways_scores <- data[which(rownames(data) %in% rownames(c_spe_pathways)),]
c_pathways_scores [1:4,1:3]

path_class_data <- data.frame("Pathway_signatures" = rownames(c_pathways_scores),
                              c_pathways_scores, check.names = F)
path_class_data[1:4,1:3]

which(!path_class_data$Pathway_signatures %in% pathway_class$Pathway_signatures)
path_class_data$Pathway_signatures[43]

path_class_data2 <- merge(pathway_class[,c(1,3)], path_class_data,
                          by = "Pathway_signatures", all=F)
path_class_data2_order  <- path_class_data2[order(path_class_data2$Pathway_TCGAclass2),]
path_class_data2_order[1:4,1:3]
rownames(path_class_data2_order) <- path_class_data2_order$Pathway_signatures
path_class_data2_order[1:4,1:3]
table(path_class_data2_order$Pathway_TCGAclass2)
#Amino acid metabolism              Carbohydrate metabolism 
#4                                   11 
#Drug metabolism                    Energy metabolism 
#3                                    1 
#Glycan biosynthesis and metabolism                     Lipid metabolism 
#5                                   17 
#Metabolism of cofactors and vitamins                     Other metabolism 
#8                                   10

temp1 <- t(path_class_data2_order[,-c(1:2)])
temp2 <- data.frame("Patients" = row.names(temp1), temp1, check.names = F)
temp_class <- merge(subclass,temp2, by = "Patients", all = F)
temp_class_order <- temp_class[order(temp_class$Class),]
temp_class_order[1:4,1:3]
rownames(temp_class_order) <- temp_class_order$Patients
temp_class_order[1:5,1:3]
temp_class_order_scale <- scale(temp_class_order[,-c(1:2)])
temp_class_order_scale[1:4,1:2]
c_spe_pathways_heatmap <- t(temp_class_order_scale)  
c_spe_pathways_heatmap[1:4,1:4]


# --- heatmap --- #
annotation_col = data.frame(Class = temp_class_order$Class)
rownames(annotation_col) = temp_class_order$Patients

annotation_row = data.frame(Pathway_class = path_class_data2_order$Pathway_TCGAclass2)
rownames(annotation_row) = path_class_data2_order$Pathway_signatures

COLOR <- brewer.pal(7,"RdBu")

# Specify colors
#display.brewer.pal(3,"Set2")
col_CLASS_COLOR <- brewer.pal(3,"Set2")

#display.brewer.pal(7,"Set3")
row_CLASS_COLOR <- brewer.pal(8,"Set3")

unique(annotation_row$Pathway_class)

col_colors <- list(
  Class = c(C1 = col_CLASS_COLOR[1],
            C2 = col_CLASS_COLOR[2],
            C3 = col_CLASS_COLOR[3])
)

row_colors <- list(
  Pathway_TCGAclass2 = c("Amino acid metabolism" = row_CLASS_COLOR[1],
                         "Carbohydrate metabolism" = row_CLASS_COLOR[2],
                         "Drug metabolism" = row_CLASS_COLOR[3],
                         "Energy metabolism" = row_CLASS_COLOR[4],
                         "Glycan biosynthesis and metabolism" = row_CLASS_COLOR[5],
                         "Lipid metabolism" = row_CLASS_COLOR[6],
                         "Metabolism of cofactors and vitamins" = row_CLASS_COLOR[7],
                         "Other metabolism" = row_CLASS_COLOR[8])
)

col_anno <- HeatmapAnnotation(Class = temp_class_order$Class,
                              col = col_colors,
                              which = "column",
                              show_legend = T,  
                              show_annotation_name = F  
)

row_anno <- HeatmapAnnotation(Pathway_TCGAclass2 = path_class_data2_order$Pathway_TCGAclass2,
                              col = row_colors,
                              which = "row",
                              show_legend = F,  
                              show_annotation_name = F  
)

class_num <- as.vector(table(temp_class_order$Class))


dev.size()
pdf(file = str_c(path,"k3\\",CLASS,"_fc0fdr05_metaPath_gsva_heatmap.pdf"),
    width = 12, height = 6)
heat <- Heatmap(c_spe_pathways_heatmap,  
                col = rev(COLOR),
                cluster_rows = T,
                show_row_dend = F,
                row_names_side = "left",
                row_names_gp = gpar(fontsize = 7),
                row_title_side = "right",
                row_split = path_class_data2_order$Pathway_TCGAclass2,
                row_title_gp = gpar(fontsize = 8),  
                row_title_ro = 0,  
                right_annotation = row_anno,  
                cluster_columns = F,  
                show_column_dend = F,
                show_column_names = F,
                column_split = temp_class_order$Class,
                column_title = c(str_c("C1 (n = ",class_num[1],")"),
                                 str_c("C2 (n = ",class_num[2],")"),
                                 str_c("C3 (n = ",class_num[3],")")),  
                column_title_side = "top",  
                column_title_gp = gpar(col = col_CLASS_COLOR),
                column_gap = unit(1, "mm"),
                top_annotation = col_anno  
)
print(heat)
dev.off()
table(path_class_data2_order$Pathway_TCGAclass2)

save.image(str_c(path,"15_TCGA_k3_metaPath_gsva_scores_complesheatmapC2.RData"))
