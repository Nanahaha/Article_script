## ------ drawing heatmap for SubMap of bevacizumab and subtype------ ##

library(stringr)
library(ComplexHeatmap)
library(circlize)  
library(RColorBrewer)  
library(pheatmap)

PATH <- getwd()
path <- str_c(PATH,"bevacizumab_exp_submap\\")

SAmatrix <- read.table(file = str_c(path,"SubMap_nominal_p_matrix_Fisher.txt"),
                       header = T, fill = T, check.names = F)
SAmatrix
rownames(SAmatrix) <- SAmatrix$Name

SAmatrix2 <- as.matrix(SAmatrix[,-1])  
rownames(SAmatrix2) <- c("C1","C2","C3")
colnames(SAmatrix2) <- c("Suboptimal",
                         "Optimal")

COLOR <- colorRamp2(c(0,  0.4, 0.6, 1),brewer.pal(9,"RdBu")[2:5])


dev.size()
pdf(file = str_c(path,"TCGA_beva_k3_submap_SAmatrix_heatmap.pdf"),
    width = 3, height = 3)
heat <- Heatmap(SAmatrix2,
                col = COLOR,
                rect_gp = gpar(col = "white", lwd = 2),  
                heatmap_legend_param = list(title = "nominal p-value",
                                            legend_height = unit(4, "cm"),
                                            title_position = "leftcenter-rot"),
                cluster_rows = F,
                show_row_dend = F,
                row_names_side = "left",
                row_names_gp = gpar(fontsize = 10),
                row_title = "TCGA Subtypes",
                row_title_gp = gpar(fontsize = 10),
                cluster_columns = F,
                show_column_dend = F,
                column_names_rot = 0,  
                column_names_gp = gpar(fontsize = 10),
                column_names_centered = T,
                column_title = "Bevacizumab Response",
                column_title_side = "bottom",
                column_title_gp = gpar(fontsize = 10),
                cell_fun = function(j, i, x, y, width, height, fill){
                  grid.text(sprintf("%0.2f", SAmatrix2[i,j]), x, y,
                            gp = gpar(fontsize = 10))
                })
print(heat)
dev.off()

save.image(file = str_c(PATH,"30_Submap_heatmap_beva.RData"))
