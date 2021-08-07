## ------ heatmap for SubMap result using ComplexHeatmap function ------ ##

library(stringr)
library(ComplexHeatmap)
library(circlize)  
library(RColorBrewer) 
library(pheatmap)

path <- getwd()

SAmatrix <- read.table(file = str_c(path,"SubMap_nominal_p_matrix_Fisher.txt"),
                       header = T, fill = T, check.names = F)
rownames(SAmatrix) <- SAmatrix$Name
SAmatrix2 <- SAmatrix[,c("Name","B2","B1","B3")]
colnames(SAmatrix2)[2:4] <- c("C1","C2","C3")
SAmatrix2
SAmatrix_t <- t(SAmatrix2[,-1]) 
colnames(SAmatrix_t) <- c("C1","C2","C3")

# 
display.brewer.pal(9,"Spectral")
COLOR <- colorRamp2(c(0,  0.4, 0.6, 1),
                    rev(brewer.pal(9,"Spectral")[6:9])) 

load(file = str_c(path,"5_Submap_heatmap.RData"))

pdf(file = str_c(path,"TCGA_ICGC_k3_submap_SAmatrix_heatmap.pdf"),
    width = 3.7, height = 3)
heat <- Heatmap(SAmatrix_t,
                col = COLOR,
                rect_gp = gpar(col = "white", lwd = 2),  
                heatmap_legend_param = list(title = "Nominal p-value",
                                            legend_height = unit(4, "cm"),
                                            title_position = "leftcenter-rot"),
                cluster_rows = F,
                show_row_dend = F,
                row_names_side = "left",
                row_names_gp = gpar(fontsize = 13),
                row_title = "ICGC",
                row_title_gp = gpar(fontsize = 13),
                cluster_columns = F,
                show_column_dend = F,
                column_names_rot = 0, 
                column_names_gp = gpar(fontsize = 13),
                column_title = "TCGA",
                column_title_side = "bottom",
                column_title_gp = gpar(fontsize = 13),
                #bottom_annotation = col_anno,  
                cell_fun = function(j, i, x, y, width, height, fill){
                  grid.text(sprintf("%0.2f", SAmatrix_t[i,j]), x, y,
                            gp = gpar(fontsize = 10))
                })
print(heat)
dev.off()

save.image(file = str_c(path,"8_Submap_heatmap.RData"))
