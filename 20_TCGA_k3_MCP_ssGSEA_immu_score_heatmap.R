## ------ drawing heatmap for immune cell infiltration ------ ##

library(pheatmap)
library(gplots)  # 函数heatmap.2
library(stringr)
library(ggplot2)
library(RColorBrewer)  # 配色包

path <- getwd()

mcp <- read.table(file = str_c(path,"k3\\TCGA_k3_MCPcounter_immu_scoreslog_class.txt"),
                  sep = "\t", header = T, fill = T, check.names = F)
mcp[1:4,1:4]
ssGSEA <- read.table(file = str_c(path,"k3\\TCGA_k3_ssGSEA_immu_scores_class.txt"),
                     sep = "\t", header = T, fill = T, check.names = F)
ssGSEA[1:4,1:4]

mcp_order <- mcp[order(mcp$Class),]
mcp_order$Class
mcp_order$Class <- factor(mcp_order$Class,
                          levels = c("C1","C2","C3"),
                          ordered = T)
rownames(mcp_order) <- mcp_order$Patients
mcp_order[1:4,1:4]
mcp_order_t <- t(mcp_order[,-c(1:2)])
mcp_order_t[1:4,1:4]


## MCPcounter
#Cancer associated fibroblast，Endothelial cell和Myeloid dendritic cell，"T cell",
cells <- c("T cell","Cancer associated fibroblast","Endothelial cell","Myeloid dendritic cell")

annotation_col = data.frame(Class = mcp_order$Class)
rownames(annotation_col) = mcp_order$Patients

col_CLASS_COLOR <- brewer.pal(3,"Set2")

ann_colors = list(
  Class = c(C1 = col_CLASS_COLOR[1],
            C2 = col_CLASS_COLOR[2],
            C3 = col_CLASS_COLOR[3])
)

pdf(file = str_c(path,"k3\\TCGA_k3_MCP_immu_score_heatmap.pdf"),
    width = 8, height = 1)
pheatmap(mcp_order_t[which(rownames(mcp_order_t) %in% cells),],  
         scale = "row",
         cluster_cols = F,
         cluster_rows = T,
         treeheight_row = 8,
         cellheight = 10,
         color = colorRampPalette(rev(brewer.pal(11,"RdBu")))(100),
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         show_colnames = F,
         fontsize_row = 7
)
dev.off()

## ------ ssGSEA ------ ##
ssGSEA[1:4,1:4]
rownames(ssGSEA) <- ssGSEA$Patients
ssGSEA_t <- t(ssGSEA[,-c(1:2)])
ssGSEA_t[1:4,1:4]
ssGSEA_t_order <- ssGSEA_t[,colnames(mcp_order_t)]
mcp_order_t[1:4,1:4]
ssGSEA_t_order[1:4,1:4]


pdf(file = str_c(path,"k3\\TCGA_k3_ssGSEA_immu_score_heatmap.pdf"),
    width = 8, height = 5)
pheatmap(ssGSEA_t_order,  
         scale = "row", 
         cluster_cols = F,
         cluster_rows = T,
         treeheight_row = 8,
         cellheight = 10,
         color = colorRampPalette(rev(brewer.pal(11,"RdBu")))(100), 
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         show_colnames = F,
         fontsize_row = 7
)
dev.off()

save.image(str_c(path,"20_TCGA_k3_MCP_ssGSEA_immu_score_heatmap.RData"))
