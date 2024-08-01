library(Seurat)
library(dittoSeq)

p1 <- dittoBarPlot(Liver, group.by = "Condition", var = "FigClustering", color.panel = my_palette_Rui_colors_B, x.reorder = c(1,3,4,5,6,9,10,11,12,13,14,7,2,8))+ggtitle("AllConditions_Rui_ppt")

p2 <- dittoBarPlot(Liver.query, group.by = "Condition", var = "FigClustering", color.panel = my_palette_Rui_colors_B, x.reorder = c(1,2,3,4,5,7,8,9,10,11,12,6))+ggtitle("Standard_AddingDll4_KO_4days")

p3 <- dittoBarPlot(Liver1, group.by = "Condition", var = "FigClustering", color.panel = my_palette_Rui_colors_B)+ggtitle("Standard")

p1/p2/p3

cairo_pdf("Plots/Comparison_Dll4KO4days_Barplot.pdf",  width = 7, height = 18, family = "Arial")
p1/p2/p3
dev.off()
