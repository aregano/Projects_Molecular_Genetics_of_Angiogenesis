library(Seurat)
library(ggplot2)
library(openxlsx)

table(Liver@meta.data$FigClustering)


C3vC4DEG <- FindMarkers(Liver, ident.1 = "C3 - Activated capillaries", ident.2 = "C4 - Endothelial tip cells", only.pos = T, min.pct = 0.5, min.diff.pct = 0.1)

write.xlsx(C3vC4DEG, "./Tables/C3vC4DEG.xlsx", rowNames = T)

C3vAllDEG <- FindMarkers(Liver, ident.1 = "C3 - Activated capillaries", only.pos = T, min.pct = 0.5, min.diff.pct = 0.2)

write.xlsx(C3vAllDEG, "./Tables/C3vAllDEG.xlsx", rowNames = T)
