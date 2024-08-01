library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(cowplot)
library(stringr)

Bestholtz_palette <- c("#DEDAD6", "#FEE392", "#FEC44E", "#FE9929", "#ED6F11", "#CC4C17", "#993411", "#65260C")

Liver <- readRDS("S:/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Liver_MacarenaFC/Figures_Paper/rds/Liver_20Jan21.Figures.All.Conditions.OldDll4Myc.rds")
Liver <- subset(Liver, subset = Condition == "Control" | Condition == "Dll4KO" | Condition == "NOTCH1KO" | Condition == "RBPJKO")



#  C1 QUiescent

DEG <- FindMarkers(Liver, group.by = "FigClustering", ident.1 = "C0 - Unspecified quiescent capillaries", min.diff.pct = 0.23, only.pos = T, logfc.threshold = 0.25)

DEG <- FindMarkers(Liver, group.by = "FigClustering", ident.1 = "C0 - Unspecified quiescent capillaries", min.pct = 0.9, only.pos = T, logfc.threshold = 0)

DEG <- FindMarkers(Liver, group.by = "FigClustering", ident.1 = "C0 - Unspecified quiescent capillaries", min.diff.pct = 0.6, only.pos = F, logfc.threshold = 0)

DEG.names <- rownames(DEG)

DEG.names <- DEG.names[1:9]

p1 <- FeaturePlot(Liver, features = DEG.names, cols = Bestholtz_palette)

p1


# C3 - Activated Capillaries

DEG <- FindMarkers(Liver, group.by = "FigClustering", ident.1 = "C3 - Activated capillaries", min.diff.pct = 0.22, only.pos = T, logfc.threshold = 0.25)

DEG.names <- rownames(DEG)

DEG.names <- append(DEG.names, c("MALAT1", "RPL36AL"))

myplots <- vector('list', 16)

DEG.names.text <- str_to_title(DEG.names)

for (i in 1:16) {
  
  p3 <- FeaturePlot(Liver, features = DEG.names[i], order = T, pt.size = 1, slot = "data", combine = F)
  
  p4 <- lapply(X = p3, FUN = function(p) p + scale_colour_gradientn(colors = Bestholtz_palette) + ggtitle(bquote(~italic(.(DEG.names.text[i])))))
  
  p4 <- Reduce( `+`, p4 )+patchwork::plot_layout( ncol = 1 )
  
  myplots[[i]] <- local({
    i <- i
    p4
  })
  
  # plegend <- VlnPlot(Liver, features = genes[i, 1], group.by = "Condition")
}

myplots[15]

p25 <- patchwork::wrap_plots(myplots, ncol = 4)

p25

cairo_pdf("Plots/Figure_4/Suplementary/Figure_FeaturePlots_C3_and_C0.pdf", width = 20, height = 20, family = "Arial")
p25
dev.off()

FeaturePlot(Liver, features = c("RPS27"), cols = Bestholtz_palette)

# Only for Malat1 and Eif3f

DEG.names <- c("MALAT1", "EIF3F")

myplots <- vector('list', 2)

DEG.names.text <- str_to_title(DEG.names)

for (i in 1:2) {
  
  p3 <- FeaturePlot(Liver, features = DEG.names[i], order = F, pt.size = 1, slot = "data", combine = F)
  
  p4 <- lapply(X = p3, FUN = function(p) p + scale_colour_gradientn(colors = Bestholtz_palette) + ggtitle(bquote(~italic(.(DEG.names.text[i])))))
  
  p4 <- Reduce( `+`, p4 )+patchwork::plot_layout( ncol = 1 )
  
  myplots[[i]] <- local({
    i <- i
    p4
  })
  
  # plegend <- VlnPlot(Liver, features = genes[i, 1], group.by = "Condition")
}

myplots[2]

p25 <- patchwork::wrap_plots(myplots, ncol = 2)

p25

cairo_pdf("Plots/Figure_4/Suplementary/Figure_FeaturePlots_Eif3f_Malat1.pdf", width = 10, height = 5, family = "Arial")
p25
dev.off()



