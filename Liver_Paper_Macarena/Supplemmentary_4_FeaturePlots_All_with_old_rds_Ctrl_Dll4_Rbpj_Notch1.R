library(Seurat)
library(ggplot2)
library(patchwork)
library("dittoSeq")
library("magick")
library("RColorBrewer")
library(cowplot)
library(dplyr)
library(extrafont)
loadfonts(device = "win")
library(openxlsx)
library(ComplexHeatmap)
library("stringr")
library("formattable")
library(readxl)
library(tidyverse)
library(showtext)
library(fgsea)

Liver <- readRDS("rds/Groups/Liver_20Jan21.Figures.rds")

table(Liver@meta.data$Condition)
p1 <- DimPlot(Liver, group.by = "Condition")

p1

Bestholtz_palette <- c("#DEDAD6", "#FEE392", "#FEC44E", "#FE9929", "#ED6F11", "#CC4C17", "#993411", "#65260C")

# C3 - Activated Capillaries

table(Liver@meta.data$FigClustering)

DEG <- FindMarkers(Liver, group.by = "FigClustering", ident.1 = "C3 - Activated capillaries", min.diff.pct = 0.22, only.pos = T, logfc.threshold = 0.25)

DEG <- FindMarkers(Liver, group.by = "FigClustering", ident.1 = "C3 - Activated capillaries", ident.2 = "C4 - Endothelial tip cells", min.diff.pct = 0.3, only.pos = T, logfc.threshold = 0.25)

DEG.names.text <- c("Gja5","Ltbp4", "Msr1", "Rspo3", "Wnt2", "Meis1", "Rpl32", "Rpl10a", "Rab3b", "Kcne3", "Odc1", "Mki67", "Top2a", "Stmn1", "Mcm2", "Mcm3", "Ctnnb1", "Cd34")

DEG.names <- rownames(DEG)

DEG.names <- append(DEG.names, c("CXCL10"))

DEG.names.text <- str_to_title(DEG.names)

myplots <- vector('list', length(DEG.names))

DEG.names <- toupper(DEG.names.text)

for (i in 1:length(DEG.names)) {
  
  p3 <- FeaturePlot(Liver, features = DEG.names[i], order = T, pt.size = 1, slot = "data", combine = F)
  
  p4 <- lapply(X = p3, FUN = function(p) p + scale_colour_gradientn(colors = Bestholtz_palette) + ggtitle(bquote(~italic(.(DEG.names.text[i])))))
  
  p4 <- Reduce( `+`, p4 )+patchwork::plot_layout( ncol = 1 )
  
  myplots[[i]] <- local({
    i <- i
    p4
  })
  
  # plegend <- VlnPlot(Liver, features = genes[i, 1], group.by = "Condition")
}

myplots[10]

p25 <- patchwork::wrap_plots(myplots, ncol = 4)

p25

cairo_pdf("Plots/Figure_4/Suplementary/FigureS4_FeaturePlots.pdf", width = 20, height = 25, family = "Arial")
p25
dev.off()

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


# FeaturePlot Cxcl10

FeaturePlot(Liver, features = "CXCL10", pt.size = 1.2, cols = Bestholtz_palette, order = T)

####################################

# CD34

p3 <- FeaturePlot(Liver, features = "CD34", order = T, pt.size = 1, slot = "data", combine = T)

p4 <- lapply(X = p3, FUN = function(p) p + scale_colour_gradientn(colors = Bestholtz_palette) + ggtitle(bquote(~italic(.("Cd34")))))


cairo_pdf("Plots/Figure_4/Suplementary/Figure_S4.1_FeaturePlot.CD34.pdf", width = 6, height = 6, family = "Arial")
p4
dev.off()

jpeg("Plots/Figure_4/Suplementary/Figure_S4.1_FeaturePlot.CD34.jpeg", width = 6, height = 6, units = 'in', res = 800)
p4
dev.off()

#################################

# Find out which one is the good rds 

Liver.Mapped <- readRDS("S:/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Liver_MacarenaFC/Figures_Paper/rds/Liver_20Jan21.Figures.All.Conditions.OldDll4Myc.rds")
Liver.Mapped <- subset(Liver.Mapped, subset = Condition == "Control" | Condition == "Dll4KO" | Condition == "NOTCH1KO" | Condition == "RBPJKO")

levels(Liver.Mapped) <- c("Control", "Dll4KO", "NOTCH1KO", "RBPJKO")

table(Liver.Mapped@metadata$Condition)

p2 <- DimPlot(Liver.Mapped)

Liver.Mapped.All <- readRDS("S:/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Liver_MacarenaFC/Figures_Paper/rds/Liver_20Jan21.Figures.All.Conditions.OldDll4Myc+Dll4KO4d.rds")
Liver.Mapped.All <- subset(Liver.Mapped.All, subset = Condition == "Control" | Condition == "Dll4KO" | Condition == "NOTCH1KO" | Condition == "RBPJKO")

levels(Liver.Mapped.All) <- c("Control", "Dll4KO", "NOTCH1KO", "RBPJKO")

table(Liver.Mapped.All@meta.data$Condition)

p3 <- DimPlot(Liver.Mapped.All)

p1+p2+p3
