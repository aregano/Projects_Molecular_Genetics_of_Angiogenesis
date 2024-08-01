# Libraries
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


# rds

Liver <- readRDS("S:/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Liver_MacarenaFC/Figures_Paper/rds/Liver_20Jan21.Figures.rds")

# Color Palettes

my_palette_Rui_colors <- c("#E95A74", "#50B6EF", "#45FF8E", "#F4A753", "#A80519", "#880088", "#E28CF4", "#C1B80C","#FC0808",  "#0E47D8")

my_palette_Rui_colors_B <- c("#F4A753", "#E95A74", "#50B6EF", "#FC0808",  "#0E47D8", "#45FF8E", "#A80519", "#E28CF4", "#880088", "#C1B80C")

custom_rdylblu_palette <- c("#005AD9", "#0070FF", "#FFE135", "#FFDA13", "#E60026", "#C60012")

BYG <- c("#0000FF", "#3838C6", "#71718D", "#AAAA55", "#E2E21C", "#FFE200", "#FFAA00", "#FF7100", "#FF3800", "#FF0000")

grey_BYG <- c( "#D0C9D7", "#0000FF", "#3838C6", "#71718D", "#AAAA55", "#E2E21C", "#FFE200", "#FFAA00", "#FF7100", "#FF3800", "#FF0000")

rui_palette <- c( "#D0C9D7", "#71718D", "#AAAA55", "#E2E21C", "#FFE200", "#FFAA00", "#FF7100", "#FF3800", "#FF0000")

macarena_palette <- c("#FDDC85", "#6E2914")

Bestholtz_palette <- c("#DEDAD6", "#FEE392", "#FEC44E", "#FE9929", "#ED6F11", "#CC4C17", "#993411", "#65260C")

palette_Maca_Vln <- c("#606060", "#F94040", "#05BE78", "#55A0FB")

BuRd <- c("#052F60", "#2066AC", "#4392C3", "#92C5DE", "#D1E4EF", "#F7F7F7", "#FCDBC6", "#F4A482", "#D55F4D", "#B1172B", "#67001F")

Modified_palette <- c("#DEDAD6", "#FEE392", "#FEC44E", "#FE9929", "#ED6F11", "#CC4C17", "#9D3A21")

Modified_palette <- c("#DEDAD6", "#FEE392", "#FEC44E", "#FE9929", "#ED6F11", "#CC4C17", "#B24221")

###############################################################################################################################

## add the Arial font
font_add("Arial", regular = "arial.ttf",
         bold = "arialbd.ttf", italic = "ariali.ttf", bolditalic = "arialbi.ttf")

# Figure 4B. Violin Plot

genes <- list(c("DLL4", "NOTCH1", "RBPJ"))
genes <- as.data.frame(genes)
gene.names <- list(c("Dll4", "Notch1", "Rbpj"))
gene.names <- as.data.frame(gene.names)

# Testing 

VlnPlot(Liver, features = "CDH5", group.by = "Condition")+ theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.title.y=element_blank())+ggtitle("Hello")

p4 <- VlnPlot(Liver, features = "CDH5", group.by = "Condition", combine = F)

p4
#  Script

myplots <- vector('list', nrow(genes))


for (i in 1:nrow(genes)) {
  
  p21 <- VlnPlot(Liver, features = genes[i, 1], group.by = "Condition", cols = palette_Maca_Vln)+ NoLegend() + theme(text = element_text(family="Arial"), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.title.y=element_blank(), plot.title = element_text(face = "bold.italic"))+ggtitle(gene.names[i, 1])
  
  myplots[[i]] <- local({
    i <- i
    p21
  })
  
  # plegend <- VlnPlot(Liver, features = genes[i, 1], group.by = "Condition")
}



plegend <- VlnPlot(Liver, features = genes[i, 1], group.by = "Condition") +
  scale_fill_manual(values= palette_Maca_Vln, labels=c("Control", expression(italic("Dll4"^"iDEC")), expression(italic(Notch1)^"iDEC"), expression(italic(Rbpj)^"iDEC")))
plegend
legend <- get_legend(plegend + theme(legend.position = "bottom", legend.justification = "center", legend.direction = "horizontal", legend.text = element_text(size = 18)))

p25 <- patchwork::wrap_plots(myplots, ncol = 1)
p25
#myplots[1:nrow(genes)] + patchwork::plot_layout(byrow = T, widths = 25, heights = 5)
p26 <- plot_grid(p25, legend, ncol = 1, rel_heights = c(1, .1))

p26

x.axis <- ggdraw() + draw_label("Expression Level", fontface='bold', angle = 90, size = 24)

p27 <- plot_grid(x.axis, p26, ncol = 2, rel_widths = c(.05, 1))

p27

cairo_pdf("Plots/Figure_4/Figure_4B_VlnPlot.pdf", width = 9, height = 12, family = "Arial")
p27
dev.off()

jpeg("Plots/Figure_4/Figure_4B_VlnPlot.jpeg", width = 8, height = 12, units = 'in', res = 800)
p27
dev.off()

# p25 <- CombinePlots(plots = myplots, ncol = 1)
# p25

###############################################################################################################################

# Figure 4B. Dot Plot

p28 <- DotPlot(Liver, features = genes, group.by = "Condition", col.min = 0, dot.scale = 20, scale = F)+ theme(strip.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank()) +
  # scale_y_discrete(labels=c("Control", expression(italic("Dll4"^"iDEC")), expression(italic(Notch1)^"iDEC"), expression(italic(Rbpj)^"iDEC")))+
  # scale_x_discrete(labels=c(bquote(italic(Dll4)), bquote(italic(Notch1)), bquote(italic(Rbpj))))+
  scale_colour_gradientn(colours = c("#5757F9","white","#C60000"))

p28
setEPS()
postscript("Plots/Figure_4/Figure_4B_DotPlot.eps", width = 8, height = 5, family = "Arial", onefile = T, useKerning = F)
p28
dev.off()

jpeg("Plots/Figure_4/Figure_4B_DotPlot.jpeg", width = 8, height = 5, units = 'in', res = 800)
p28
dev.off()

# Trying the other palette

p28 <- DotPlot(Liver, features = genes, group.by = "Condition", col.min = 0, dot.scale = 20, scale = F)+ theme(strip.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank()) +
  scale_y_discrete(labels=c("Control", bquote(italic(Dll4)^iDEC), bquote(italic(Notch1)^iDEC), bquote(italic(Rbpj)^iDEC)))+
  scale_x_discrete(labels=c(bquote(italic(Dll4)), bquote(italic(Notch1)), bquote(italic(Rbpj))))+
  scale_colour_gradientn(colours = Bestholtz_palette)

cairo_pdf("Plots/Figure_4/Figure_4B_DotPlot_Palette.pdf", width = 8, height = 5, family = "Arial")
p28
dev.off()

jpeg("Plots/Figure_4/Figure_4B_DotPlot_Palette.jpeg", width = 8, height = 5, units = 'in', res = 800)
p28
dev.off()

###############################################################################################################################

# Figure 4C. UMAP

# Rename clusters and order them

Idents(Liver) <- Liver@meta.data$NewClustering

# Liver@active.ident <- Liver@meta.data$NewClustering

New_Cluster_Labels <- c("C1a - Arterial capillaries", "C1v - Venous capillaries", "C3 - Activated capillaries", "C0 - Unspecified quiescent capillaries",
                        "C4 - Endothelial tip cells", "C5ip - Proliferating S-phase", "C5p - G2-M", "C6 - X", "C2a - Large arteries", "C2v - Large veins", "C1")

# names(New_Cluster_Labels) <- levels(Liver)

# Liver <- RenameIdents(Liver, New_Cluster_Labels)

# levels(Liver) <- c("C0 - Unspecified quiescent capillaries", "C1a - Arterial capillaries", "C1v - Venous capillaries", "C2a - Large arteries", "C2v - Large veins",  "C3 - Activated capillaries", 
# "C4 - Endothelial tip cells", "C5ip - Proliferating S-phase", "C5p - G2-M", "C6 - X","C1")

DimPlot(Liver, label = T)

# Liver@active.ident -> Liver@meta.data$FigClustering

jpeg("Plots/Figure_4/Figure_4C_Axes.jpeg", width = 10, height = 6, units = 'in', res = 800)
DimPlot(Liver, cols = my_palette_Rui_colors_B, pt.size = 1.2)+theme(legend.text = element_text(size = 14))
dev.off()

p1 <- DimPlot(Liver, cols = my_palette_Rui_colors_B, pt.size = 1.2, split.by = "Condition", combine = T)


new_labels <- c("Control" = "Control", "Dll4KO" = "italic(Dll4)^iDEC", "NOTCH1KO" = "italic(Notch1)^iDEC", "RBPJKO" = "italic(Rbpj)^iDEC")
p1

jpeg("Plots/Figure_4/Figure_4C_Axes_Conditions.jpeg", width = 28, height = 6, units = 'in', res = 800)
p1+facet_grid(~Condition, labeller = as_labeller(new_labels, default = label_parsed))+theme(strip.text.x = element_text(size = 24), legend.text = element_text(size = 14))
dev.off()

p21 <- DimPlot(Liver, cols = my_palette_Rui_colors_B, pt.size = 1.2)+theme(legend.text = element_text(size = 14), axis.title.x = element_blank())+NoLegend()
p22 <- p1+facet_grid(~Condition, labeller = as_labeller(new_labels, default = label_parsed))+
  theme(strip.text.x = element_text(size = 24), legend.text = element_text(size = 14), axis.line.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks = element_blank())

p21
p22

p23 <- list(p21, p22)

design <- c(patchwork::area(1, 1, 1, 1), patchwork::area(1, 2, 1, 5.5))

p24 <- Reduce( `+`,  p23)+patchwork::plot_layout(design = design)


cairo_pdf("Plots/Figure_4/Figure_4C_Axes_Conditions_All.pdf",  width = 34, height = 6, family = "Arial")
p24
dev.off()


jpeg("Plots/Figure_4/Figure_4C_Axes_Conditions_All.jpeg", width = 34, height = 6, units = 'in', res = 800)
p24
dev.off()

###############################################################################################################################

# Figure 4D. Bar Plot

p1 <- dittoBarPlot(Liver, var = "FigClustering", group.by = "Condition", color.panel = my_palette_Rui_colors_B, xlab = NULL, main = NULL, y.breaks = c(0,0.25,0.5, 0.75, 1), min = 0, max = 1)+
  theme(axis.text=element_text(size=18), axis.title = element_text(size = 20), legend.text = element_text(size = 16), axis.text.x = element_text(angle = 25))+
  scale_x_discrete(labels=c("Control", bquote(italic(Dll4)^iDEC), bquote(italic(Notch1)^iDEC), bquote(italic(Rbpj)^iDEC)))+
  # scale_y_discrete(labels=c("0", "25", "50", "75", "100"))
  scale_y_continuous(labels = c("0", "25", "50", "75", "100"))
p1

cairo_pdf("Plots/Figure_4/Figure_4D_BarPlot.pdf",  width = 8, height = 6, family = "Arial")
p1
dev.off()


jpeg("Plots/Figure_4/Figure_4D_BarPlot.jpeg", width = 8, height = 6, units = 'in', res = 800)
p1
dev.off()

###############################################################################################################################

# NewFigure. Heatmap with Top10 Upregulated genes. I will keep going with the wilcox test for the DEG analysis, as it gives a better idea of the markers that give rise to the Tip cell cluster

DEG_FigClustering <- FindAllMarkers(Liver, logfc.threshold = 0.5, test.use = "wilcox", only.pos = T, verbose = T, min.diff.pct = -Inf)
DEG_FigClustering$order <- 1:nrow(DEG_FigClustering)
write.xlsx(DEG_FigClustering, "Tables/DEG_FigClustering.wilcox.xlsx", rowNames = T)
DEG_FigClustering.wilcox <- read_excel("Tables/DEG_FigClustering.wilcox.xlsx")


DEG_FigClustering.0.25diff <- FindAllMarkers(Liver, logfc.threshold = 0.5, test.use = "MAST", only.pos = T, verbose = T, min.diff.pct = 0.25)
DEG_FigClustering.0.25diff$order <- 1:nrow(DEG_FigClustering.0.25diff)
write.xlsx(DEG_FigClustering.0.25diff, "Tables/DEG_FigClustering.0.25diff.MAST.xlsx", rowNames = T)
DEG_FigClustering.0.25diff <- read_excel("Tables/DEG_FigClustering.0.25diff.MAST.xlsx")

DEG_FigClustering.0.1diff <- FindAllMarkers(Liver, logfc.threshold = 0.5, test.use = "MAST", only.pos = T, verbose = T, min.diff.pct = 0.1)
DEG_FigClustering.0.1diff$order <- 1:nrow(DEG_FigClustering.0.1diff)
write.xlsx(DEG_FigClustering.0.1diff, "Tables/DEG_FigClustering.0.1diff.MAST.xlsx", rowNames = T)
DEG_FigClustering.0.1diff <- read_excel("Tables/DEG_FigClustering.0.1diff.MAST.xlsx")

#  Doing DEG of only C4 v each cluster. GOing for wilcox

C4.markers <- FindMarkers(object = Liver, ident.1 = "C4 - Endothelial tip cells", ident.2 = c("C0 - Unspecified quiescent capillaries", "C1a - Arterial capillaries", "C1v - Venous capillaries", "C2a - Large arteries", "C2v - Large veins",  "C3 - Activated capillaries", "C5ip - Proliferating S-phase", "C5p - G2-M", "C6 - X"), min.pct = 0.25, min.diff.pct = 0.25, test.use = "wilcox")
C4.markers$order <- 1:nrow(C4.markers)
write.xlsx(C4.markers, "Tables/C4.markers.wilcox.xlsx", rowNames = T)
C4_markers_wilcox <- read_excel("Tables/C4.markers.wilcox.xlsx")

C4.markers.0.3diff <- FindMarkers(object = Liver, ident.1 = "C4 - Endothelial tip cells", ident.2 = c("C0 - Unspecified quiescent capillaries", "C1a - Arterial capillaries", "C1v - Venous capillaries", "C2a - Large arteries", "C2v - Large veins",  "C3 - Activated capillaries", "C5ip - Proliferating S-phase", "C5p - G2-M", "C6 - X"), min.pct = 0.25, min.diff.pct = 0.3)
C4.markers.0.3diff$order <- 1:nrow(C4.markers.0.3diff)
write.xlsx(C4.markers.0.3diff, "Tables/C4.markers.0.3diff.MAST.xlsx", rowNames = T)
C4_markers_0.3diff_MAST <- read_excel("Tables/C4.markers.0.3diff.MAST.xlsx")

DEG_FigClustering <- as.data.frame(DEG_FigClustering)
DEG_FigClustering$order <- 1:nrow(DEG_FigClustering)

Top10_logFC <- DEG_FigClustering %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
Top10 <- DEG_FigClustering %>% group_by(cluster) %>% top_n(n = 10, wt = -order)

write.xlsx(Top10, "Tables/Top10_per_cluster.wilcox.xlsx", rowNames = T)

# Once it is produced

Top10_per_cluster_wilcox <- read_excel("Tables/Top10_per_cluster.wilcox.xlsx")

Top10_per_cluster_wilcox -> Top10

Top10genes <- Top10$gene

Top10genes <- unique(Top10genes)

Idents(Liver) <- Liver@meta.data$FigClustering

avgexp <- AverageExpression(Liver, return.seurat = T)

avgexp@meta.data$FigClustering <- avgexp@active.ident

Top10genesR <- rev(Top10genes) 

Top10genesR

levels(avgexp) <- c("C0 - Unspecified quiescent capillaries", "C1a - Arterial capillaries", "C1v - Venous capillaries", "C2a - Large arteries", "C2v - Large veins",  "C3 - Activated capillaries", 
                    "C4 - Endothelial tip cells", "C5ip - Proliferating S-phase", "C5p - G2-M", "C6 - X","C1")

p2 <- DoHeatmap(avgexp, draw.lines = F, features = Top10genes, combine = T, label = F, group.bar = F,
                group.colors = my_palette_Rui_colors_B, size = 5, slot = "scale.data", angle = 60, group.bar.height = 0.01)+scale_fill_gradientn(colours = BuRd)+
  scale_y_discrete(label = Top10genesR)+theme(axis.text.x.bottom = element_text(size = 10, angle = 90, hjust = 1), axis.text.y = element_text(face = "italic"), legend.text = element_text(size = 12), legend.key.size = unit(0.5, "cm"))+guides(col = FALSE)


p2

cairo_pdf("Plots/Figure_4/Figure_4New_HeatMap.pdf",  width = 4, height = 14, family = "Arial")
p2
dev.off()


jpeg("Plots/Figure_4/Figure_4New_HeatMap.jpeg", width = 4, height = 14, units = 'in', res = 800)
p2
dev.off()



# p1 <- dittoHeatmap(avgexp, genes = Top10genes, cell.names.meta = "FigClustering",  cluster_rows = F, heatmap.colors = BuRd, highlight.features = "MALAT1", complex = T)
# annot.by = "FigClustering",

#  DotPlot Clusters

Top10genesR <- rev(Top10genes)

Top10genesDotPlot <- unique(Top10genes) %>% tolower() %>% str_to_title() %>% rev()

Top10genesDotPlotH <- unique(Top10genes) %>% tolower() %>% str_to_title()

p3 <- DotPlot(Liver, features = Top10genesR, group.by = "FigClustering", col.min = 0, dot.scale = 5, scale = T)+
  theme(strip.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_text(face = "italic"), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_x_discrete(label = Top10genesDotPlot)+
  scale_colour_gradientn(colours = Bestholtz_palette)+coord_flip()

p3

cairo_pdf("Plots/Figure_4/Figure_4New_DotPlot_Horizontal.pdf", width = 22, height = 5, family = "Arial")
DotPlot(Liver, features = Top10genes, group.by = "FigClustering", col.min = 0, dot.scale = 5, scale = T)+
  theme(strip.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic")) +
  scale_x_discrete(label = Top10genesDotPlotH)+
  scale_colour_gradientn(colours = Bestholtz_palette)
dev.off()


jpeg("Plots/Figure_4/Figure_4New_DotPlot_Horizontal.jpeg", width = 22, height = 5, units = 'in', res = 800)
DotPlot(Liver, features = Top10genes, group.by = "FigClustering", col.min = 0, dot.scale = 5, scale = T)+
  theme(strip.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic")) +
  scale_x_discrete(label = Top10genesDotPlotH)+
  scale_colour_gradientn(colours = Bestholtz_palette)
dev.off()


cairo_pdf("Plots/Figure_4/Figure_4New_DotPlot.pdf", width = 6, height = 22, family = "Arial")
p3
dev.off()


jpeg("Plots/Figure_4/Figure_4New_DotPlot.jpeg", width = 6, height = 22, units = 'in', res = 800)
p3
dev.off()

###############################################################################################################################

# FeaturePlots New Figure

genesFeaturePlots <- c("GJA5", "LTBP4", "MSR1", "RSPO3", "WNT2", "MKI67", "TOP2A", "STMN1", "RPL32", "RPL10A", "MALAT1", "MEIS1", "KCNE3", "ODC1", "CTNNB1", "RPL21")
Artery_Vein <- c("GJA5", "LTBP4", "MSR1", "RSPO3", "WNT2")
Artery_Vein <- as.data.frame(Artery_Vein)
Artery_Vein_text <- c("Gja5", "Ltbp4", "Msr1", "Rspo3", "Wnt2")

myplots <- vector('list', nrow(Artery_Vein))

for (i in 1:nrow(Artery_Vein)) {
  
  p3 <- FeaturePlot(Liver, features = Artery_Vein[i, 1], order = T, pt.size = 1.2, slot = "data", combine = F)
  
  p4 <- lapply(X = p3, FUN = function(p) p + scale_colour_gradientn(colors = Modified_palette) + ggtitle(bquote(~italic(.(Artery_Vein_text[i])))))
  
  p4 <- Reduce( `+`, p4 )+patchwork::plot_layout( ncol = 1 )
  
  myplots[[i]] <- local({
    i <- i
    p4
  })
  
  # plegend <- VlnPlot(Liver, features = genes[i, 1], group.by = "Condition")
}


myplots[5]

p25 <- patchwork::wrap_plots(myplots, ncol = 5)

p25


cairo_pdf("Plots/Figure_4/Figure_4NewFigure_FeaturePlot.Artery_Vein.pdf", width = 28, height = 6, family = "Arial")
p25
dev.off()

jpeg("Plots/Figure_4/Figure_4NewFigure_FeaturePlot.Artery_Vein.jpeg", width = 28, height = 6, units = 'in', res = 800)
p25
dev.off()

# Tip Cells + Proliferating

TipCell_Proliferating <- c("KCNE3", "ODC1", "MKI67", "TOP2A", "STMN1")
TipCell_Proliferating <- as.data.frame(TipCell_Proliferating)
TipCell_Proliferating_text <- c("Kcne3", "Odc1", "Mki67", "Top2a", "Stmn1")

myplots <- vector('list', nrow(TipCell_Proliferating))

for (i in 1:nrow(TipCell_Proliferating)) {
  
  p3 <- FeaturePlot(Liver, features = TipCell_Proliferating[i, 1], order = T, pt.size = 1.2, slot = "data", combine = F)
  
  p4 <- lapply(X = p3, FUN = function(p) p + scale_colour_gradientn(colors = Modified_palette) + ggtitle(bquote(~italic(.(TipCell_Proliferating_text[i])))))
  
  p4 <- Reduce( `+`, p4 )+patchwork::plot_layout( ncol = 1 )
  
  myplots[[i]] <- local({
    i <- i
    p4
  })
  
  # plegend <- VlnPlot(Liver, features = genes[i, 1], group.by = "Condition")
}


myplots[5]

p25 <- patchwork::wrap_plots(myplots, ncol = 5)

p25


cairo_pdf("Plots/Figure_4/Figure_4NewFigure_FeaturePlot.TipCell_Proliferating.pdf", width = 28, height = 6, family = "Arial")
p25
dev.off()

jpeg("Plots/Figure_4/Figure_4NewFigure_FeaturePlot.TipCell_Proliferating.jpeg", width = 28, height = 6, units = 'in', res = 800)
p25
dev.off()

# Quiescent+Activated

Quiescent_Activated <- c("MALAT1", "MEIS1", "RPL32", "RPL10A")
Quiescent_Activated <- as.data.frame(Quiescent_Activated)
Quiescent_Activated_text <- c("Malat1", "Meis1", "Rpl32", "Rpl10a")

myplots <- vector('list', nrow(Quiescent_Activated))

for (i in 1:nrow(Quiescent_Activated)) {
  
  p3 <- FeaturePlot(Liver, features = Quiescent_Activated[i, 1], order = T, pt.size = 1.2, slot = "data", combine = F)
  
  p4 <- lapply(X = p3, FUN = function(p) p + scale_colour_gradientn(colors = Modified_palette) + ggtitle(bquote(~italic(.(Quiescent_Activated_text[i])))))
  
  p4 <- Reduce( `+`, p4 )+patchwork::plot_layout( ncol = 1 )
  
  myplots[[i]] <- local({
    i <- i
    p4
  })
  
  # plegend <- VlnPlot(Liver, features = genes[i, 1], group.by = "Condition")
}


myplots[3]

p25 <- patchwork::wrap_plots(myplots, ncol = 5)

p25


cairo_pdf("Plots/Figure_4/Figure_4NewFigure_FeaturePlot.Quiescent_Activated.pdf", width = 28, height = 6, family = "Arial")
p25
dev.off()

jpeg("Plots/Figure_4/Figure_4NewFigure_FeaturePlot.Quiescent_Activated.jpeg", width = 28, height = 6, units = 'in', res = 800)
p25
dev.off()

# C6 - X
C6 <- c("CTNNB1", "RPL21")
C6 <- as.data.frame(C6)
C6_text <- c("Ctnnb1", "Rpl21")

myplots <- vector('list', nrow(C6))

for (i in 1:nrow(C6)) {
  
  p3 <- FeaturePlot(Liver, features = C6[i, 1], order = T, pt.size = 1.2, slot = "data", combine = F)
  
  p4 <- lapply(X = p3, FUN = function(p) p + scale_colour_gradientn(colors = Modified_palette) + ggtitle(bquote(~italic(.(C6_text[i])))))
  
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


cairo_pdf("Plots/Figure_4/Figure_4NewFigure_FeaturePlot.C6.pdf", width = 28, height = 6, family = "Arial")
p25
dev.off()

jpeg("Plots/Figure_4/Figure_4NewFigure_FeaturePlot.C6.jpeg", width = 28, height = 6, units = 'in', res = 800)
p25
dev.off()

###############################################################################################################################

# Figure 4E Tip Cell Cluster Top 20 or Top 50 genes

DEG_FigClustering.wilcox <- read_excel("Tables/DEG_FigClustering.wilcox.xlsx")

# DEG_FigClustering_MAST$order <- 1:nrow(DEG_FigClustering_MAST)

Top50_C4_Tip_Cells <- filter(DEG_FigClustering.wilcox, cluster == "C4 - Endothelial tip cells") %>% top_n(n = 50, wt = -order)

write.xlsx(Top50_C4_Tip_Cells, "Tables/Top50_C4_Tip_Cells.wilcox.xlsx", rowNames = T)

# Figure 4E DotPlot. Top 50

Top50_C4_Tip_Cells_genes <- Top50_C4_Tip_Cells$gene %>% rev()

Top50genesC4DotPlot <- unique(Top50_C4_Tip_Cells_genes) %>% tolower() %>% str_to_title()

p3 <- DotPlot(Liver, features = Top50_C4_Tip_Cells_genes, group.by = "FigClustering", col.min = 0, dot.scale = 5, scale = T)+ 
  theme(strip.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_text(face = "italic"), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_x_discrete(label = Top50genesC4DotPlot)+
  scale_colour_gradientn(colours = Bestholtz_palette)+coord_flip()

p3

cairo_pdf("Plots/Figure_4/Figure_4E_C4_tip_cells_DotPlot_Top50.pdf", width = 5, height = 14, family = "Arial")
p3
dev.off()

jpeg("Plots/Figure_4/Figure_4E_C4_tip_cells_DotPlot_Top50.jpeg", width = 5, height = 14, units = 'in', res = 800)
p3
dev.off()



# Figure 4E DotPlot. Top 20

Top20_C4_Tip_Cells_genes <- Top50_C4_Tip_Cells_genes[34:50] %>% rev()

Top20_C4_Tip_Cells_genes <- append(Top20_C4_Tip_Cells_genes, c("KCNE3", "ODC1", "PEG10", "ESM1", "APLN", "VEGFA", "TPI1" , "MIF", "ENO1"), after = 0) %>% unique()

Top20_C4_Tip_Cells_genes

Top20genesC4DotPlot <- unique(Top20_C4_Tip_Cells_genes) %>% tolower() %>% str_to_title() %>% rev()

p3 <- DotPlot(Liver, features = Top20_C4_Tip_Cells_genes, group.by = "FigClustering", col.min = 0, dot.scale = 5, scale = T)+
  theme(strip.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_text(face = "italic"), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  scale_x_discrete(label = Top20genesC4DotPlot)+
  scale_colour_gradientn(colours = Bestholtz_palette)+coord_flip()

p3

cairo_pdf("Plots/Figure_4/Figure_4E_C4_tip_cells_DotPlot_Top20_Modified_palette.pdf", width = 5, height = 7, family = "Arial")
DotPlot(Liver, features = Top20_C4_Tip_Cells_genes, group.by = "FigClustering", col.min = 0, dot.scale = 5, scale = T)+
  theme(strip.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_text(face = "italic"), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  scale_x_discrete(label = Top20genesC4DotPlot)+
  scale_colour_gradientn(colours = Modified_palette)+coord_flip()
dev.off()

cairo_pdf("Plots/Figure_4/Figure_4E_C4_tip_cells_DotPlot_Top20_Gray_Red.pdf", width = 5, height = 7, family = "Arial")
DotPlot(Liver, features = Top20_C4_Tip_Cells_genes, group.by = "FigClustering", col.min = 0, dot.scale = 5, scale = T)+
  theme(strip.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_text(face = "italic"), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  scale_x_discrete(label = Top20genesC4DotPlot)+
  scale_colour_gradientn(colours = c("gray", "#C60000"))+coord_flip()
dev.off()

cairo_pdf("Plots/Figure_4/Figure_4E_C4_tip_cells_DotPlot_Top20_Horizontal.pdf", width = 9, height = 4, family = "Arial")
DotPlot(Liver, features = Top20_C4_Tip_Cells_genes, group.by = "FigClustering", col.min = 0, dot.scale = 5, scale = T)+
  theme(strip.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x.bottom = element_text(face = "italic"), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  scale_x_discrete(label = Top20genesC4DotPlot)+
  scale_colour_gradientn(colours = Bestholtz_palette)
dev.off()

cairo_pdf("Plots/Figure_4/Figure_4E_C4_tip_cells_DotPlot_Top20.pdf", width = 5, height = 7, family = "Arial")
p3
dev.off()


jpeg("Plots/Figure_4/Figure_4E_C4_tip_cells_DotPlot_Top20.jpeg", width = 5, height = 7, units = 'in', res = 800)
p3
dev.off()

# Figure 4E Heatmap Top 50 C4 Endothelial Tip Cell

Idents(Liver) <- Liver@meta.data$FigClustering

avgexp <- AverageExpression(Liver, return.seurat = T)

avgexp@meta.data$FigClustering <- avgexp@active.ident

Top50_C4_Tip_Cells_genes <- Top50_C4_Tip_Cells$gene

Top50_C4_Tip_Cells_genes_Heatmap <- rev(Top50_C4_Tip_Cells_genes) %>% tolower() %>% str_to_title()

Top50_C4_Tip_Cells_genes_Heatmap

levels(avgexp) <- c("C0 - Unspecified quiescent capillaries", "C1a - Arterial capillaries", "C1v - Venous capillaries", "C2a - Large arteries", "C2v - Large veins",  "C3 - Activated capillaries", 
                    "C4 - Endothelial tip cells", "C5ip - Proliferating S-phase", "C5p - G2-M", "C6 - X","C1")

p2 <- DoHeatmap(avgexp, draw.lines = F, features = Top50_C4_Tip_Cells_genes, combine = T, label = F, group.bar = F,
                size = 5, slot = "scale.data", angle = 60, group.bar.height = 0.01)+scale_fill_gradientn(colours = BuRd)+
  scale_y_discrete(label = Top50_C4_Tip_Cells_genes_Heatmap)+
  theme(axis.text.x.bottom = element_text(size = 10, angle = 90, hjust = 1), axis.text.y = element_text(face = "italic"), legend.text = element_text(size = 12), legend.key.size = unit(0.5, "cm"))+guides(col = FALSE)  


p2

cairo_pdf("Plots/Figure_4/Figure_4E_HeatMap_C4_tip_cell_Top50.pdf", width = 4, height = 13, family = "Arial")
p2
dev.off()

jpeg("Plots/Figure_4/Figure_4E_HeatMap_C4_tip_cell_Top50.jpeg", width = 4, height = 10, units = 'in', res = 800)
p2
dev.off()

# Figure 4E Heatmap Top 20 C4 Endothelial Tip Cell

Top20_C4_Tip_Cells_genes <- Top50_C4_Tip_Cells_genes[34:50] %>% rev()

Top20_C4_Tip_Cells_genes <- append(Top20_C4_Tip_Cells_genes, c("KCNE3", "ODC1", "PEG10", "ESM1", "APLN", "VEGFA", "TPI1" , "MIF", "ENO1"), after = 0) %>% unique()

Top20_C4_Tip_Cells_genes

Top20_C4_Tip_Cells_genes_Heatmap <- rev(Top20_C4_Tip_Cells_genes) %>% tolower() %>% str_to_title()

Top20_C4_Tip_Cells_genes_Heatmap

p2 <- DoHeatmap(avgexp, draw.lines = F, features = Top20_C4_Tip_Cells_genes, combine = T, label = F, group.bar = F,
                size = 5, slot = "scale.data", angle = 60, group.bar.height = 0.01)+scale_fill_gradientn(colours = BuRd)+
  scale_y_discrete(label = Top20_C4_Tip_Cells_genes_Heatmap)+
  theme(axis.text.x.bottom = element_text(size = 10, angle = 90, hjust = 1), axis.text.y = element_text(face = "italic"), legend.text = element_text(size = 12), legend.key.size = unit(0.5, "cm"))+guides(col = FALSE)  


p2

cairo_pdf("Plots/Figure_4/Figure_4E_HeatMap_C4_tip_cell_Top20.pdf", width = 4, height = 6, family = "Arial")
p2
dev.off()


jpeg("Plots/Figure_4/Figure_4E_HeatMap_C4_tip_cell_Top20.jpeg", width = 4, height = 6, units = 'in', res = 800)
p2
dev.off()

#  Bar graph top20 or top50 genes C4 Tip cells v All

DEG_C4_v_all <- FindMarkers(Liver, ident.1 = "C4 - Endothelial tip cells", logfc.threshold = 0.5, test.use = "wilcox", only.pos = T, verbose = T)

write.xlsx(DEG_C4_v_all, "Tables/DEG_C4vAll.wilcox.xlsx", rowNames = T)

DEG_C4vAll_wilcox <- read_excel("Tables/DEG_C4vAll.wilcox.xlsx")

# Bar Graph Top 20

Top20_DEG_C4vAll <- head(DEG_C4vAll_wilcox, 20)

Top20_DEG_C4vAll$order <- 1:nrow(Top20_DEG_C4vAll)

Top20_DEG_C4vAll$gene <- tolower(Top20_DEG_C4vAll$gene) %>% str_to_title()

Top20_DEG_C4vAll$gene = fct_rev(fct_inorder(Top20_DEG_C4vAll$gene))

Top20_DEG_C4vAll$avg_log2FC = round(Top20_DEG_C4vAll$avg_log2FC, 2)

p5 <- Top20_DEG_C4vAll %>%
  group_by(order) %>%
  ggplot(., aes(x = avg_log2FC, y = gene)) + 
  geom_col(fill=c("#C60000"), colour="black") + 
  geom_text(aes(label=avg_log2FC),  hjust=-0.2, fontface = "bold", size = 3) +
  xlab(expression(paste(Log[2], "FC"))) +
  scale_x_continuous(limits = c(0, 3.2), expand = c(0,0))+
  theme_classic()+
  theme(axis.text.y = element_text(face = "italic"), axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5))+
  ggtitle("C4 - Endothelial Tip Cell vs. All")

p5

jpeg("Plots/Figure_4/Figure_4E_Bar_Graph_C4_tip_cells_Top20.jpeg", width = 7, height = 5, units = 'in', res = 800)
p5
dev.off()

# Now with the Top 50 genes

# DEG_C4vAll_wilcox <- append(DEG_C4vAll_wilcox, c("KCNE3", "ODC1", "PEG10", "ESM1", "APLN", "VEGFA", "TPI1" , "MIF", "ENO1"), after = 0) %>% unique()

Top50_DEG_C4vAll <- head(DEG_C4vAll_wilcox, 50)

Top50_DEG_C4vAll$...1 <- tolower(Top50_DEG_C4vAll$...1) %>% str_to_title()

Top50_DEG_C4vAll$...1 = fct_rev(fct_inorder(Top50_DEG_C4vAll$...1))

# Top50_DEG_C4vAll$avg_log2FC = fct_rev(fct_inorder(Top50_DEG_C4vAll$avg_log2FC))

Top50_DEG_C4vAll <- Top50_DEG_C4vAll[order(Top50_DEG_C4vAll$avg_log2FC, decreasing = T),]

Top50_DEG_C4vAll$avg_log2FC = round(Top50_DEG_C4vAll$avg_log2FC, 2)

Top50_DEG_C4vAll$order <- 1:nrow(Top50_DEG_C4vAll)

Top50_DEG_C4vAll$...1 <- factor(Top50_DEG_C4vAll$...1, levels = Top50_DEG_C4vAll$...1[order(Top50_DEG_C4vAll$avg_log2FC)])

Top50_DEG_C4vAll

p5 <- Top50_DEG_C4vAll %>%
  group_by(order) %>%
  ggplot(., aes(x = avg_log2FC, y = ...1)) + 
  geom_col(fill=c("#C60000"), colour="black") + 
  geom_text(aes(label=avg_log2FC),  hjust=-0.2, fontface = "bold", size = 3) +
  xlab(expression(paste(Log[2], "FC"))) +
  scale_x_continuous(limits = c(0, 3.2), expand = c(0,0))+
  theme_classic()+
  theme(axis.text.y = element_text(face = "italic"), axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5))+
  ggtitle("C4 - Endothelial Tip Cell vs. All")

p5

cairo_pdf("Plots/Figure_4/Figure_4E_Bar_Graph_C4_tip_cells_Top50.pdf", width = 7, height = 11, family = "Arial")
p5
dev.off()

jpeg("Plots/Figure_4/Figure_4E_Bar_Graph_C4_tip_cells_Top50.jpeg", width = 7, height = 11, units = 'in', res = 800)
p5
dev.off()


saveRDS(Liver, "./rds/Liver_20Jan21.Figures.rds")
