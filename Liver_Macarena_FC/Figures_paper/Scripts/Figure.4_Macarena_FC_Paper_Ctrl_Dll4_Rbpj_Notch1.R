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
library(fgsea)


# rds

Liver <- readRDS("S:/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Liver_MacarenaFC/Figures_Paper/rds/Liver_20Jan21.Figures.All.Conditions.rds")
Liver <- subset(Liver, subset = Condition == "Control" | Condition == "Dll4KO" | Condition == "NOTCH1KO" | Condition == "RBPJKO")

levels(Liver) <- c("Control", "Dll4KO", "NOTCH1KO", "RBPJKO")

Liver@active.ident -> Liver@meta.data$Condition

# Color Palettes

my_palette_Rui_colors <- c("#E95A74", "#50B6EF", "#45FF8E", "#F4A753", "#A80519", "#880088", "#E28CF4", "#C1B80C","#FC0808",  "#0E47D8")

my_palette_Rui_colors_B <- c("#F4A753", "#E95A74", "#50B6EF", "#FC0808",  "#0E47D8", "#45FF8E", "#A80519", "#E28CF4", "#880088", "#C1B80C")

custom_rdylblu_palette <- c("#005AD9", "#0070FF", "#FFE135", "#FFDA13", "#E60026", "#C60012")

macarena_palette <- c("#FDDC85", "#6E2914")

Bestholtz_palette <- c("#DEDAD6", "#FEE392", "#FEC44E", "#FE9929", "#ED6F11", "#CC4C17", "#993411", "#65260C")

palette_Maca_Vln <- c("#606060", "#F94040", "#05BE78", "#55A0FB")

BuRd <- c("#052F60", "#2066AC", "#4392C3", "#92C5DE", "#D1E4EF", "#F7F7F7", "#FCDBC6", "#F4A482", "#D55F4D", "#B1172B", "#67001F")

###############################################################################################################################

## add the Arial font
font_add("Arial", regular = "arial.ttf",
         bold = "arialbd.ttf", italic = "ariali.ttf", bolditalic = "arialbi.ttf")

###############################################################################################################################

# Figure 4B. Violin Plot

genes <- list(c("DLL4", "NOTCH1", "RBPJ"))
genes <- as.data.frame(genes)
gene.names <- list(c("Dll4", "Notch1", "Rbpj"))
gene.names <- as.data.frame(gene.names)

#  Script

myplots <- vector('list', nrow(genes))


for (i in 1:nrow(genes)) {
  
  p21 <- VlnPlot(Liver, features = genes[i, 1], group.by = "Condition", cols = palette_Maca_Vln)+ NoLegend() + theme(text = element_text(family="Arial"), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.title.y=element_blank(), plot.title = element_text(face = "bold.italic"))+ggtitle(gene.names[i, 1])
  
  myplots[[i]] <- local({
    i <- i
    p21
  })
  
}

plegend <- VlnPlot(Liver, features = genes[i, 1], group.by = "Condition") +
   scale_fill_manual(values= palette_Maca_Vln, labels=c("Control", expression(italic("Dll4"^"iDEC")), expression(italic(Notch1)^"iDEC"), expression(italic(Rbpj)^"iDEC")))
plegend
legend <- get_legend(plegend + theme(legend.position = "bottom", legend.justification = "center", legend.direction = "horizontal", legend.text = element_text(size = 18)))

p25 <- patchwork::wrap_plots(myplots, ncol = 1)

p25

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

###############################################################################################################################

# Figure 4C. UMAP

DimPlot(Liver, label = T)

# jpeg("Plots/Figure_4/Figure_4C_Axes.jpeg", width = 10, height = 6, units = 'in', res = 800)
# DimPlot(Liver, cols = my_palette_Rui_colors_B, pt.size = 1.2)+theme(legend.text = element_text(size = 14))
# dev.off()

p1 <- DimPlot(Liver, cols = my_palette_Rui_colors_B, pt.size = 1.2, split.by = "Condition", group.by = "FigClustering", combine = T)


new_labels <- c("Control" = "Control", "Dll4KO" = "italic(Dll4)^iDEC", "NOTCH1KO" = "italic(Notch1)^iDEC", "RBPJKO" = "italic(Rbpj)^iDEC")
p1

# jpeg("Plots/Figure_4/Figure_4C_Axes_Conditions.jpeg", width = 28, height = 6, units = 'in', res = 800)
# p1+facet_grid(~Condition, labeller = as_labeller(new_labels, default = label_parsed))+theme(strip.text.x = element_text(size = 24), legend.text = element_text(size = 14))
# dev.off()

p21 <- p1+facet_grid(~Condition, labeller = as_labeller(new_labels, default = label_parsed))+
  theme(strip.text.x = element_text(size = 24), plot.title = element_blank())+NoLegend()+xlab("UMAP_1")+ylab("UMAP_2")

p22 <- DimPlot(Liver, cols = my_palette_Rui_colors_B, group.by = "FigClustering", pt.size = 1.2)+theme(legend.text = element_text(size = 14), axis.title.x = element_blank(), axis.line.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks = element_blank(), plot.title = element_blank())

p21
p22

p23 <- list(p21, p22)

design <- c(patchwork::area(1, 1, 1, 4), patchwork::area(1, 5, 1, 5.5))

p24 <- Reduce( `+`,  p23)+patchwork::plot_layout(design = design)


cairo_pdf("Plots/Figure_4/Figure_4C1_UMAP.pdf",  width = 34, height = 6, family = "Arial")
p24
dev.off()


jpeg("Plots/Figure_4/Figure_4C1_UMAP.jpeg", width = 34, height = 6, units = 'in', res = 800)
p24
dev.off()

###############################################################################################################################

# Figure 4D. Bar Plot

p1 <- dittoBarPlot(Liver, var = "FigClustering", group.by = "Condition", color.panel = my_palette_Rui_colors_B, xlab = NULL, main = NULL, y.breaks = c(0,0.25,0.5, 0.75, 1), min = 0, max = 1)+
  theme(axis.text=element_text(size=15), axis.title = element_text(size = 18), legend.text = element_text(size = 16), axis.text.x = element_text(angle = 25))+
  scale_x_discrete(labels=c("Control", bquote(italic(Dll4)^iDEC), bquote(italic(Notch1)^iDEC), bquote(italic(Rbpj)^iDEC)))+
  geom_col(width = 0.1)+
  # scale_y_discrete(labels=c("0", "25", "50", "75", "100"))
  scale_y_continuous(labels = c("0", "25", "50", "75", "100"))
p1

cairo_pdf("Plots/Figure_4/Figure_4C2_BarPlot.pdf",  width = 7, height = 6, family = "Arial")
p1
dev.off()


jpeg("Plots/Figure_4/Figure_4C2_BarPlot.jpeg", width = 7, height = 6, units = 'in', res = 800)
p1
dev.off()

###############################################################################################################################

# Figure4D.  DotPlot with Top10 Upregulated genes per cluster. I will keep going with the wilcox test for the DEG analysis, as it gives a better idea of the markers that give rise to the Tip cell cluster

DEG_FigClustering <- FindAllMarkers(Liver, logfc.threshold = 0.5, test.use = "wilcox", only.pos = T, verbose = T, min.diff.pct = -Inf)
DEG_FigClustering$order <- 1:nrow(DEG_FigClustering)
write.xlsx(DEG_FigClustering, "Tables/DEG_FigClustering.wilcox.xlsx", rowNames = T)
DEG_FigClustering.wilcox <- read_excel("Tables/DEG_FigClustering.wilcox.xlsx")


DEG_FigClustering <- as.data.frame(DEG_FigClustering.wilcox)
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

Top10genesR <- rev(Top10genes) 

Top10genesR

figclusteringR <- c("C0 - Unspecified quiescent capillaries", "C1a - Arterial capillaries", "C1v - Venous capillaries", "C2a - Large arteries", "C2v - Large veins",  "C3 - Activated capillaries", 
                    "C4 - Endothelial tip cells", "C5ip - Proliferating S-phase", "C5p - G2-M", "C6 - X") %>% rev()

figclusteringR



levels(Liver) <- figclusteringR


Top10genesR <- rev(Top10genes)

Top10genesDotPlot <- unique(Top10genes) %>% tolower() %>% str_to_title() %>% rev()

Top10genesDotPlotH <- unique(Top10genes) %>% tolower() %>% str_to_title()

cairo_pdf("Plots/Figure_4/Figure_4D_DotPlot_Horizontal.pdf", width = 22, height = 3, family = "Arial")
DotPlot(Liver, features = Top10genes, col.min = 0, dot.scale = 5, scale = T)+
  theme(strip.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"),
        legend.title = element_text(size = 10), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm")) +
  scale_x_discrete(label = Top10genesDotPlotH)+
  scale_colour_gradientn(colours = Bestholtz_palette)
dev.off()


jpeg("Plots/Figure_4/Figure_4D_DotPlot_Horizontal.jpeg", width = 22, height = 3, units = 'in', res = 800)
DotPlot(Liver, features = Top10genes, col.min = 0, dot.scale = 5, scale = T)+
  theme(strip.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"),
        legend.title = element_text(size = 10), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm")) +
  scale_x_discrete(label = Top10genesDotPlotH)+
  scale_colour_gradientn(colours = Bestholtz_palette)
dev.off()

#  DotPlot with Macarena selected genes

Dotplot_genes <- read.csv2("DotplotGenes_MFC.csv", sep = ",", header = F)

Dotplot_genes <- Dotplot_genes$V1

Idents(Liver) <- Liver@meta.data$FigClustering

Dotplot_genesR <- rev(Dotplot_genes) 

Dotplot_genesR

figclusteringR <- c("C0 - Unspecified quiescent capillaries", "C1a - Arterial capillaries", "C1v - Venous capillaries", "C2a - Large arteries", "C2v - Large veins",  "C3 - Activated capillaries", 
                    "C4 - Endothelial tip cells", "C5ip - Proliferating S-phase", "C5p - G2-M", "C6 - X") %>% rev()

figclusteringR



levels(Liver) <- figclusteringR


DotPlotGenesH <- unique(Dotplot_genes) %>% tolower() %>% str_to_title()

df <- data.frame(
  x = rep(c(5.5, 15.5, 25.5, 35.5, 45.5, 55.5, 65.5, 75.5, 85.5, 95.5)),
  y = rep(c(11), 1),
  z = factor(rep(1:10)),
  color = my_palette_Rui_colors)

jpeg("Plots/Figure_4/Figure_4D_DotPlot_Maca_genes.Final.jpeg", width = 22, height = 3, units = 'in', res = 800)
DotPlot(Liver, features = Dotplot_genes, col.min = 0, dot.scale = 5, scale = T)+
  theme(strip.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"),
        legend.title = element_text(size = 10), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"), 
        plot.margin = unit(c(0, 0, 0, 0), "cm")
  ) +
  scale_x_discrete(label = DotPlotGenesH)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  geom_vline(xintercept = c(10.5, 20.5, 30.5, 40.5, 50.5, 60.5, 70.5, 80.5, 90.5),linetype = 2 )+
  geom_raster(df, mapping = aes(x, y, fill = factor(z)), show.legend = F)+
  scale_fill_manual(values = my_palette_Rui_colors_B)
dev.off()

cairo_pdf("Plots/Figure_4/Figure_4D_DotPlot_Maca_genes.Final.pdf", width = 22, height = 3, family = "Arial")
DotPlot(Liver, features = Dotplot_genes, col.min = 0, dot.scale = 5, scale = T)+
  theme(strip.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"),
        legend.title = element_text(size = 10), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"), 
        plot.margin = unit(c(0, 0, 0, 0), "cm")
  ) +
  scale_x_discrete(label = DotPlotGenesH)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  geom_vline(xintercept = c(10.5, 20.5, 30.5, 40.5, 50.5, 60.5, 70.5, 80.5, 90.5),linetype = 2 )+
  geom_raster(df, mapping = aes(x, y, fill = factor(z)), show.legend = F)+
  scale_fill_manual(values = my_palette_Rui_colors_B)
dev.off()


###############################################################################################################################

# Figure 4E. Violin Plot Hes1


p21 <- VlnPlot(Liver, features = "HES1", group.by = "Condition", cols = palette_Maca_Vln)+ NoLegend() + theme(text = element_text(family="Arial"), 
                                  axis.title.x=element_blank(), axis.text.x=element_blank(), plot.title = element_text(face = "bold.italic"), axis.title.y = element_blank())+
                                  ggtitle("Hes1")

plegend <- VlnPlot(Liver, features = "HES1", group.by = "Condition") +
  scale_fill_manual(values= palette_Maca_Vln, labels=c("Control", expression(italic("Dll4"^"iDEC")), expression(italic(Notch1)^"iDEC"), expression(italic(Rbpj)^"iDEC")))
plegend
legend <- get_legend(plegend + theme(legend.position = "bottom", legend.justification = "center", legend.direction = "horizontal", legend.text = element_text(size = 16)))

p26 <- plot_grid(p21, legend, ncol = 1, rel_heights = c(1, .1))

x.axis <- ggdraw() + draw_label("Expression Level", fontface='bold', angle = 90, size = 18)

p27 <- plot_grid(x.axis, p26, ncol = 2, rel_widths = c(.05, 1))

p27
cairo_pdf("Plots/Figure_4/Figure_4E_VlnPlot_Hes1.pdf", width = 6, height = 4, family = "Arial")
p27
dev.off()

jpeg("Plots/Figure_4/Figure_4E_VlnPlot_Hes1.jpeg", width = 6, height = 4, units = 'in', res = 800)
p27
dev.off()



###############################################################################################################################

# Figure 4F. Violin Plot Stmn1

p21 <- VlnPlot(Liver, features = "STMN1", group.by = "Condition", cols = palette_Maca_Vln)+ NoLegend() + theme(text = element_text(family="Arial"), 
                axis.title.x=element_blank(), axis.text.x=element_blank(), plot.title = element_text(face = "bold.italic"), axis.title.y = element_blank())+
               ggtitle("Stmn1")

plegend <- VlnPlot(Liver, features = "HES1", group.by = "Condition") +
  scale_fill_manual(values= palette_Maca_Vln, labels=c("Control", expression(italic("Dll4"^"iDEC")), expression(italic(Notch1)^"iDEC"), expression(italic(Rbpj)^"iDEC")))
plegend
legend <- get_legend(plegend + theme(legend.position = "bottom", legend.justification = "center", legend.direction = "horizontal", legend.text = element_text(size = 16)))

p26 <- plot_grid(p21, legend, ncol = 1, rel_heights = c(1, .1))

x.axis <- ggdraw() + draw_label("Expression Level", fontface='bold', angle = 90, size = 18)

p27 <- plot_grid(x.axis, p26, ncol = 2, rel_widths = c(.05, 1))

p27
cairo_pdf("Plots/Figure_4/Figure_4F_VlnPlot_Stmn1.pdf", width = 6, height = 4, family = "Arial")
p27
dev.off()

jpeg("Plots/Figure_4/Figure_4F_VlnPlot_Stmn1.jpeg", width = 6, height = 4, units = 'in', res = 800)
p27
dev.off()



###############################################################################################################################

# Figure 4G. Violin Plot

genes <- list(c("CDKN1A", "TRP53", "CDKN2A"))
genes <- as.data.frame(genes)
gene.names <- list(c("Cdkn1a", "Trp53", "Cdkn2a"))
gene.names <- as.data.frame(gene.names)

#  Script

myplots <- vector('list', nrow(genes))


for (i in 1:nrow(genes)) {
  
  p21 <- VlnPlot(Liver, features = genes[i, 1], group.by = "Condition", cols = palette_Maca_Vln)+ NoLegend() + theme(text = element_text(family="Arial"), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.title.y=element_blank(), plot.title = element_text(face = "bold.italic"))+ggtitle(gene.names[i, 1])
  
  myplots[[i]] <- local({
    i <- i
    p21
  })
  
}

plegend <- VlnPlot(Liver, features = genes[i, 1], group.by = "Condition") +
  scale_fill_manual(values= palette_Maca_Vln, labels=c("Control", expression(italic("Dll4"^"iDEC")), expression(italic(Notch1)^"iDEC"), expression(italic(Rbpj)^"iDEC")))
plegend
legend <- get_legend(plegend + theme(legend.position = "bottom", legend.justification = "center", legend.direction = "horizontal", legend.text = element_text(size = 16)))

p25 <- patchwork::wrap_plots(myplots, ncol = 1)

p26 <- plot_grid(p25, legend, ncol = 1, rel_heights = c(1, .1))

x.axis <- ggdraw() + draw_label("Expression Level", fontface='bold', angle = 90, size = 24)

p27 <- plot_grid(x.axis, p26, ncol = 2, rel_widths = c(.05, 1))

cairo_pdf("Plots/Figure_4/Figure_4G_VlnPlot.pdf", width = 6, height = 9, family = "Arial")
p27
dev.off()

jpeg("Plots/Figure_4/Figure_4G_VlnPlot.jpeg", width = 6, height = 9, units = 'in', res = 800)
p27
dev.off()



###############################################################################################################################

saveRDS(Liver, "./rds/Liver_20Jan21.Figures.rds")

###############################################################################################################################

# SUPLEMENTARY

# FeaturePlots Figure.S4.1

genesFeaturePlots <- c("GJA5", "LTBP4", "MSR1", "RSPO3", "WNT2", "MKI67", "TOP2A", "STMN1", "RPL32", "RPL10A", "MALAT1", "MEIS1", "KCNE3", "ODC1", "CTNNB1", "RPL21")
Artery_Vein <- c("GJA5", "LTBP4", "MSR1", "RSPO3", "WNT2")
Artery_Vein <- as.data.frame(Artery_Vein)
Artery_Vein_text <- c("Gja5", "Ltbp4", "Msr1", "Rspo3", "Wnt2")

myplots <- vector('list', nrow(Artery_Vein))

for (i in 1:nrow(Artery_Vein)) {
  
  p3 <- FeaturePlot(Liver, features = Artery_Vein[i, 1], order = T, pt.size = 1.2, slot = "data", combine = F)
  
  p4 <- lapply(X = p3, FUN = function(p) p + scale_colour_gradientn(colors = Bestholtz_palette) + ggtitle(bquote(~italic(.(Artery_Vein_text[i])))))
  
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


cairo_pdf("Plots/Figure_4/Suplementary/Figure_S4.1_FeaturePlot.Artery_Vein.pdf", width = 30, height = 6, family = "Arial")
p25
dev.off()

jpeg("Plots/Figure_4/Suplementary/Figure_S4.1_FeaturePlot.Artery_Vein.jpeg", width = 30, height = 6, units = 'in', res = 800)
p25
dev.off()

# Tip Cells + Proliferating

TipCell_Proliferating <- c("KCNE3", "ODC1", "MKI67", "TOP2A", "STMN1")
TipCell_Proliferating <- as.data.frame(TipCell_Proliferating)
TipCell_Proliferating_text <- c("Kcne3", "Odc1", "Mki67", "Top2a", "Stmn1")

myplots <- vector('list', nrow(TipCell_Proliferating))

for (i in 1:nrow(TipCell_Proliferating)) {
  
  p3 <- FeaturePlot(Liver, features = TipCell_Proliferating[i, 1], order = T, pt.size = 1.2, slot = "data", combine = F)
  
  p4 <- lapply(X = p3, FUN = function(p) p + scale_colour_gradientn(colors = Bestholtz_palette) + ggtitle(bquote(~italic(.(TipCell_Proliferating_text[i])))))
  
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


cairo_pdf("Plots/Figure_4/Suplementary/Figure_S4.1_FeaturePlot.TipCell_Proliferating.pdf", width = 30, height = 6, family = "Arial")
p25
dev.off()

jpeg("Plots/Figure_4/Suplementary/Figure_S4.1_FeaturePlot.TipCell_Proliferating.jpeg", width = 30, height = 6, units = 'in', res = 800)
p25
dev.off()

# Quiescent+Activated

Quiescent_Activated <- c("MALAT1", "MEIS1", "RPL32", "RPL10A")
Quiescent_Activated <- as.data.frame(Quiescent_Activated)
Quiescent_Activated_text <- c("Malat1", "Meis1", "Rpl32", "Rpl10a")

myplots <- vector('list', nrow(Quiescent_Activated))

for (i in 1:nrow(Quiescent_Activated)) {
  
  p3 <- FeaturePlot(Liver, features = Quiescent_Activated[i, 1], order = T, pt.size = 1.2, slot = "data", combine = F)
  
  p4 <- lapply(X = p3, FUN = function(p) p + scale_colour_gradientn(colors = Bestholtz_palette) + ggtitle(bquote(~italic(.(Quiescent_Activated_text[i])))))
  
  p4 <- Reduce( `+`, p4 )+patchwork::plot_layout( ncol = 1 )
  
  myplots[[i]] <- local({
    i <- i
    p4
  })
  
  # plegend <- VlnPlot(Liver, features = genes[i, 1], group.by = "Condition")
}


myplots[3]

p25 <- patchwork::wrap_plots(myplots, ncol = 4)

p25


cairo_pdf("Plots/Figure_4/Suplementary/Figure_S4.1_FeaturePlot.Quiescent_Activated.pdf", width = 24, height = 6, family = "Arial")
p25
dev.off()

jpeg("Plots/Figure_4/Suplementary/Figure_S4.1_FeaturePlot.Quiescent_Activated.jpeg", width = 24, height = 6, units = 'in', res = 800)
p25
dev.off()

# C6 - X
C6 <- c("CTNNB1", "RPL21")
C6 <- as.data.frame(C6)
C6_text <- c("Ctnnb1", "Rpl21")

myplots <- vector('list', nrow(C6))

for (i in 1:nrow(C6)) {
  
  p3 <- FeaturePlot(Liver, features = C6[i, 1], order = T, pt.size = 1.2, slot = "data", combine = F)
  
  p4 <- lapply(X = p3, FUN = function(p) p + scale_colour_gradientn(colors = Bestholtz_palette) + ggtitle(bquote(~italic(.(C6_text[i])))))
  
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


cairo_pdf("Plots/Figure_4/Suplementary/Figure_S4.1_FeaturePlot.C6.pdf", width = 12, height = 6, family = "Arial")
p25
dev.off()

jpeg("Plots/Figure_4/Suplementary/Figure_S4.1_FeaturePlot.C6.jpeg", width = 12, height = 6, units = 'in', res = 800)
p25
dev.off()

###############################################################################################################################

# Figure 4S.2 Tip Cell Cluster Top 20 genes Horizontal

DEG_FigClustering.wilcox <- read_excel("Tables/DEG_FigClustering.wilcox.xlsx")

Top20_C4_Tip_Cells <- filter(DEG_FigClustering.wilcox, cluster == "C4 - Endothelial tip cells") %>% top_n(n = 20, wt = -order)

write.xlsx(Top20_C4_Tip_Cells, "Tables/Top20_C4_Tip_Cells.wilcox.xlsx", rowNames = T)

Top20_C4_Tip_Cells_wilcox <- read_excel("Tables/Top20_C4_Tip_Cells.wilcox.xlsx")

Top20_C4_Tip_Cells_genes <- Top20_C4_Tip_Cells_wilcox$...2[1:16] %>% rev()

Top20_C4_Tip_Cells_genes <- append(Top20_C4_Tip_Cells_genes, c("KCNE3", "ODC1", "PEG10", "ESM1", "APLN", "VEGFA", "TPI1" , "MIF", "ENO1"), after = 0) %>% unique()

Top20_C4_Tip_Cells_genes

Top20genesC4DotPlot <- unique(Top20_C4_Tip_Cells_genes) %>% tolower() %>% str_to_title()

figclusteringR <- c("C0 - Unspecified quiescent capillaries", "C1a - Arterial capillaries", "C1v - Venous capillaries", "C2a - Large arteries", "C2v - Large veins",  "C3 - Activated capillaries", 
                    "C4 - Endothelial tip cells", "C5ip - Proliferating S-phase", "C5p - G2-M", "C6 - X") %>% rev()

figclusteringR



levels(Liver) <- figclusteringR


cairo_pdf("Plots/Figure_4/Suplementary/Figure_S4.2_C4_tip_cells_DotPlot_Top20_Horizontal.pdf", width = 9, height = 3, family = "Arial")
DotPlot(Liver, features = Top20_C4_Tip_Cells_genes, col.min = 0, dot.scale = 5, scale = T)+
  theme(strip.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x.bottom = element_text(face = "italic"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.title = element_text(size = 10), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"))+
  scale_x_discrete(label = Top20genesC4DotPlot)+
  scale_colour_gradientn(colours = Bestholtz_palette)
dev.off()

jpeg("Plots/Figure_4/Suplementary/Figure_S4.2_C4_tip_cells_DotPlot_Top20_Horizontal.jpeg", width = 9, height = 3, units = "in", family = "Arial", res = 800)
DotPlot(Liver, features = Top20_C4_Tip_Cells_genes, col.min = 0, dot.scale = 5, scale = T)+
  theme(strip.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x.bottom = element_text(face = "italic"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.title = element_text(size = 10), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"))+
  scale_x_discrete(label = Top20genesC4DotPlot)+
  scale_colour_gradientn(colours = Bestholtz_palette)
dev.off()
