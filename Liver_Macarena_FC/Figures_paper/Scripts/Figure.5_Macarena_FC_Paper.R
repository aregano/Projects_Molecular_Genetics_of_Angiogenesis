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

Liver <- readRDS("S:/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Liver_MacarenaFC/Figures_Paper/rds/Liver_20Jan21.Figures.All.Conditions.OldDll4Myc.rds")

table(Liver@meta.data$Condition)

LiverC <- subset(Liver, subset = Condition == "Control")

Liver <- subset(Liver, subset = Condition == "Control" | Condition == "Jag1/Jag2/Dll1KO")

# Color Palettes

my_palette_Rui_colors <- c("#E95A74", "#50B6EF", "#45FF8E", "#F4A753", "#A80519", "#880088", "#E28CF4", "#C1B80C","#FC0808",  "#0E47D8")

my_palette_Rui_colors_B <- c("#F4A753", "#E95A74", "#50B6EF", "#FC0808",  "#0E47D8", "#45FF8E", "#A80519", "#E28CF4", "#880088", "#C1B80C")

custom_rdylblu_palette <- c("#005AD9", "#0070FF", "#FFE135", "#FFDA13", "#E60026", "#C60012")

macarena_palette <- c("#FDDC85", "#6E2914")

Bestholtz_palette <- c("#DEDAD6", "#FEE392", "#FEC44E", "#FE9929", "#ED6F11", "#CC4C17", "#993411", "#65260C")

palette_Maca_Vln <- c("#606060", "#DBB49A")

palette_Maca_Vln_gene <- c("#A00000", "#FCB4B4", "#F2F0F0", "#F2F0F0")

BuRd <- c("#052F60", "#2066AC", "#4392C3", "#92C5DE", "#D1E4EF", "#F7F7F7", "#FCDBC6", "#F4A482", "#D55F4D", "#B1172B", "#67001F")

###############################################################################################################################

## add the Arial font
font_add("Arial", regular = "arial.ttf",
         bold = "arialbd.ttf", italic = "ariali.ttf", bolditalic = "arialbi.ttf")


###############################################################################################################################

# Figure 5A. Violin Plot

genes <- list(c("DLL4", "DLL1", "JAG1", "JAG2"))
genes <- as.data.frame(genes)
gene.names <- list(c("Dll4", "Dll1", "Jag1", "Jag2"))
gene.names <- as.data.frame(gene.names)

#  Script

myplots <- vector('list', nrow(genes))


for (i in 1:nrow(genes)) {
  
  p21 <- VlnPlot(LiverC, features = genes[i, 1], group.by = "Condition", cols = palette_Maca_Vln_gene[i])+ NoLegend() + theme(text = element_text(family="Arial"), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.title.y=element_blank(), plot.title = element_text(face = "bold.italic"))+ggtitle(gene.names[i, 1])
  
  myplots[[i]] <- local({
    i <- i
    p21
  })
  
}

plegend <- VlnPlot(LiverC, features = genes[i, 1], group.by = "Condition") +
  scale_fill_manual(values= palette_Maca_Vln, labels=c("Control"))
plegend
legend <- get_legend(plegend + theme(legend.position = "bottom", legend.justification = "center", legend.direction = "horizontal", legend.text = element_text(size = 14)))

p25 <- patchwork::wrap_plots(myplots, ncol = 1)

p26 <- plot_grid(p25, legend, ncol = 1, rel_heights = c(1, .07))

x.axis <- ggdraw() + draw_label("Expression Level", fontface='bold', angle = 90, size = 24)

p27 <- plot_grid(x.axis, p26, ncol = 2, rel_widths = c(.3, 1))

cairo_pdf("Plots/Figure_5/Figure_5A_VlnPlot.pdf", width = 3, height = 12, family = "Arial")
p27
dev.off()

jpeg("Plots/Figure_5/Figure_5A_VlnPlot.jpeg", width = 3, height = 12, units = 'in', res = 800)
p27
dev.off()

###############################################################################################################################

# Figure 5E. Violin Plot

genes <- list(c("DLL4", "DLL1", "JAG1", "JAG2"))
genes <- as.data.frame(genes)
gene.names <- list(c("Dll4", "Dll1", "Jag1", "Jag2"))
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
  scale_fill_manual(values= palette_Maca_Vln, labels=c("Control", expression(italic(Jag1/2)^"iDEC")))
plegend
legend <- get_legend(plegend + theme(legend.position = "bottom", legend.justification = "center", legend.direction = "horizontal", legend.text = element_text(size = 14)))

p25 <- patchwork::wrap_plots(myplots, ncol = 1)

p26 <- plot_grid(p25, legend, ncol = 1, rel_heights = c(1, .07))

x.axis <- ggdraw() + draw_label("Expression Level", fontface='bold', angle = 90, size = 24)

p27 <- plot_grid(x.axis, p26, ncol = 2, rel_widths = c(.2, 1))

cairo_pdf("Plots/Figure_5/Figure_5E.1_VlnPlot.pdf", width = 4, height = 12, family = "Arial")
p27
dev.off()

jpeg("Plots/Figure_5/Figure_5E.1_VlnPlot.jpeg", width = 4, height = 12, units = 'in', res = 800)
p27
dev.off()

###############################################################################################################################

# Figure 5E.2. Violin Plot

genes <- list(c("EFNB2", "WNT2", "CD34", "HES1"))
genes <- as.data.frame(genes)
gene.names <- list(c("Efnb2", "Wnt2", "CD34", "Hes1"))
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
  scale_fill_manual(values= palette_Maca_Vln, labels=c("Control", expression(italic(Jag1/2)^"iDEC")))
plegend
legend <- get_legend(plegend + theme(legend.position = "bottom", legend.justification = "center", legend.direction = "horizontal", legend.text = element_text(size = 14)))

p25 <- patchwork::wrap_plots(myplots, ncol = 1)

p26 <- plot_grid(p25, legend, ncol = 1, rel_heights = c(1, .07))

x.axis <- ggdraw() + draw_label("Expression Level", fontface='bold', angle = 90, size = 24)

p27 <- plot_grid(x.axis, p26, ncol = 2, rel_widths = c(.1, 1))

cairo_pdf("Plots/Figure_5/Figure_5E.2_VlnPlot.pdf", width = 4, height = 12, family = "Arial")
p27
dev.off()

jpeg("Plots/Figure_5/Figure_5E.2_VlnPlot.jpeg", width = 4, height = 12, units = 'in', res = 800)
p27
dev.off()

###############################################################################################################################
