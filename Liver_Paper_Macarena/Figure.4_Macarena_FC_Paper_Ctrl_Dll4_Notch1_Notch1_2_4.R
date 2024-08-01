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

Liver <- subset(Liver, subset = Condition == "Control" | Condition == "Dll4KO" | Condition == "NOTCH1KO" | Condition == "N1N4HET")

levels(Liver) <- c("Control", "Dll4KO", "NOTCH1KO", "N1N4HET")

Liver@active.ident -> Liver@meta.data$Condition

table(Liver@meta.data$Condition)

# Color Palettes

my_palette_Rui_colors <- c("#E95A74", "#50B6EF", "#45FF8E", "#F4A753", "#A80519", "#880088", "#E28CF4", "#C1B80C","#FC0808",  "#0E47D8")

my_palette_Rui_colors_B <- c("#F4A753", "#E95A74", "#50B6EF", "#FC0808",  "#0E47D8", "#45FF8E", "#A80519", "#E28CF4", "#880088", "#C1B80C")

custom_rdylblu_palette <- c("#005AD9", "#0070FF", "#FFE135", "#FFDA13", "#E60026", "#C60012")

macarena_palette <- c("#FDDC85", "#6E2914")

Bestholtz_palette <- c("#DEDAD6", "#FEE392", "#FEC44E", "#FE9929", "#ED6F11", "#CC4C17", "#993411", "#65260C")

palette_Maca_Vln <- c("#606060", "#F94040", "#05BE78", "#FFA040")

BuRd <- c("#052F60", "#2066AC", "#4392C3", "#92C5DE", "#D1E4EF", "#F7F7F7", "#FCDBC6", "#F4A482", "#D55F4D", "#B1172B", "#67001F")

###############################################################################################################################

## add the Arial font
font_add("Arial", regular = "arial.ttf",
         bold = "arialbd.ttf", italic = "ariali.ttf", bolditalic = "arialbi.ttf")


###############################################################################################################################

# Figure 4M.1 Bar Plot

table(Liver@meta.data$Condition)

p1 <- dittoBarPlot(Liver, var = "FigClustering", group.by = "Condition", x.reorder = c(1,2,4,3), color.panel = my_palette_Rui_colors_B, xlab = NULL, main = NULL, y.breaks = c(0,0.25,0.5, 0.75, 1), min = 0, max = 1)+
  theme(axis.text=element_text(size=15), axis.title = element_text(size = 18), legend.text = element_text(size = 16), axis.text.x = element_text(angle = 25))+
  scale_x_discrete(labels=c("Control", bquote(italic(Dll4)^iDEC), bquote(italic(Notch1)^iDEC), bquote(italic(Notch1/2/4)^iDEC)))+
  geom_col(width = 0.1)+
  # scale_y_discrete(labels=c("0", "25", "50", "75", "100"))
  scale_y_continuous(labels = c("0", "25", "50", "75", "100"))
p1

cairo_pdf("Plots/Figure_4/Figure_4M.1_BarPlot.pdf",  width = 7, height = 6, family = "Arial")
p1
dev.off()


jpeg("Plots/Figure_4/Figure_4M.1_BarPlot.jpeg", width = 7, height = 6, units = 'in', res = 800)
p1
dev.off()

###############################################################################################################################

# Figure 4M.2 UMAP

DimPlot(Liver, label = T)

Idents(Liver) <- Liver@meta.data$FigClustering

# jpeg("Plots/Figure_4/Figure_4C_Axes.jpeg", width = 10, height = 6, units = 'in', res = 800)
# DimPlot(Liver, cols = my_palette_Rui_colors_B, pt.size = 1.2)+theme(legend.text = element_text(size = 14))
# dev.off()

p1 <- DimPlot(Liver, cols = my_palette_Rui_colors_B, pt.size = 1.2, group.by = "FigClustering", split.by = "Condition", combine = T)


new_labels <- c("Control" = "Control", "Dll4KO" = "italic(Dll4)^iDEC", "NOTCH1KO" = "italic(Notch1)^iDEC", "N1N4HET" = "italic(Notch1/2/4)^iDEC")
p1

# jpeg("Plots/Figure_4/Figure_4C_Axes_Conditions.jpeg", width = 28, height = 6, units = 'in', res = 800)
# p1+facet_grid(~Condition, labeller = as_labeller(new_labels, default = label_parsed))+theme(strip.text.x = element_text(size = 24), legend.text = element_text(size = 14))
# dev.off()

p21 <- p1+facet_grid(~Condition, labeller = as_labeller(new_labels, default = label_parsed))+
  theme(strip.text.x = element_text(size = 24), plot.title = element_blank())+NoLegend()

p22 <- DimPlot(Liver, group.by = "FigClustering", cols = my_palette_Rui_colors_B, pt.size = 1.2)+theme(legend.text = element_text(size = 14), axis.title.x = element_blank(), axis.line.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks = element_blank(), plot.title = element_blank())

p21
p22

p23 <- list(p21, p22)

design <- c(patchwork::area(1, 1, 1, 4), patchwork::area(1, 5, 1, 5.5))

p24 <- Reduce( `+`,  p23)+patchwork::plot_layout(design = design)


cairo_pdf("Plots/Figure_4/Figure_4M.2_UMAP.pdf",  width = 34, height = 6, family = "Arial")
p24
dev.off()


jpeg("Plots/Figure_4/Figure_4M.2_UMAP.jpeg", width = 34, height = 6, units = 'in', res = 800)
p24
dev.off()


###############################################################################################################################

# Figure 4N. Violin Plot

genes <- list(c("NOTCH1", "NOTCH2", "NOTCH4"))
genes <- as.data.frame(genes)
gene.names <- list(c("Notch1", "Notch2", "Notch4"))
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
  scale_fill_manual(values= palette_Maca_Vln, labels=c("Control", expression(italic("Dll4"^"iDEC")), expression(italic(Notch1)^"iDEC"), expression(italic(Notch1/2/4)^"iDEC")))
plegend
legend <- get_legend(plegend + theme(legend.position = "bottom", legend.justification = "center", legend.direction = "horizontal", legend.text = element_text(size = 14)))

p25 <- patchwork::wrap_plots(myplots, ncol = 1)

p26 <- plot_grid(p25, legend, ncol = 1, rel_heights = c(1, .07))

x.axis <- ggdraw() + draw_label("Expression Level", fontface='bold', angle = 90, size = 24)

p27 <- plot_grid(x.axis, p26, ncol = 2, rel_widths = c(.1, 1))

cairo_pdf("Plots/Figure_4/Figure_4N_VlnPlot.pdf", width = 6, height = 9, family = "Arial")
p27
dev.off()

jpeg("Plots/Figure_4/Figure_4N_VlnPlot.jpeg", width = 6, height = 9, units = 'in', res = 800)
p27
dev.off()

###############################################################################################################################

# Figure 4O. Violin Plot

genes <- list(c("ESM1", "STMN1"))
genes <- as.data.frame(genes)
gene.names <- list(c("Esm1", "Stmn1"))
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
  scale_fill_manual(values= palette_Maca_Vln, labels=c("Control", expression(italic("Dll4"^"iDEC")), expression(italic(Notch1)^"iDEC"), expression(italic(Notch1/2/4)^"iDEC")))
plegend
legend <- get_legend(plegend + theme(legend.position = "bottom", legend.justification = "center", legend.direction = "horizontal", legend.text = element_text(size = 14)))

p25 <- patchwork::wrap_plots(myplots, ncol = 1)

p26 <- plot_grid(p25, legend, ncol = 1, rel_heights = c(1, .07))

x.axis <- ggdraw() + draw_label("Expression Level", fontface='bold', angle = 90, size = 24)

p27 <- plot_grid(x.axis, p26, ncol = 2, rel_widths = c(.06, 1))

cairo_pdf("Plots/Figure_4/Figure_4O_VlnPlot.pdf", width = 6, height = 6, family = "Arial")
p27
dev.off()

jpeg("Plots/Figure_4/Figure_4O_VlnPlot.jpeg", width = 6, height = 6, units = 'in', res = 800)
p27
dev.off()


###############################################################################################################################

# Figure 4Q. Violin Plot

p21 <- VlnPlot(Liver, features = "KCNE3", group.by = "Condition", cols = palette_Maca_Vln)+ NoLegend() + theme(text = element_text(family="Arial"), 
                                                                                                              axis.title.x=element_blank(), axis.text.x=element_blank(), plot.title = element_text(face = "bold.italic"), axis.title.y = element_blank())+
  ggtitle("Kcne3")

plegend <- VlnPlot(Liver, features = "HES1", group.by = "Condition") +
  scale_fill_manual(values= palette_Maca_Vln, labels=c("Control", expression(italic("Dll4"^"iDEC")), expression(italic(Notch1)^"iDEC"), expression(italic(Notch1/2/4)^"iDEC")))
plegend
legend <- get_legend(plegend + theme(legend.position = "bottom", legend.justification = "center", legend.direction = "horizontal", legend.text = element_text(size = 14)))

p26 <- plot_grid(p21, legend, ncol = 1, rel_heights = c(1, .1))

x.axis <- ggdraw() + draw_label("Expression Level", fontface='bold', angle = 90, size = 18)

p27 <- plot_grid(x.axis, p26, ncol = 2, rel_widths = c(.05, 1))

p27
cairo_pdf("Plots/Figure_4/Figure_4Q_v1_VlnPlot_Kcne3.pdf", width = 6, height = 4, family = "Arial")
p27
dev.off()

jpeg("Plots/Figure_4/Figure_4Q_v1_VlnPlot_Kcne3.jpeg", width = 6, height = 4, units = 'in', res = 800)
p27
dev.off()

###############################################################################################################################

# Figure 4O_v2. Violin Plot

# Myc, Odc1, Cdkn1a, MiKi67, Ltbp4 and Msr1.

genes <- list(c("MYC", "ODC1", "CDKN1A", "MKI67", "LTBP4", "MSR1"))
genes <- as.data.frame(genes)
gene.names <- list(c("Myc", "Odc1", "Cdkn1a", "Mki67", "Ltbp4", "Msr1"))
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
  scale_fill_manual(values= palette_Maca_Vln, labels=c("Control", expression(italic("Dll4"^"iDEC")), expression(italic(Notch1)^"iDEC"), expression(italic(Notch1/2/4)^"iDEC")))
plegend
legend <- get_legend(plegend + theme(legend.position = "bottom", legend.justification = "center", legend.direction = "horizontal", legend.text = element_text(size = 14)))

p25 <- patchwork::wrap_plots(myplots, ncol = 1)

p26 <- plot_grid(p25, legend, ncol = 1, rel_heights = c(1, .02))

x.axis <- ggdraw() + draw_label("Expression Level", fontface='bold', angle = 90, size = 24)

p27 <- plot_grid(x.axis, p26, ncol = 2, rel_widths = c(.06, 1))

cairo_pdf("Plots/Figure_4/Figure_4Q_v2_VlnPlot.pdf", width = 6, height = 18, family = "Arial")
p27
dev.off()

jpeg("Plots/Figure_4/Figure_4Q_v2_VlnPlot.jpeg", width = 6, height = 18, units = 'in', res = 800)
p27
dev.off()



###############################################################################################################################

# saveRDS(Liver, "./rds/Liver_20Jan21.Figures.Ctrl.Dll4Ko.Notch1KO.Notch1,2,4KO.rds")

#########################################SUPPLEMENTARY######################################################################################

# Figure S4.4E Violin Plot

genes <- list(c("HES1"))
genes <- as.data.frame(genes)
gene.names <- list(c("Hes1"))
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
  scale_fill_manual(values= palette_Maca_Vln, labels=c("Control", expression(italic("Dll4"^"iDEC")), expression(italic(Notch1)^"iDEC"), expression(italic(Notch1/2/4)^"iDEC")))
plegend
legend <- get_legend(plegend + theme(legend.position = "bottom", legend.justification = "center", legend.direction = "horizontal", legend.text = element_text(size = 14)))

p25 <- patchwork::wrap_plots(myplots, ncol = 1)

p26 <- plot_grid(p25, legend, ncol = 1, rel_heights = c(1, .1))

x.axis <- ggdraw() + draw_label("Expression Level", fontface='bold', angle = 90, size = 12)

p27 <- plot_grid(x.axis, p26, ncol = 2, rel_widths = c(.05, 1))

p27

cairo_pdf("Plots/Figure_4/Suplementary/Figure_S4.4D_VlnPlot.pdf", width = 6, height = 3, family = "Arial")
p27
dev.off()

jpeg("Plots/Figure_4/Suplementary/Figure_S4.4D_VlnPlot.jpeg", width = 6, height = 3, units = 'in', res = 800)
p27
dev.off()

