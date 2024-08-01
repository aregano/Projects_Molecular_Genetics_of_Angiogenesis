# Libraries
library(Seurat)
library(ggplot2)
library(patchwork)
library("dittoSeq")
library("magick")
library("RColorBrewer")
library(cowplot)
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

Liver <- readRDS("S:/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Liver_MacarenaFC/Figures_Paper/rds/Liver_20Jan21.Figures.All.Conditions.rds")

Liver4d <- readRDS("S:/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Liver_MacarenaFC/Figures_Paper/rds/Liver4d.Mapping.rds")

table(Liver@meta.data$Condition)

table(Liver4d@meta.data$Condition)

# Color Palettes

my_palette_Rui_colors <- c("#E95A74", "#50B6EF", "#45FF8E", "#F4A753", "#A80519", "#880088", "#E28CF4", "#C1B80C","#FC0808",  "#0E47D8")

my_palette_Rui_colors_B <- c("#F4A753", "#E95A74", "#50B6EF", "#FC0808",  "#0E47D8", "#45FF8E", "#A80519", "#E28CF4", "#880088", "#C1B80C")

custom_rdylblu_palette <- c("#005AD9", "#0070FF", "#FFE135", "#FFDA13", "#E60026", "#C60012")

macarena_palette <- c("#FDDC85", "#6E2914")

Bestholtz_palette <- c("#DEDAD6", "#FEE392", "#FEC44E", "#FE9929", "#ED6F11", "#CC4C17", "#993411", "#65260C")

palette_Maca_Vln <- c("#606060", "#F94040", "#05BE78", "#FFA040")

palette_Maca_Vln_2 <- c("#606060", "lightyellow3", "#F94040", "lightpink2", "olivedrab2")

BuRd <- c("#052F60", "#2066AC", "#4392C3", "#92C5DE", "#D1E4EF", "#F7F7F7", "#FCDBC6", "#F4A482", "#D55F4D", "#B1172B", "#67001F")

## add the Arial font
font_add("Arial", regular = "arial.ttf",
         bold = "arialbd.ttf", italic = "ariali.ttf", bolditalic = "arialbi.ttf")


##############################################################################################################################

# Figure UMAP

DimPlot(Liver, label = T)

Idents(Liver) <- Liver@meta.data$FigClustering

p1 <- DimPlot(Liver, cols = my_palette_Rui_colors_B, pt.size = 1.2, group.by = "FigClustering", split.by = "Condition", combine = T)

p1

new_labels <- c("Control" = "Control", "Dll4KO" = "italic(Dll4)^iDEC", "Dll4/MycKO" = "italic(Dll4/Myc)^iDEC")

p21 <- p1+facet_grid(~Condition, labeller = as_labeller(new_labels, default = label_parsed))+
  theme(strip.text.x = element_text(size = 24), plot.title = element_blank())+NoLegend()+xlab("UMAP_1")+ylab("UMAP_2")

p22 <- DimPlot(Liver, group.by = "FigClustering", cols = my_palette_Rui_colors_B, pt.size = 1.2)+theme(legend.text = element_text(size = 14), axis.title.x = element_blank(), axis.line.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks = element_blank(), plot.title = element_blank())

p21 <- p1+theme(strip.text.x = element_text(size = 20), plot.title = element_blank())+NoLegend()+xlab("UMAP_1")+ylab("UMAP_2")


p21
p22

p23 <- list(p21, p22)

design <- c(patchwork::area(1, 1, 1, 11), patchwork::area(1, 12, 1, 12.5))

p24 <- Reduce( `+`,  p23)+patchwork::plot_layout(design = design)


cairo_pdf("Plots/Figure_Rui_ppt/UMAP_All_Conditions.pdf",  width = 75, height = 6, family = "Arial")
p24
dev.off()


jpeg("Plots/Figure_Rui_ppt/UMAP_All_Conditions.jpeg", width = 75, height = 6, units = 'in', res = 800)
p24
dev.off()

##############################################################################################################################

# Figure UMAP All Conditions Including 4 days deletion

DimPlot(Liver, label = T)

Idents(Liver) <- Liver@meta.data$FigClustering

p1 <- DimPlot(Liver, cols = my_palette_Rui_colors_B, pt.size = 1.2, group.by = "FigClustering", split.by = "Condition", combine = T)

p1

# p2 <- DimPlot(Liver4d, cols = my_palette_Rui_colors_B, pt.size = 1.2, group.by = "FigClustering", reduction = "ref.umap", split.by = "Condition", combine = T)

p2

# new_labels <- c("Control" = "Control", "Dll4KO" = "italic(Dll4)^iDEC", "Dll4/MycKO" = "italic(Dll4/Myc)^iDEC")

p21 <- p1+theme(strip.text.x = element_text(size = 24), plot.title = element_blank())+NoLegend()+xlab("UMAP_1")+ylab("UMAP_2")

p23 <- p2+theme(strip.text.x = element_text(size = 24), plot.title = element_blank(),axis.title.x = element_blank(), axis.line.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks = element_blank())+NoLegend()+xlab("UMAP_1")+ylab("UMAP_2")


p22 <- DimPlot(Liver, group.by = "FigClustering", cols = my_palette_Rui_colors_B, pt.size = 1.2)+theme(legend.text = element_text(size = 14), axis.title.x = element_blank(), axis.line.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks = element_blank(), plot.title = element_blank())

p21
p22
p23

p24 <- list(p21, p22)

design <- c(patchwork::area(1, 1, 1, 14), patchwork::area(1, 15, 1, 15.5))

p25 <- Reduce( `+`,  p24)+patchwork::plot_layout(design = design)

p25

cairo_pdf("Plots/Figure_Rui_ppt/UMAP_All_Conditions_New.pdf",  width = 93, height = 6, family = "Arial")
p25
dev.off()


jpeg("Plots/Figure_Rui_ppt/UMAP_All_Conditions_New.jpeg", width = 93, height = 6, units = 'in', res = 600)
p25
dev.off()

###############################################################################################################################

# Figure ppt Rui Bar Plot

table(Liver@meta.data$Condition)

p1 <- dittoBarPlot(Liver, var = "FigClustering", group.by = "Condition",  x.reorder = c(1,3,4,5,6,9,10,11,12,13,14,2,8,7),color.panel = my_palette_Rui_colors_B, xlab = NULL, main = NULL, y.breaks = c(0,0.25,0.5, 0.75, 1), min = 0, max = 1)+
  theme(axis.text=element_text(size=15), axis.title = element_text(size = 18), legend.text = element_text(size = 16), axis.text.x = element_text(angle = 25))+
  # scale_x_discrete(labels=c("Control", bquote(italic(Dll4)^iDEC), bquote(italic(Dll4Het)^iDEC), bquote(italic(DBZ)^iDEC)))+
  geom_col(width = 0.1)+
  # scale_y_discrete(labels=c("0", "25", "50", "75", "100"))
  scale_y_continuous(labels = c("0", "25", "50", "75", "100"))
p1

cairo_pdf("Plots/Figure_Rui_ppt/Figure_BarPlot_All_Conditions.pdf",  width = 12, height = 6, family = "Arial")
p1
dev.off()


jpeg("Plots/Figure_Rui_ppt/Figure_BarPlot_All_Conditions.jpeg", width = 12, height = 6, units = 'in', res = 800)
p1
dev.off()

###############################################################################################################################

# Figure ppt Rui. Violin Plot Hes1. WT, Ctrl2w, Ctrl4d, Dll4KO, D1D4Het, WT_DBZ

Liver <- subset(Liver, subset = Condition == "Control" | Condition == "WT" | Condition == "WT_DBZ" | Condition == "Dll4KO" | Condition == "D1D4HET")

levels(Liver) <- c("Control", "WT", "Dll4KO", "D1D4HET", "WT_DBZ")

Liver@active.ident -> Liver@meta.data$Condition

p27 <- VlnPlot(Liver, features = "HES1", group.by = "Condition") +  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

p27
p21 <- VlnPlot(Liver, features = "HES1", group.by = "Condition", cols = palette_Maca_Vln_2)+ NoLegend() + theme(text = element_text(family="Arial"), axis.title.x=element_blank(), axis.text.x=element_blank(), plot.title = element_text(face = "bold.italic"), axis.title.y = element_blank())+ggtitle("Hes1")

p21

plegend <- VlnPlot(Liver, features = "HES1", group.by = "Condition") +
  scale_fill_manual(values= palette_Maca_Vln_2, labels=c("Control", "Control(wt)", expression(italic("Dll4")^"iDEC"), expression(italic(Dll4Het)^"iDEC"), expression(italic(DBZ))))+
  guides(colour = guide_legend(nrow = 1, byrow = T))
plegend
legend <- get_legend(plegend + guides(fill = guide_legend(nrow = 1, byrow = T)) + theme(legend.box.background = element_blank(), legend.background = element_blank(), legend.justification = "center", legend.direction = "horizontal", legend.text = element_text(size = 12, hjust = 0)))

legend

p26 <- plot_grid(p21, legend, ncol = 1, rel_heights = c(1, .1))

x.axis <- ggdraw() + draw_label("Expression Level", fontface='bold', angle = 90, size = 18)

p27 <- plot_grid(x.axis, p26, ncol = 2, rel_widths = c(.05, 1))

p27
cairo_pdf("Plots/Figure_Rui_ppt/Figure_ppt_VlnPlot_Hes1_new.pdf", width = 7, height = 4, family = "Arial")
p27
dev.off()

jpeg("Plots/Figure_Rui_ppt/Figure_ppt_VlnPlot_Hes1_new.jpeg", width = 7, height = 4, units = 'in', res = 800)
p27
dev.off()
