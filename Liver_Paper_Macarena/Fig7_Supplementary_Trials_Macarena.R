library(Seurat)
library(dittoSeq)
library(ggplot2)

# Color Palettes

my_palette_Rui_colors <- c("#E95A74", "#50B6EF", "#45FF8E", "#F4A753", "#A80519", "#880088", "#E28CF4", "#C1B80C","#FC0808",  "#0E47D8")

my_palette_Rui_colors_B <- c("#F4A753", "#E95A74", "#50B6EF", "#FC0808",  "#0E47D8", "#45FF8E", "#A80519", "#E28CF4", "#880088", "#C1B80C")

custom_rdylblu_palette <- c("#005AD9", "#0070FF", "#FFE135", "#FFDA13", "#E60026", "#C60012")

macarena_palette <- c("#FDDC85", "#6E2914")

Bestholtz_palette <- c("#DEDAD6", "#FEE392", "#FEC44E", "#FE9929", "#ED6F11", "#CC4C17", "#993411", "#65260C")

palette_Maca_Vln <- c("#606060", "#F94040", "#BB005E")

BuRd <- c("#052F60", "#2066AC", "#4392C3", "#92C5DE", "#D1E4EF", "#F7F7F7", "#FCDBC6", "#F4A482", "#D55F4D", "#B1172B", "#67001F")


palette_Maca_GSEA <- c("#F94040", "mediumpurple1", "#BB005E")

####################################################################

# Supplementary 7.1i Ctrl 4d, Dll4KO4d, Dll4KO2w

Liver <- readRDS("S:/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Liver_MacarenaFC/Figures_Paper/rds/Liver_20Jan21.Figures.Fig7Sup.1.rds")

levels(Liver@meta.data$Condition)

Liver <- subset(Liver, subset = Condition == "Control(4d)" | Condition == "Dll4KO(4d)+Vehicle" | Condition == "Dll4KO")

levels(Liver) <- c("Control(4d)", "Dll4KO(4d)+Vehicle", "Dll4KO")

Liver@meta.data$Condition <- Liver@active.ident

DimPlot(Liver)

p1 <- DimPlot(Liver, cols = my_palette_Rui_colors_B, pt.size = 1.2, split.by = "Condition", group.by = "FigClustering", combine = T)


new_labels <- c("Control(4d)" = "Control4d", "Dll4KO(4d)+Vehicle" = "italic(Dll4)^iDEC4d+Vehicle", "Dll4KO" = "italic(Dll4)^iDEC2w")
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

design <- c(patchwork::area(1, 1, 1, 3), patchwork::area(1, 4, 1, 4.5))

p24 <- Reduce( `+`,  p23)+patchwork::plot_layout(design = design)


cairo_pdf("Plots/Figure_7/Supplementary/Trials/Figure_S7_1i_UMAP.pdf",  width = 28, height = 6, family = "Arial")
p24
dev.off()


jpeg("Plots/Figure_7/Supplementary/Trials/Figure_S7_1i_UMAP.jpeg", width = 28, height = 6, units = 'in', res = 800)
p24
dev.off()

###################################################################


# Supplementary S7.1i BarPlot

p1 <- dittoBarPlot(Liver, var = "FigClustering", group.by = "Condition", color.panel = my_palette_Rui_colors_B, xlab = NULL, main = NULL, x.reorder = c(1,3,2), y.breaks = c(0,0.25,0.5, 0.75, 1), min = 0, max = 1)+
  theme(axis.text=element_text(size=15), axis.title = element_text(size = 18), legend.text = element_text(size = 16), axis.text.x = element_text(angle = 25))+
  scale_x_discrete(labels=c("Control4d", bquote(italic(Dll4)^iDEC4d+Veh), bquote(italic(Dll4)^iDEC2w)))+
  geom_col(width = 0.1)+
  # scale_y_discrete(labels=c("0", "25", "50", "75", "100"))
  scale_y_continuous(labels = c("0", "25", "50", "75", "100"))
p1

cairo_pdf("Plots/Figure_7/Supplementary/Trials/Figure_S7_1i_BarPlot.pdf",  width = 7, height = 6, family = "Arial")
p1
dev.off()


jpeg("Plots/Figure_7/Supplementary/Trials/Figure_S7_1i_BarPlot.jpeg", width = 7, height = 6, units = 'in', res = 800)
p1
dev.off()


####################################################################

# Supplementary 7.1ii Ctrl 4d, Dll4KO4d, Dll4KO2w

Liver <- readRDS("S:/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Liver_MacarenaFC/Figures_Paper/rds/Liver_20Jan21.Figures.Fig7Sup.1.rds")

levels(Liver@meta.data$Condition)

Liver <- subset(Liver, subset = Condition == "Control" | Condition == "Dll4KO(4d)+Vehicle" | Condition == "Dll4KO")

levels(Liver) <- c("Control", "Dll4KO(4d)+Vehicle", "Dll4KO")

Liver@meta.data$Condition <- Liver@active.ident

DimPlot(Liver)

p1 <- DimPlot(Liver, cols = my_palette_Rui_colors_B, pt.size = 1.2, split.by = "Condition", group.by = "FigClustering", combine = T)


new_labels <- c("Control" = "Control", "Dll4KO(4d)+Vehicle" = "italic(Dll4)^iDEC4d+Vehicle", "Dll4KO" = "italic(Dll4)^iDEC2w")
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

design <- c(patchwork::area(1, 1, 1, 3), patchwork::area(1, 4, 1, 4.5))

p24 <- Reduce( `+`,  p23)+patchwork::plot_layout(design = design)


cairo_pdf("Plots/Figure_7/Supplementary/Trials/Figure_S7_1i_UMAP.pdf",  width = 28, height = 6, family = "Arial")
p24
dev.off()


jpeg("Plots/Figure_7/Supplementary/Trials/Figure_S7_1i_UMAP.jpeg", width = 28, height = 6, units = 'in', res = 800)
p24
dev.off()

###################################################################


# Supplementary S7.1ii BarPlot

p1 <- dittoBarPlot(Liver, var = "FigClustering", group.by = "Condition", color.panel = my_palette_Rui_colors_B, xlab = NULL, main = NULL, x.reorder = c(1,3,2), y.breaks = c(0,0.25,0.5, 0.75, 1), min = 0, max = 1)+
  theme(axis.text=element_text(size=15), axis.title = element_text(size = 18), legend.text = element_text(size = 16), axis.text.x = element_text(angle = 25))+
  scale_x_discrete(labels=c("Control", bquote(italic(Dll4)^iDEC4d+Veh), bquote(italic(Dll4)^iDEC2w)))+
  geom_col(width = 0.1)+
  # scale_y_discrete(labels=c("0", "25", "50", "75", "100"))
  scale_y_continuous(labels = c("0", "25", "50", "75", "100"))
p1

cairo_pdf("Plots/Figure_7/Supplementary/Trials/Figure_S7_1ii_BarPlot.pdf",  width = 7, height = 6, family = "Arial")
p1
dev.off()


jpeg("Plots/Figure_7/Supplementary/Trials/Figure_S7_1ii_BarPlot.jpeg", width = 7, height = 6, units = 'in', res = 800)
p1
dev.off()

###################################################################


# Supplementary 7.2i Ctrl 4d, Dll4KO4d, Dll4KO2w

Liver <- readRDS("S:/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Liver_MacarenaFC/Figures_Paper/rds/Liver_20Jan21.Figures.Fig7Sup.2.rds")

levels(Liver@meta.data$Condition)

Liver <- subset(Liver, subset = Condition == "Control(4d)" | Condition == "Dll4KO(4d)+Vehicle" | Condition == "Dll4KO")

levels(Liver) <- c("Control(4d)", "Dll4KO(4d)+Vehicle", "Dll4KO")

Liver@meta.data$Condition <- Liver@active.ident

DimPlot(Liver)

p1 <- DimPlot(Liver, cols = my_palette_Rui_colors_B, pt.size = 1.2, split.by = "Condition", group.by = "FigClustering", combine = T)


new_labels <- c("Control(4d)" = "Control4d", "Dll4KO(4d)+Vehicle" = "italic(Dll4)^iDEC4d+Vehicle", "Dll4KO" = "italic(Dll4)^iDEC2w")
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

design <- c(patchwork::area(1, 1, 1, 3), patchwork::area(1, 4, 1, 4.5))

p24 <- Reduce( `+`,  p23)+patchwork::plot_layout(design = design)


cairo_pdf("Plots/Figure_7/Supplementary/Trials/Figure_S7_2i_UMAP.pdf",  width = 28, height = 6, family = "Arial")
p24
dev.off()


jpeg("Plots/Figure_7/Supplementary/Trials/Figure_S7_2i_UMAP.jpeg", width = 28, height = 6, units = 'in', res = 800)
p24
dev.off()

###################################################################


# Supplementary S7.2i BarPlot

p1 <- dittoBarPlot(Liver, var = "FigClustering", group.by = "Condition", color.panel = my_palette_Rui_colors_B, xlab = NULL, main = NULL, x.reorder = c(1,3,2), y.breaks = c(0,0.25,0.5, 0.75, 1), min = 0, max = 1)+
  theme(axis.text=element_text(size=15), axis.title = element_text(size = 18), legend.text = element_text(size = 16), axis.text.x = element_text(angle = 25))+
  scale_x_discrete(labels=c("Control4d", bquote(italic(Dll4)^iDEC4d+Veh), bquote(italic(Dll4)^iDEC2w)))+
  geom_col(width = 0.1)+
  # scale_y_discrete(labels=c("0", "25", "50", "75", "100"))
  scale_y_continuous(labels = c("0", "25", "50", "75", "100"))
p1

cairo_pdf("Plots/Figure_7/Supplementary/Trials/Figure_S7_2i_BarPlot.pdf",  width = 7, height = 6, family = "Arial")
p1
dev.off()


jpeg("Plots/Figure_7/Supplementary/Trials/Figure_S7_2i_BarPlot.jpeg", width = 7, height = 6, units = 'in', res = 800)
p1
dev.off()
