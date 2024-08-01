#  https://davemcg.github.io/post/lets-plot-scrna-dotplots/

library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(readxl)
library(gridExtra)
library(gridGraphics)
library(grid)
library(tidyverse)
library(ggpubr)
require(reshape2)

Liver <- readRDS("S:/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Liver_MacarenaFC/Figures_Paper/rds/Liver_20Jan21.Figures.All.Conditions.rds")
Liver <- subset(Liver, subset = Condition == "Control" | Condition == "Dll4KO" | Condition == "NOTCH1KO" | Condition == "RBPJKO")


Bestholtz_palette <- c("#DEDAD6", "#FEE392", "#FEC44E", "#FE9929", "#ED6F11", "#CC4C17", "#993411", "#65260C")
my_palette_Rui_colors_B <- c("#F4A753", "#E95A74", "#50B6EF", "#FC0808",  "#0E47D8", "#45FF8E", "#A80519", "#E28CF4", "#880088", "#C1B80C")


levels(Liver) <- c("Control", "Dll4KO", "NOTCH1KO", "RBPJKO")

Liver@active.ident -> Liver@meta.data$Condition

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

cairo_pdf("Plots/Figure_4/Figure_4D_DotPlot_Maca_genes.pdf", width = 22, height = 3, family = "Arial")
DotPlot(Liver, features = Dotplot_genes, col.min = 0, dot.scale = 5, scale = T)+
  theme(strip.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"),
        legend.title = element_text(size = 10), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm")) +
  scale_x_discrete(label = DotPlotGenesH)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  geom_vline(xintercept = c(10.5, 20.5, 30.5, 40.5, 50.5, 60.5, 70.5, 80.5, 90.5),linetype = 2 )

dev.off()

plot.new()
jpeg("Plots/Figure_4/Figure_4D_DotPlot_Maca_genes.jpeg", width = 22, height = 3, units = 'in', res = 800)
DotPlot(Liver, features = Dotplot_genes, col.min = 0, dot.scale = 5, scale = T)+
  theme(strip.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"),
        legend.title = element_text(size = 10), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm")) +
  scale_x_discrete(label = DotPlotGenesH)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  geom_vline(xintercept = c(10.5, 20.5, 30.5, 40.5, 50.5, 60.5, 70.5, 80.5, 90.5),linetype = 2 )+
  rect(100, 400, 125, 450, col = "green", border = "blue")
dev.off()


plot(c(0, 10), c(0, 1), type = "n", ylab="",yaxt="n", xlab = "", xaxt = "n", axes = FALSE)
i <- 1*(0:9)
## draw rectangles with bottom left (100, 300)+i  and top right (150, 380)+i
rect(0+i, 0, 1+i, 1, col=my_palette_Rui_colors_B)

grid.echo()
clusters_label <- grid.grab()
clusters_label

plot.new()

plot(c(0, 10), c(0, 1), type = "n", ylab="",yaxt="n", xlab = "", xaxt = "n", axes = FALSE)
blank <- grid.grab()

plot.new()

df <- data.frame(
  x = rep(c(2, 5, 7, 9, 12), 2),
  y = rep(c(1, 2), each = 5))

df <- data.frame(
  x = rep(c(5.5, 15.5, 25.5, 35.5, 45.5, 55.5, 65.5, 75.5, 85.5, 95.5)),
  y = rep(c(11), 1),
  z = factor(rep(1:10)),
  color = my_palette_Rui_colors)

jpeg("Plots/Figure_4/Figure_4D_DotPlot_Maca_genes.TRIAL.jpeg", width = 22, height = 3, units = 'in', res = 800)
DotPlot(Liver, features = Dotplot_genes, col.min = 0, dot.scale = 5, scale = T)+
  theme(strip.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"),
        legend.title = element_text(size = 10), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"), 
        plot.margin = unit(c(0, 0, 0, 0), "cm")
        ) +
  scale_x_discrete(label = DotPlotGenesH)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  geom_vline(xintercept = c(10.5, 20.5, 30.5, 40.5, 50.5, 60.5, 70.5, 80.5, 90.5),linetype = 2 )+
  # geom_rect(df, mapping=aes(xmin=x, xmax=x,fill=t), color="black", alpha=0.5)
  geom_raster(df, mapping = aes(x, y, fill = factor(z)), show.legend = F)+
  # scale_fill_brewer(type = "seq", palette = 1)
  # scale_fill_discrete(palette (my_palette_Rui_colors_B))
  scale_fill_manual(values = my_palette_Rui_colors_B)
  # scale_fill_hue(h = my_palette_Rui_colors_B)
  # scale_fill_gradientn(colours=c("#0000FFFF","#FFFFFFFF","#FF0000FF"))
dev.off()

dotplot
  i <- 10*(0:9)
## draw rectangles with bottom left (100, 300)+i  and top right (150, 380)+i
rect(0+i, 0, 1+i, 1, col=my_palette_Rui_colors_B)
ggdraw()

p27 <- plot_grid(blank, clusters_label, dotplot, ncol = 1, rel_heights = c(.02, 1), rel_widths = c(0.3, 1), axis = "l", hjust =)
p27

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
