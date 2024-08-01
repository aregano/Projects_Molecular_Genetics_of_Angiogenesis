library(Seurat)
library(ggplot2)
library("knitr")
library("rmarkdown")
library("yaml")
library("patchwork")
library("dittoSeq")
library("magick")
library("RColorBrewer")
library(cowplot)
library("stringr")
library("formattable")
library(readxl)
library(tidyverse)
library(showtext)
library(fgsea)

palette_Vln_Dll4 <- c("#606060", "#F94040", "#C19EDD")

Bestholtz_palette <- c("#DEDAD6", "#FEE392", "#FEC44E", "#FE9929", "#ED6F11", "#CC4C17", "#993411", "#65260C")

Liver <- readRDS("//Tierra/SC/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Liver_MacarenaFC/Figures_Paper/rds/Liver_20Jan21.Figures.All.Conditions.OldDll4Myc.rds")

table(Liver@meta.data$Condition)

Liver <- subset(Liver, subset = Condition == "Control" | Condition == "Dll4KO" | Condition == "Dll4/MycKO")

Idents(Liver) <- "Condition"

levels(Liver) <- c("Control", "Dll4KO", "Dll4/MycKO")

Liver@meta.data$Condition <- Liver@active.ident

genes_text <- c("Cxcr4", "Cxcl12", "Kcne3", "Esm1", "Apln")
genes <- genes_text %>% toupper()

genes_text <- as.data.frame(genes_text)
genes <- as.data.frame(genes)



#  Script

myplots <- vector('list', nrow(genes))


for (i in 1:nrow(genes)) {
  
  p21 <- VlnPlot(Liver, features = genes[i, 1], split.by = "Condition", cols = palette_Vln_Dll4) + NoLegend()+theme(axis.title.y = element_blank(), axis.text.x = element_blank(), plot.title = element_text(face = "bold.italic"))+ ggtitle(genes_text[i, 1])
  
  myplots[[i]] <- local({
    i <- i
    p21
  })
  
  plegend <- VlnPlot(Liver, features = genes[i, 1], split.by = "Condition", cols = palette_Vln_Dll4)+ ggtitle(genes_text[i, 1])
}

plegend

legend <- get_legend(plegend + theme(legend.position = "bottom", legend.justification = "center", text = element_text(size = 18)) +
  scale_fill_manual(values= palette_Vln_Dll4, labels=c("Control", expression(italic(Dll4)^"iDEC"), expression(italic(Dll4/Myc)^"iDEC"))))

legend

p25 <- patchwork::wrap_plots(myplots, ncol = 1)
p25
#myplots[1:nrow(genes)] + patchwork::plot_layout(byrow = T, widths = 25, heights = 5)
p26 <- plot_grid(p25, legend, ncol = 1, rel_heights = c(1, .05))

p26

x.axis <- ggdraw() + draw_label("Expression Level", fontface='bold', angle = 90, size = 24)

p27 <- plot_grid(x.axis, p26, ncol = 2, rel_widths = c(.05, 1))

p27

cairo_pdf("Plots/Rui_requests/VlnPlot_Liver_Tip_Cell_markers.pdf", width = 6, height = 15, family = "Arial")
p27
dev.off()

jpeg("Plots/Rui_requests/VlnPlot_Liver_Tip_Cell_markers.jpeg", width = 6, height = 15, units = 'in', res = 800)
p27
dev.off()

###################################################################################

# FeaturePlots same palette genes

myplots <- vector('list', nrow(genes))

for (i in 1:nrow(genes)) {
  
  p21 <- FeaturePlot(Liver, pt.size = 1, features = genes[i, 1], order = T, combine = F)
  
  p22 <- FeaturePlot(Liver, split.by = "Condition", order = T,features = genes[i, 1], combine = F)
  
  p21 <- lapply(X = p21, FUN = function(p) p + NoAxes()+ scale_colour_gradientn(colors = Bestholtz_palette)+ ggtitle(genes_text[i, 1])+ theme(plot.title = element_text(face = "bold.italic")))
  
  p22 <- lapply(X = p22, FUN = function(p) p + NoLegend() + NoAxes()+ scale_colour_gradientn(colors = Bestholtz_palette)+ theme(plot.title = element_blank()))
  
  p21 <- Reduce( `+`, p21 )+patchwork::plot_layout( ncol = 1 )
  
  p22 <- Reduce( `+`, p22 )+patchwork::plot_layout( ncol = 3 )
  
  p23 <- list(p21, p22)
  
  
  design <- c(patchwork::area(1, 1, 1, 1), patchwork::area(1, 2, 1, 4))
  
  p24 <- Reduce( `+`,  p23)+patchwork::plot_layout(design = design)
  
  
  myplots[[i]] <- local({
    i <- i
    p24
  })
  
  
}

p25 <- patchwork::wrap_plots(myplots, ncol = 1)

#myplots[1:nrow(genes)] + patchwork::plot_layout(byrow = T, widths = 25, heights = 5)
p25


cairo_pdf("Plots/Rui_requests/FeaturePlot.Liver.TipCell.pdf", width = 16, height = 20, family = "Arial")
p25
dev.off()

jpeg("Plots/Rui_requests/FeaturePlot.Liver.TipCell.jpeg", width = 16, height = 20, units = 'in', res = 800)
p25
dev.off()

