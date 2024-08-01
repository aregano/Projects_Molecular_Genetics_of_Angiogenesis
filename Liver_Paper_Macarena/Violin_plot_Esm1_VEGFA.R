# Violin plot Esm1 and VegfA

library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(dplyr)
library(patchwork)

Liver <- readRDS("//Tierra/SC/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Liver_MacarenaFC/Figures_Paper/rds/Liver_20Jan21.Figures.All.Conditions.ppt_Rui.rds")

table(Liver@meta.data$Condition)

Liver <- subset(Liver, subset = Condition == "Control" | Condition == "Dll4KO" | Condition == "Dll4/MycKO" | Condition == "NOTCH1KO" | Condition == "RBPJKO" | Condition == "D4KO_aVEGF" | Condition == "Dll4KO(4d)+SL327")

Idents(Liver) <- "Condition"

levels(Liver) <- c("Control", "Dll4KO", "D4KO_aVEGF", "Dll4KO(4d)+SL327", "Dll4/MycKO", "NOTCH1KO", "RBPJKO")

Liver@active.ident -> Liver@meta.data$Condition

table(Liver@meta.data$Condition)

palette_Maca_Vln <- c("#606060", "#F94040", "#BB005E", "#CCF9E8", "mediumpurple1", "#05BE78", "#55A0FB")

# Figure 4N. Violin Plot

genes <- list(c("ESM1", "VEGFA"))
genes <- as.data.frame(genes)
gene.names <- list(c("Esm1", "Vegfa"))
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
  scale_fill_manual(values= palette_Maca_Vln, labels=c("Control", expression(italic("Dll4"^"iDEC")), expression(italic(Dl4)^"iDEC"+"aVEGF"), expression(italic(Dll4)^"iDEC"+"SL327"), expression(italic(Dll4/Myc)^"iDEC"), expression(italic(Notch1)^"iDEC"), expression(italic(Rbpj)^"iDEC")))+
  guides(colour = guide_legend(ncol = 7))

plegend
plegend+theme(legend.position = "bottom", legend.justification = "center", legend.key.size = unit(0.5, 'cm'), legend.direction = "horizontal", legend.text = element_text(size = 10))+guides(colour = guide_legend(nrow = 1))
plegend+theme(legend.position = "bottom", legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size = 10))+guides(fill = guide_legend(ncol = 7))
legend <- get_legend(plegend + theme(legend.position = "bottom",
                                     legend.justification = "center", 
                                     legend.key.size = unit(0.5, 'cm'), legend.direction = "horizontal", legend.text = element_text(size = 10))+
                       guides(fill = guide_legend(ncol = 7)))

p25 <- patchwork::wrap_plots(myplots, ncol = 1)

p26 <- plot_grid(p25, legend, ncol = 1, rel_heights = c(1, .07))+guides(colour = guide_legend(nrow = 1))

x.axis <- ggdraw() + draw_label("Expression Level", fontface='bold', angle = 90, size = 24)

p27 <- plot_grid(x.axis, p26, ncol = 2, rel_widths = c(.1, 1))

p27

cairo_pdf("Plots/Hypoxia/Figure_Esm1_VegfA_VlnPlot.pdf", width = 9, height = 6, family = "Arial")
p27
dev.off()

jpeg("Plots/Hypoxia/Figure_Esm1_VegfA_VlnPlot.jpeg", width = 9, height = 6, units = 'in', res = 800)
p27
dev.off()
