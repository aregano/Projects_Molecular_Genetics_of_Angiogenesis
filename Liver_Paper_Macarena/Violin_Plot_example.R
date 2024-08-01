# Libraries
library(Seurat)
library(ggplot2)
library(showtext)

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
  
  # plegend <- VlnPlot(Liver, features = genes[i, 1], group.by = "Condition")
}



plegend <- VlnPlot(Liver, features = genes[i, 1], group.by = "Condition") + scale_fill_manual(values= palette_Maca_Vln, labels=c("Control", expression(italic("Dll4"^"iDEC")), expression(italic(Notch1)^"iDEC"), expression(italic(Rbpj)^"iDEC")))
plegend
legend <- get_legend(plegend + theme(legend.position = "bottom", legend.justification = "center", legend.direction = "horizontal", legend.text = element_text(size = 18)))

p25 <- patchwork::wrap_plots(myplots, ncol = 1)
p25

p26 <- plot_grid(p25, legend, ncol = 1, rel_heights = c(1, .1))

p26

x.axis <- ggdraw() + draw_label("Expression Level", fontface='bold', angle = 90, size = 24)

p27 <- plot_grid(x.axis, p26, ncol = 2, rel_widths = c(.05, 1))

p27

setEPS()
showtext_begin()
postscript("Plots/Figure_4/Figure_4B_VlnPlot.eps", width = 8, height = 12, family = "Arial")
p27
dev.off()

jpeg("Plots/Figure_4/Figure_4B_VlnPlot.jpeg", width = 8, height = 12, units = 'in', res = 800)
p27
dev.off()
