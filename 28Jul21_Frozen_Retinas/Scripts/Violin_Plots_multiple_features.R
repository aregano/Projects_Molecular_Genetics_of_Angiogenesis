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

Liver <- subset(Liver, subset = Condition == "Control" | Condition == "Jag1/Jag2/Dll1KO(2w)")

Idents(Liver) <- Liver@meta.data$predicted.id

genes <- list(c("CD34", "CD36", "MSR1", "NTN4", "LTBP4", "EFNB2", "WNT2"))
genes <- as.data.frame(genes)


#  Script

myplots <- vector('list', nrow(genes))


for (i in 1:nrow(genes)) {
  
  p21 <- VlnPlot(Liver, features = genes[i, 1], idents = c("C0a", "C0v"), split.by = "Condition") + NoLegend()
  
  myplots[[i]] <- local({
    i <- i
    p21
  })
  
  plegend <- VlnPlot(Liver, features = genes[i, 1], idents = c("C0a", "C0v"), split.by = "Condition")
}

legend <- get_legend(plegend + theme(legend.position = "bottom", legend.justification = "center"))

p25 <- patchwork::wrap_plots(myplots, ncol = 4)
p25
#myplots[1:nrow(genes)] + patchwork::plot_layout(byrow = T, widths = 25, heights = 5)
p26 <- plot_grid(p25, legend, ncol = 1, rel_heights = c(1, .1))

p26

# p25 <- CombinePlots(plots = myplots, ncol = 1)
# p25

p21 <- VlnPlot(Liver, features = "CDH5", idents = c("C0a", "C0v"), split.by = "Condition")
p21
