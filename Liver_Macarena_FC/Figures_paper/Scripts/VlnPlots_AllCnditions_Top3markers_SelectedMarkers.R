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

# rds

Liver <- readRDS("S:/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Liver_MacarenaFC/Figures_Paper/rds/Liver_20Jan21.Figures.All.Conditions.VlnPlots.rds")

#########################################################################

#  Order COnditions logically

table(Liver@meta.data$Condition)

levels(Liver) <- c("Control", "WT", "Control(4d)", "Dll4KO", "D4KO_aVEGF", "D1D4HET", "Dll4/MycLOF", 
                   "Dll4/MycKO", "Dll4KO(4d)+Vehicle", "Dll4KO(4d)+SL327", "NOTCH1KO", "N1N4HET", "RBPJKO", "Jag1/Jag2/Dll1KO", "WT_DBZ")

Liver@active.ident -> Liver@meta.data$Condition

#########################################################################

# Color palette

color_Vln <- c("firebrick1", "coral", "darkgoldenrod", "darkolivegreen", "chartreuse", "green", "palegreen", "mediumturquoise",
               "skyblue", "steelblue", "mediumblue", "darkorchid1", "orchid", "maroon1", "palevioletred")

#########################################################################

# Violin Plot Selected genes



genes <- c("DLL4", "DLL1", "JAG1", "NOTCH1", "RBPJ", "JAG2", "HES1", "HEY2", "HEY2", "HEYL", "ODC1", "KCNE3", "MKI67", 
                "WNT2", "MSR1", "LTBP4", "RSPO3", "GJA5", "GJA4", "CDKN1A", "STMN1", "VEGFA", "ESM1", "MYC", "CDKN1A", "CDKN1B", 
                "CDKN2A", "TRP53", "APLNR", "APLN", "NR2F2", "CD34", "EFNB2")

gene.names <- tolower(genes) %>% str_to_title()
genes <- as.data.frame(genes)
gene.names <- as.data.frame(gene.names)

# Testing 

VlnPlot(Liver, features = "CDH5", group.by = "Condition")+ theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.title.y=element_blank())+ggtitle("Hello")

p4 <- VlnPlot(Liver, features = "CDH5", group.by = "Condition", combine = F)

p4
#  Script

myplots <- vector('list', nrow(genes))


for (i in 1:nrow(genes)) {
  
  p21 <- VlnPlot(Liver, features = genes[i, 1], group.by = "Condition")+ NoLegend() + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), plot.title = element_text(face = "bold.italic"))+ggtitle(gene.names[i, 1])+
                  scale_x_discrete(labels=c("Control", "Control(wt)", "Control(4d)", expression(italic("Dll4"^"iDEC")), expression(italic(Dll4)^"iDEC"~"+aVEGF"), expression(italic(Dll4Het)^"iDEC"), expression(italic(Dll4/Myc)^"iDEC"~"New"), expression(italic(Dll4/Myc)^"iDEC"~"Old"),
                  expression(italic(Dll4)^"iDEC 4 days"~"+Veh"), expression(italic(Dll4)^"iDEC 4 days"~"+SL327"), expression(italic(Notch1)^"iDEC"), expression(italic(Notch1/2/4)^"iDEC"), expression(italic(Rbpj)^"iDEC"), expression(italic(Jag1/2)^"iDEC"), expression(DBZ)))
                  
  
  myplots[[i]] <- local({
    i <- i
    p21
  })
  
  # plegend <- VlnPlot(Liver, features = genes[i, 1], group.by = "Condition")
}



# plegend <- VlnPlot(Liver, features = "CDH5", group.by = "Condition") +
#   scale_fill_manual(values = color_Vln, labels=c("Control", "Control(wt)", "Control(4d)", expression(italic("Dll4"^"iDEC")), expression(italic(Dll4)^"iDEC"~"+aVEGF"), expression(italic(Dll4Het)^"iDEC"), expression(italic(Dll4/Myc)^"iDEC"~"New"), expression(italic(Dll4/Myc)^"iDEC"~"Old"),
#                              expression(italic(Dll4)^"iDEC 4 days"~"+Veh"), expression(italic(Dll4)^"iDEC 4 days"~"+SL327"), expression(italic(Notch1)^"iDEC"), expression(italic(Notch1/2/4)^"iDEC"), expression(italic(Rbpj)^"iDEC"), expression(italic(Jag1/2)^"iDEC"), expression(DBZ)))
# plegend
# 
# plegend <- VlnPlot(Liver, features = "CDH5", group.by = "Condition") + ggtitle(gene.names[3, 1])+
#   scale_x_discrete(labels=c("Control", "Control(wt)", "Control(4d)", expression(italic("Dll4"^"iDEC")), expression(italic(Dll4)^"iDEC"~"+aVEGF"), expression(italic(Dll4Het)^"iDEC"), expression(italic(Dll4/Myc)^"iDEC"~"New"), expression(italic(Dll4/Myc)^"iDEC"~"Old"),
#                                                  expression(italic(Dll4)^"iDEC 4 days"~"+Veh"), expression(italic(Dll4)^"iDEC 4 days"~"+SL327"), expression(italic(Notch1)^"iDEC"), expression(italic(Notch1/2/4)^"iDEC"), expression(italic(Rbpj)^"iDEC"), expression(italic(Jag1/2)^"iDEC"), expression(DBZ)))
# plegend
# 
# legend <- get_legend(plegend + theme(legend.position = "bottom", legend.justification = "center", legend.direction = "horizontal", legend.text = element_text(size = 18)))

p25 <- patchwork::wrap_plots(myplots, ncol = 6)

p25

cairo_pdf("Plots/Vln_plot_General/VlnPlot_Selected_Genes.pdf", width = 36, height = 36, family = "Arial")
p25
dev.off()


###############################################################################################################################


# Violin Plot top 3 Cluster genes

# Differential Analysis

Idents(Liver) <- "FigClustering"

levels(Liver) <- c("C0 - Unspecified quiescent capillaries", "C1a - Arterial capillaries", "C1v - Venous capillaries", "C2a - Large arteries", "C2v - Large veins",  "C3 - Activated capillaries", 
                    "C4 - Endothelial tip cells", "C5ip - Proliferating S-phase", "C5p - G2-M", "C6 - X")

DimPlot(Liver, label = T)

DEG_FigClustering <- FindAllMarkers(Liver, logfc.threshold = 0.5, test.use = "wilcox", only.pos = T, verbose = T, min.diff.pct = 0.1)
DEG_FigClustering$order <- 1:nrow(DEG_FigClustering)
write.xlsx(DEG_FigClustering, "Tables/DEG_FigClustering.wilcox.VlnPlots.AllConditions.xlsx", rowNames = T)
DEG_FigClustering.wilcox <- read_excel("Tables/DEG_FigClustering.wilcox.VlnPlots.AllConditions.xlsx")

Top3 <- DEG_FigClustering.wilcox %>% group_by(cluster) %>% top_n(n = 3, wt = -order)
write.xlsx(Top3, "Tables/DEG_FigClustering.wilcox.VlnPlots.AllConditions.Top3.xlsx", rowNames = T)


genes <- Top3$gene
gene.names <- tolower(genes) %>% str_to_title()

genes <- as.data.frame(genes)
gene.names <- as.data.frame(gene.names)

#  Script

myplots <- vector('list', nrow(genes))


for (i in 1:nrow(genes)) {
  
  p21 <- VlnPlot(Liver, features = genes[i, 1], group.by = "Condition")+ NoLegend() + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), plot.title = element_text(face = "bold.italic"))+ggtitle(gene.names[i, 1])+
    scale_x_discrete(labels=c("Control", "Control(wt)", "Control(4d)", expression(italic("Dll4"^"iDEC")), expression(italic(Dll4)^"iDEC"~"+aVEGF"), expression(italic(Dll4Het)^"iDEC"), expression(italic(Dll4/Myc)^"iDEC"~"New"), expression(italic(Dll4/Myc)^"iDEC"~"Old"),
                              expression(italic(Dll4)^"iDEC 4 days"~"+Veh"), expression(italic(Dll4)^"iDEC 4 days"~"+SL327"), expression(italic(Notch1)^"iDEC"), expression(italic(Notch1/2/4)^"iDEC"), expression(italic(Rbpj)^"iDEC"), expression(italic(Jag1/2)^"iDEC"), expression(DBZ)))
  
  
  myplots[[i]] <- local({
    i <- i
    p21
  })
  
  # plegend <- VlnPlot(Liver, features = genes[i, 1], group.by = "Condition")
}


p25 <- patchwork::wrap_plots(myplots, ncol = 6)

p25

cairo_pdf("Plots/Vln_plot_General/VlnPlot_Top3_Clusters.pdf", width = 36, height = 36, family = "Arial")
p25
dev.off()
