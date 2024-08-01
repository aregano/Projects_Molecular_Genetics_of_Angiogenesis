# Inflammatory Response

library(Seurat)
library(msigdbr)
library(fgsea)
library(dplyr)
library(ggplot2)
library(presto)
library(tidyverse)
library(openxlsx)
library(SingleCellExperiment)
library(dittoSeq)
library(escape)
library(reshape2)
library(clusterProfiler)
library(pathview)
library(enrichplot)
library(rlist)
library(cowplot)
library(dplyr)
library(patchwork)

Bestholtz_palette <- c("#DEDAD6", "#FEE392", "#FEC44E", "#FE9929", "#ED6F11", "#CC4C17", "#993411", "#65260C")


###############################################################################################################################

#  Comparing ranks Inflammatory Response genes

Liver <- readRDS("//Tierra/SC/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Liver_MacarenaFC/Figures_Paper/rds/Liver_20Jan21.Figures.All.Conditions.OldDll4Myc.rds")
table(Liver@meta.data$Condition)

Liver <- subset(Liver, subset = Condition == "Control" | Condition == "Dll4KO" | Condition == "Dll4/MycKO")

Idents(Liver) <- Liver@meta.data$Condition

levels(Liver) <- c("Dll4/MycKO", "Dll4KO", "Control")

Liver@active.ident -> Liver@meta.data$Condition

table(Liver@meta.data$Condition)

ranks_Dll4Myc <- read.csv2("Tables/GSEA_Fig6/ranks.Dll4MycvsCtrl.csv", sep = ",")

ranks_Dll4 <- read.csv2("Tables/GSEA_Fig6/ranks.Dll4vsCtrl.csv", sep = ",")

msigdbr_species()

m_df<- msigdbr(species = "Mus musculus", category = 13)

h_df<- msigdbr(species = "Homo sapiens", category = 7)


table(m_df$gs_subcat)

table(h_gene_sets$gs_cat)

# Using Hallmark genes only

h_gene_sets = msigdbr(species = "mouse", category = "H")

head(m_df)

fgsea_sets<- h_gene_sets %>% split(x = .$gene_symbol, f = .$gs_name)

inflammatory_response_gene <- fgsea_sets$HALLMARK_INFLAMMATORY_RESPONSE

inflammatory_response_gene <- as.data.frame(inflammatory_response_gene)

library(tidyverse)

inflammatory_response_gene <- t(inflammatory_response_gene)

inflammatory_response_gene <- as.vector(inflammatory_response_gene)

inflammatory_rank_Dll4 <- ranks_Dll4 %>% select(one_of(inflammatory_response_gene))

inflammatory_rank_Dll4 <- ranks_Dll4[ ranks_Dll4$gene %in% inflammatory_response_gene, ]

inflammatory_rank_Dll4Myc <- ranks_Dll4Myc %>% select(one_of(inflammatory_response_gene))

inflammatory_rank_Dll4Myc <- ranks_Dll4Myc[ ranks_Dll4Myc$X %in% inflammatory_response_gene, ]

# Top 10 Hypoxia genes

Idents(Liver) <- "Condition"

Top10_Inflammatory_Dll4 <- inflammatory_rank_Dll4$gene[c(1:4, 6:11)] %>% toupper()

Top_10_Inflammatory_H_Dll4 <- inflammatory_rank_Dll4$gene[c(1:4, 6:11)]

Top10_Inflammatory_Dll4Myc <- inflammatory_rank_Dll4Myc$X[c(1:3, 6:8, 10, 13:15, 17)] %>% toupper()

Top_10_Inflammatory_H_Dll4Myc <- inflammatory_rank_Dll4Myc$X[c(1:3, 6:8, 10, 13:15, 17)]

Top20_Inflammatory <- c(Top10_Inflammatory_Dll4, Top10_Inflammatory_Dll4Myc)

Top20_Inflammatory <- sort(Top20_Inflammatory)

Top20_Inflammatory_H <- c(Top_10_Inflammatory_H_Dll4, Top_10_Inflammatory_H_Dll4Myc)

Top20_Inflammatory_H <- sort(Top20_Inflammatory_H)

levels(Liver) <- c("Dll4/MycKO", "Dll4KO", "Control")

Liver@active.ident -> Liver@meta.data$Condition

table(Liver@meta.data$Condition)

Top10_Inflammatory_Dll4

Top10_Inflammatory_Dll4Myc

Top20_Inflammatory

cairo_pdf("Plots/Inflammatory_Response/Figure_DotPlot_Maca_Inflammatory_Response_genes.Top20.Ctl.Dll4KO.Dll4Myc.pdf", width = 8, height = 2.2, family = "Arial")
DotPlot(Liver, features = Top20_Inflammatory, col.min = 0, dot.scale = 6, scale = T, scale.max = 80, scale.min = 0)+
  theme(strip.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"),
        legend.title = element_text(size = 10), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"), 
        plot.margin = unit(c(0, 0, 0, 0), "cm")
  ) +
  scale_y_discrete(labels=c(bquote(italic(Dll4/Myc)^iDEC), bquote(italic(Dll4)^iDEC), "Control"))+
  scale_x_discrete(label = Top20_Inflammatory_H)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  # geom_vline(xintercept = c(10.5, 21.5, 29.5, 39.5, 47.5),linetype = 2 )+
  theme(legend.box = "horizontal")
dev.off()



cairo_pdf("Plots/Inflammatory_Response/Figure_DotPlot_Maca_Inflammatory_genes.Dll4MycvCtl.Ctl.Dll4KO.Dll4Myc.pdf", width = 6, height = 2, family = "Arial")
DotPlot(Liver, features = Top10_Inflammatory_Dll4Myc, col.min = 0, dot.scale = 5, scale = T)+
  theme(strip.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"),
        legend.title = element_text(size = 10), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"), 
        plot.margin = unit(c(0, 0, 0, 0), "cm")
  ) +
  scale_y_discrete(labels=c(bquote(italic(Dll4/Myc)^iDEC), bquote(italic(Dll4)^iDEC), "Control"))+
  scale_x_discrete(label = Top_10_Inflammatory_H_Dll4Myc)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  # geom_vline(xintercept = c(10.5, 21.5, 29.5, 39.5, 47.5),linetype = 2 )+
  theme(legend.box = "horizontal")
dev.off()

cairo_pdf("Plots/Inflammatory_Response/Figure_DotPlot_Maca_Inflammatory_genes.Dll4MycvCtl.Ctl.Dll4KO.Dll4Myc.Unscaled.pdf", width = 6, height = 2, family = "Arial")
DotPlot(Liver, features = Top10_Inflammatory_Dll4, col.min = 0, dot.scale = 5, scale = F)+
  theme(strip.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"),
        legend.title = element_text(size = 10), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"), 
        plot.margin = unit(c(0, 0, 0, 0), "cm")
  ) +
  scale_y_discrete(labels=c(bquote(italic(Dll4/Myc)^iDEC), bquote(italic(Dll4)^iDEC), "Control"))+
  scale_x_discrete(label = Top_10_Inflammatory_H_Dll4)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  # geom_vline(xintercept = c(10.5, 21.5, 29.5, 39.5, 47.5),linetype = 2 )+
  theme(legend.box = "horizontal")
dev.off()
