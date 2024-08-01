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

Liver <- readRDS("//Tierra/SC/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Liver_MacarenaFC/Figures_Paper/rds/Liver_20Jan21.Figures.All.Conditions.ppt_Rui.rds")

table(Liver@meta.data$Condition)

Liver <- subset(Liver, subset = Condition == "Control" | Condition == "Dll4KO" | Condition == "Dll4/MycKO" | Condition == "NOTCH1KO" | Condition == "RBPJKO" | Condition == "D4KO_aVEGF" | Condition == "Dll4KO(4d)+SL327")

Idents(Liver) <- "Condition"

levels(Liver) <- c("Control", "Dll4KO", "D4KO_aVEGF", "Dll4KO(4d)+SL327", "Dll4/MycKO", "NOTCH1KO", "RBPJKO")

Liver@active.ident -> Liver@meta.data$Condition

table(Liver@meta.data$Condition)

Liver4d <- readRDS("S:/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Liver_MacarenaFC/Figures_Paper/rds/Liver4d.Mapping.rds")


palette_Maca_Vln <- c("#606060", "#F94040", "#BB005E", "#CCF9E8", "mediumpurple1", "#05BE78", "#55A0FB")

palette_Maca_GSEA <- c("#F94040", "#05BE78", "#55A0FB")

BuRd <- c("#052F60", "#2066AC", "#4392C3", "#92C5DE", "#D1E4EF", "#F7F7F7", "#FCDBC6", "#F4A482", "#D55F4D", "#B1172B", "#67001F")


#####################################################################################################################

# Figure Violin Plot

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

######################################################################################

#  VlnPlot SL327

palette_Maca_Vln <- c("#606060", "#F94040", "#BB005E")

genes <- Top10_Hypoxia
genes <- as.data.frame(genes)
gene.names <- Top10_Hypoxia_H
gene.names <- as.data.frame(gene.names)

#  Script

myplots <- vector('list', nrow(genes))


for (i in 1:nrow(genes)) {
  
  p21 <- VlnPlot(Liver4d, features = genes[i, 1], group.by = "Condition", cols = palette_Maca_Vln)+ NoLegend() + theme(text = element_text(family="Arial"), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.title.y=element_blank(), plot.title = element_text(face = "bold.italic"))+ggtitle(gene.names[i, 1])
  
  myplots[[i]] <- local({
    i <- i
    p21
  })
  
}

plegend <- VlnPlot(Liver4d, features = genes[i, 1], group.by = "Condition") +
  scale_fill_manual(values= palette_Maca_Vln, labels=c("Control4d", expression(italic(Dll4)^"iDEC4d"), expression(italic(Dll4)^"iDEC4d"+"SL327")))
plegend
legend <- get_legend(plegend + theme(legend.position = "bottom", legend.justification = "center", legend.direction = "horizontal", legend.text = element_text(size = 14)))

p25 <- patchwork::wrap_plots(myplots, ncol = 1)

p26 <- plot_grid(p25, legend, ncol = 1, rel_heights = c(1, .07))

x.axis <- ggdraw() + draw_label("Expression Level", fontface='bold', angle = 90, size = 24)

p27 <- plot_grid(x.axis, p26, ncol = 2, rel_widths = c(.1, 1))

p27

cairo_pdf("Plots/Hypoxia/Figure_VlnPlot_Hypoxia_genes_SL327.pdf", width = 5, height = 33, family = "Arial")
p27
dev.off()

###############################################################################################################################
#  Comparing ranks Hypoxia genes

Liver <- readRDS("//Tierra/SC/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Liver_MacarenaFC/Figures_Paper/rds/Liver_20Jan21.Figures.All.Conditions.OldDll4Myc.rds")

table(Liver@meta.data$Condition)

Liver <- subset(Liver, subset = Condition == "Control" | Condition == "Dll4KO" | Condition == "Dll4/MycKO" | Condition == "D4KO_aVEGF")

Idents(Liver) <- "Condition"

levels(Liver) <- c("Control", "Dll4KO", "Dll4/MycKO", "D4KO_aVEGF")

Liver@active.ident -> Liver@meta.data$Condition

table(Liver@meta.data$Condition)

ranks <- read.csv2("Tables/GSEA_Fig6/ranks.Dll4vsCtrl.csv", sep = ",")

hypoxia_genes <- fgsea_sets$HALLMARK_HYPOXIA

hypoxia_genes <- as.data.frame(hypoxia_genes)

library(tidyverse)

hypoxia_genes <- t(hypoxia_genes)

hypoxia_genes <- as.vector(hypoxia_genes)

hypoxia_rank <- ranks[ ranks$gene %in% hypoxia_genes, ]

# Top 10 Hypoxia genes

Idents(Liver) <- "Condition"

levels(Liver) <- c("D4KO_aVEGF", "Dll4/MycKO", "Dll4KO", "Control")

Top10_Hypoxia <- hypoxia_rank$gene[c(1:10, 32)] %>% toupper()

Top10_Hypoxia_H <- hypoxia_rank$gene[c(1:10, 32)]

cairo_pdf("Plots/Hypoxia/Figure_DotPlot_Maca_Hypoxia_genes.Final.pdf", width = 7, height = 2, family = "Arial")
DotPlot(Liver, features = Top10_Hypoxia, col.min = 0, dot.scale = 5, scale = T)+
  theme(strip.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"),
        legend.title = element_text(size = 10), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"), 
        plot.margin = unit(c(0, 0, 0, 0), "cm")
  ) +
  scale_y_discrete(labels=c(bquote(italic(Dll4)^iDEC+aVEGF), bquote(italic(Dll4/Myc)^iDEC), bquote(italic(Dll4)^iDEC), "Control"))+
  scale_x_discrete(label = Top10_Hypoxia_H)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
# geom_vline(xintercept = c(10.5, 21.5, 29.5, 39.5, 47.5),linetype = 2 )+
  theme(legend.box = "horizontal")
dev.off()
#

###############################################################################################################################
#  Comparing ranks Hypoxia genes

Liver <- readRDS("//Tierra/SC/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Liver_MacarenaFC/Figures_Paper/rds/Liver_20Jan21.Figures.All.Conditions.OldDll4Myc.rds")

table(Liver@meta.data$Condition)

Liver <- subset(Liver, subset = Condition == "Control" | Condition == "Dll4KO" | Condition == "Dll4/MycKO" | Condition == "D4KO_aVEGF")

Idents(Liver) <- "Condition"

levels(Liver) <- c("Control", "Dll4KO", "Dll4/MycKO", "D4KO_aVEGF")

Liver@active.ident -> Liver@meta.data$Condition

table(Liver@meta.data$Condition)

ranks <- read.csv2("Tables/GSEA_Fig6/ranks.Dll4MycvsCtrl.csv", sep = ",")

hypoxia_genes <- fgsea_sets$HALLMARK_HYPOXIA

hypoxia_genes <- as.data.frame(hypoxia_genes)

library(tidyverse)

hypoxia_genes <- t(hypoxia_genes)

hypoxia_genes <- as.vector(hypoxia_genes)

hypoxia_rank <- ranks[ ranks$X %in% hypoxia_genes, ]

# Top 10 Hypoxia genes

Idents(Liver) <- "Condition"

levels(Liver) <- c("D4KO_aVEGF", "Dll4/MycKO", "Dll4KO", "Control")

Top10_Hypoxia <- hypoxia_rank$X[c(1:10, 148)] %>% toupper()

Top10_Hypoxia_H <- hypoxia_rank$X[c(1:10, 148)]

Top10_Hypoxia

cairo_pdf("Plots/Hypoxia/Figure_DotPlot_Maca_Hypoxia_genes.Dll4MycvCtl.Final.pdf", width = 7, height = 2, family = "Arial")
DotPlot(Liver, features = Top10_Hypoxia, col.min = 0, dot.scale = 5, scale = T)+
  theme(strip.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"),
        legend.title = element_text(size = 10), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"), 
        plot.margin = unit(c(0, 0, 0, 0), "cm")
  ) +
  scale_y_discrete(labels=c(bquote(italic(Dll4)^iDEC+aVEGF), bquote(italic(Dll4/Myc)^iDEC), bquote(italic(Dll4)^iDEC), "Control"))+
  scale_x_discrete(label = Top10_Hypoxia_H)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  # geom_vline(xintercept = c(10.5, 21.5, 29.5, 39.5, 47.5),linetype = 2 )+
  theme(legend.box = "horizontal")
dev.off()
#


############################################################


# Heatmap Hallmark_Hypoxia genes

msigdbr_species()

m_df<- msigdbr(species = "Mus musculus", category = 13)

h_df<- msigdbr(species = "Homo sapiens", category = 7)


table(m_df$gs_subcat)

table(h_gene_sets$gs_cat)

# Using Hallmark genes only

h_gene_sets = msigdbr(species = "mouse", category = "H")

head(m_df)

fgsea_sets<- h_gene_sets %>% split(x = .$gene_symbol, f = .$gs_name)

hypoxia_genes <- fgsea_sets$HALLMARK_HYPOXIA


Liver <- readRDS("S:/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Liver_MacarenaFC/Figures_Paper/rds/Liver_20Jan21.Figures.All.Conditions.OldDll4Myc.rds")

Idents(Liver) <- "Condition"

table(Liver@meta.data$Condition)

# Standard Ctrl, Dll4KO, Notch1 and Rbpj

Liver <- subset(Liver, subset = Condition == "Control" | Condition == "Dll4KO" | Condition == "NOTCH1KO" | Condition == "RBPJKO")

levels(Liver) <- c("Control", "Dll4KO", "NOTCH1KO", "RBPJKO")

palette_Maca_Vln <- c("#606060", "#F94040", "#05BE78", "#55A0FB")

Liver@active.ident -> Liver@meta.data$Condition

avgexp <- AverageExpression(Liver, return.seurat = T)

avgexp@meta.data$Condition <- avgexp@active.ident

levels(avgexp) <- c("Control", "Dll4KO", "NOTCH1KO", "RBPJKO")

# Ctrl, Dll4KO, Dll4/MycKO, Dll4KO+antiVEGFA

Liver <- subset(Liver, subset = Condition == "Control" | Condition == "Dll4KO" | Condition == "Dll4/MycKO" | Condition == "D4KO_aVEGF")

levels(Liver) <- c("Control", "Dll4KO", "Dll4/MycKO", "D4KO_aVEGF")

palette_Maca_Vln <- c("#606060", "#F94040", "#C19EDD", "#BB005E")

Liver@active.ident -> Liver@meta.data$Condition

avgexp <- AverageExpression(Liver, return.seurat = T)

avgexp@meta.data$Condition <- avgexp@active.ident

levels(avgexp) <- c("Control", "Dll4KO", "Dll4/MycKO", "D4KO_aVEGF")


# Liver 4 days: Ctrol4d, Dll4KO4days, Dll4+SL327


Liver4d <- readRDS("S:/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Liver_MacarenaFC/Figures_Paper/rds/Liver4d.Mapping.rds")

table(Liver4d@meta.data$Condition)

Idents(Liver4d) <- "Condition"

levels(Liver4d) <- c("Control(4d)", "Dll4KO(4d)+Vehicle", "Dll4KO(4d)+SL327")

palette_Maca_Vln <- c("#606060", "#F94040", "#CCF9E8")

Liver4d@active.ident -> Liver4d@meta.data$Condition

avgexp <- AverageExpression(Liver4d, return.seurat = T)

avgexp@meta.data$Condition <- avgexp@active.ident

levels(avgexp) <- c("Control(4d)", "Dll4KO(4d)+Vehicle", "Dll4KO(4d)+SL327")


# Same for all conditions


ranks <- read.csv2("Tables/GSEA_Fig6/ranks.Dll4vsCtrl.csv", sep = ",")

hypoxia_genes <- fgsea_sets$HALLMARK_HYPOXIA

hypoxia_genes <- as.data.frame(hypoxia_genes)

library(tidyverse)

hypoxia_genes <- t(hypoxia_genes)

hypoxia_genes <- as.vector(hypoxia_genes)

hypoxia_rank <- ranks[ ranks$gene %in% hypoxia_genes, ]

hypoxia <- hypoxia_rank[1]

Hypoxia_genes <- as.character(hypoxia)

Hypoxia_genes <- parse(text = Hypoxia_genes)
Hypoxia_genes <- eval(Hypoxia_genes)

# Option 2 , only 173 genes but they do work

Hypoxia_genes <- read.csv2("Tables/Hypoxia_Hallmark_genes.csv", sep = ",", header = F)

Hypoxia_genes <- as.character(Hypoxia_genes)


#  Same in both

Hypoxia_genes <- toupper(Hypoxia_genes)



Hypoxia_genes

# Take it from the fgsea dataset in mgsidbr

Idents(Liver) <- "Condition"






hypoxia_genes




p1 <- dittoHeatmap(avgexp, genes = Hypoxia_genes, 
                   annot.by = "Condition",
                   annot.colors = palette_Maca_Vln,
                   fontsize = 7,
                   heatmap.colors = BuRd,
                   cluster_cols = F,
                   cluster_rows = F,
                   scaled.to.max = F,
                   scale = "none",
                   complex = T, data.out = F)


p1

cairo_pdf("Plots/Hypoxia/GSEA_Hypoxia_Genes_Heatmap_unscaled.pdf",  width = 4, height = 30, family = "Arial")
p1
dev.off()

p1 <- dittoHeatmap(avgexp, genes = Hypoxia_genes, 
                   annot.by = "Condition",
                   annot.colors = palette_Maca_Vln,
                   fontsize = 7,
                   heatmap.colors = BuRd,
                   cluster_cols = F,
                   cluster_rows = F,
                   scaled.to.max = F,
                   complex = T, data.out = F)


cairo_pdf("Plots/Hypoxia/GSEA_Hypoxia_Genes_Heatmap_scaled.Ctrl_Dll4KO(4d)_Dll4KO(4d)+SL327.pdf",  width = 3.5, height = 30, family = "Arial")
p1
dev.off()


###############################################################################################################################

#  Comparing ranks Hypoxia genes

Liver <- readRDS("//Tierra/SC/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Liver_MacarenaFC/Figures_Paper/rds/Liver_20Jan21.Figures.All.Conditions.OldDll4Myc.rds")
table(Liver@meta.data$Condition)

Liver <- subset(Liver, subset = Condition == "Control" | Condition == "Dll4KO" | Condition == "Dll4/MycKO")

Idents(Liver) <- Liver@meta.data$Condition

levels(Liver) <- c("Control", "Dll4KO", "Dll4/MycKO")

Liver@active.ident -> Liver@meta.data$Condition

table(Liver@meta.data$Condition)

# Dll4 v Ctl

ranks <- read.csv2("Tables/GSEA_Fig6/ranks.Dll4vsCtrl.csv", sep = ",")

# Dll4/Myc v Ctl

ranks <- read.csv2("Tables/GSEA_Fig6/ranks.Dll4MycvsCtrl.csv", sep = ",")

hypoxia_genes <- fgsea_sets$HALLMARK_HYPOXIA

hypoxia_genes <- as.data.frame(hypoxia_genes)

library(tidyverse)

hypoxia_genes <- t(hypoxia_genes)

hypoxia_genes <- as.vector(hypoxia_genes)

hypoxia_rank <- ranks %>% select(one_of(hypoxia_genes))

hypoxia_rank_Dll4 <- ranks[ ranks$gene %in% hypoxia_genes, ]

hypoxia_rank_Dll4Myc <- ranks[ ranks$X %in% hypoxia_genes, ]

# Top 10 Hypoxia genes

Idents(Liver) <- "Condition"

Top10_Hypoxia_Dll4 <- hypoxia_rank_Dll4$gene[c(1:10)] %>% toupper()

Top10_Hypoxia_H_Dll4 <- hypoxia_rank_Dll4$gene[c(1:10)]

Top10_Hypoxia_Dll4Myc <- hypoxia_rank_Dll4Myc$X[c(1:6, 8, 10:11, 14)] %>% toupper()

Top10_Hypoxia_H_Dll4Myc <- hypoxia_rank_Dll4Myc$X[c(1:6, 8, 10:11, 14)]

Top20_Hypoxia <- c(Top10_Hypoxia_Dll4, Top10_Hypoxia_Dll4Myc)

Top20_Hypoxia <- sort(Top20_Hypoxia)

Top20_Hypoxia_H <- c(Top10_Hypoxia_H_Dll4, Top10_Hypoxia_H_Dll4Myc)

Top20_Hypoxia_H <- sort(Top20_Hypoxia_H)

levels(Liver) <- c("Dll4/MycKO", "Dll4KO", "Control")

Liver@active.ident -> Liver@meta.data$Condition

table(Liver@meta.data$Condition)

cairo_pdf("Plots/Hypoxia/Figure_DotPlot_Maca_Hypoxia_genes.Top20.Ctl.Dll4KO.Dll4Myc.pdf", width = 8, height = 2.2, family = "Arial")
DotPlot(Liver, features = Top20_Hypoxia, col.min = 0, dot.scale = 6, scale = T, scale.max = 60, scale.min = 0)+
  theme(strip.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"),
        legend.title = element_text(size = 10), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"), 
        plot.margin = unit(c(0, 0, 0, 0), "cm")
  ) +
  scale_y_discrete(labels=c(bquote(italic(Dll4/Myc)^iDEC), bquote(italic(Dll4)^iDEC), "Control"))+
  scale_x_discrete(label = Top20_Hypoxia_H)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  # geom_vline(xintercept = c(10.5, 21.5, 29.5, 39.5, 47.5),linetype = 2 )+
  theme(legend.box = "horizontal")
dev.off()



cairo_pdf("Plots/Hypoxia/Figure_DotPlot_Maca_Hypoxia_genesDll4vCtl.Ctl.Dll4KO.Dll4Myc.pdf", width = 6, height = 2, family = "Arial")
DotPlot(Liver, features = Top10_Hypoxia, col.min = 0, dot.scale = 5, scale = T)+
  theme(strip.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"),
        legend.title = element_text(size = 10), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"), 
        plot.margin = unit(c(0, 0, 0, 0), "cm")
  ) +
  scale_y_discrete(labels=c(bquote(italic(Dll4/Myc)^iDEC), bquote(italic(Dll4)^iDEC), "Control"))+
  scale_x_discrete(label = Top10_Hypoxia_H)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
# geom_vline(xintercept = c(10.5, 21.5, 29.5, 39.5, 47.5),linetype = 2 )+
  theme(legend.box = "horizontal")
dev.off()

cairo_pdf("Plots/Hypoxia/Figure_DotPlot_Maca_Hypoxia_genes.Dll4MycvCtl.Ctl.Dll4KO.Dll4Myc.pdf", width = 6, height = 2, family = "Arial")
DotPlot(Liver, features = Top10_Hypoxia, col.min = 0, dot.scale = 5, scale = T)+
  theme(strip.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"),
        legend.title = element_text(size = 10), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"), 
        plot.margin = unit(c(0, 0, 0, 0), "cm")
  ) +
  scale_y_discrete(labels=c(bquote(italic(Dll4/Myc)^iDEC), bquote(italic(Dll4)^iDEC), "Control"))+
  scale_x_discrete(label = Top10_Hypoxia_H)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  # geom_vline(xintercept = c(10.5, 21.5, 29.5, 39.5, 47.5),linetype = 2 )+
  theme(legend.box = "horizontal")
dev.off()
###############################################################################################################################

#  Comparing ranks Hypoxia genes

Liver <- readRDS("//Tierra/SC/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Liver_MacarenaFC/Figures_Paper/rds/Liver_20Jan21.Figures.All.Conditions.OldDll4Myc.rds")
table(Liver@meta.data$Condition)

Liver <- subset(Liver, subset = Condition == "Control" | Condition == "Dll4KO" | Condition == "Dll4/MycKO")

Idents(Liver) <- Liver@meta.data$Condition

levels(Liver) <- c("Dll4/MycKO", "Dll4KO", "Control")

Liver@active.ident -> Liver@meta.data$Condition

table(Liver@meta.data$Condition)

ranks <- read.csv2("Tables/GSEA_Fig6/ranks.Dll4MycvsCtrl.csv", sep = ",")

hypoxia_genes <- fgsea_sets$HALLMARK_HYPOXIA

hypoxia_genes <- as.data.frame(hypoxia_genes)

library(tidyverse)

hypoxia_genes <- t(hypoxia_genes)

hypoxia_genes <- as.vector(hypoxia_genes)

hypoxia_rank <- ranks %>% select(one_of(hypoxia_genes))

hypoxia_rank <- ranks[ ranks$X %in% hypoxia_genes, ]

# Top 10 Hypoxia genes

Idents(Liver) <- "Condition"

Top10_Hypoxia <- hypoxia_rank$X[c(1:10)] %>% toupper()

Top10_Hypoxia_H <- hypoxia_rank$X[c(1:10)]

Top10_Hypoxia

cairo_pdf("Plots/Hypoxia/Figure_DotPlot_Maca_Hypoxia_genes.Dll4MycvCtl.Ctl.Dll4KO.Dll4Myc.pdf", width = 6, height = 2, family = "Arial")
DotPlot(Liver, features = Top10_Hypoxia, col.min = 0, dot.scale = 5, scale = T)+
  theme(strip.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"),
        legend.title = element_text(size = 10), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"), 
        plot.margin = unit(c(0, 0, 0, 0), "cm")
  ) +
  scale_y_discrete(labels=c(bquote(italic(Dll4/Myc)^iDEC), bquote(italic(Dll4)^iDEC), "Control"))+
  scale_x_discrete(label = Top10_Hypoxia_H)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  # geom_vline(xintercept = c(10.5, 21.5, 29.5, 39.5, 47.5),linetype = 2 )+
  theme(legend.box = "horizontal")
dev.off()


###############################################################################################################################

# Figure S4.4C Violin Plot

genes <- list(c("KCNE3", "ESM1", "APLN", "ANGPT2", "CD34", "MYC", "ODC1", "VEGFA"))
genes <- as.data.frame(genes)
gene.names <- list(c("Kcne3", "Esm1", "Apln", "Angpt2", "Cd34", "Myc", "Odc1", "Vegfa"))
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
  scale_fill_manual(values= palette_Maca_Vln, labels=c("Control", expression(italic("Dll4"^"iDEC")), expression(italic(Notch1)^"iDEC"), expression(italic(Rbpj)^"iDEC")))
plegend
legend <- get_legend(plegend + theme(legend.position = "bottom", legend.justification = "center", legend.direction = "horizontal", legend.text = element_text(size = 18)))

p25 <- patchwork::wrap_plots(myplots, ncol = 1)

p25

p26 <- plot_grid(p25, legend, ncol = 1, rel_heights = c(1, .03))

p26

x.axis <- ggdraw() + draw_label("Expression Level", fontface='bold', angle = 90, size = 48)

p27 <- plot_grid(x.axis, p26, ncol = 2, rel_widths = c(.1, 1))

p27

cairo_pdf("Plots/Figure_4/Suplementary/Figure_S4_4C_VlnPlot.pdf", width = 9, height = 30, family = "Arial")
p27
dev.off()

jpeg("Plots/Figure_4/Suplementary/Figure_S4_4C_VlnPlot.jpeg", width = 8, height = 30, units = 'in', res = 800)
p27
dev.off()

########################################################################################
Liver <- readRDS("S:/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Liver_MacarenaFC/Figures_Paper/rds/Liver_20Jan21.Figures.All.Conditions.OldDll4Myc.rds")
Liver <- subset(Liver, subset = Condition == "Control" | Condition == "Dll4KO" | Condition == "NOTCH1KO" | Condition == "RBPJKO")

palette_Maca_Vln <- c("#606060", "#F94040", "#05BE78", "#55A0FB")

levels(Liver) <- c("Control", "Dll4KO", "NOTCH1KO", "RBPJKO")

Liver@active.ident -> Liver@meta.data$Condition

# Trying Hypoxia genes

genes <- c("AK4", "P4HA2", "ADM", "VEGFA", "PFKL", "ERO1L", "PFKP", "P4HA1", "PGM2", "HK1", "SELENBP1", "CCNG2", "GYS1", "XPNPEP1", "KDM3A", "SLC2A1")

gene.names <- genes %>% str_to_title()

genes <- as.data.frame(genes)

gene.names <- as.data.frame(gene.names)

myplots <- vector('list', nrow(genes))


for (i in 1:nrow(genes)) {
  
  p21 <- VlnPlot(Liver, features = genes[i, 1], group.by = "Condition", cols = palette_Maca_Vln)+ NoLegend() + theme(text = element_text(family="Arial"), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.title.y=element_blank(), plot.title = element_text(face = "bold.italic"))+ggtitle(gene.names[i, 1])
  
  myplots[[i]] <- local({
    i <- i
    p21
  })
  
}

plegend <- VlnPlot(Liver, features = genes[i, 1], group.by = "Condition") +
  scale_fill_manual(values= palette_Maca_Vln, labels=c("Control", expression(italic("Dll4"^"iDEC")), expression(italic(Notch1)^"iDEC"), expression(italic(Rbpj)^"iDEC")))
plegend
legend <- get_legend(plegend + theme(legend.position = "bottom", legend.justification = "center", legend.direction = "horizontal", legend.text = element_text(size = 18)))

p25 <- patchwork::wrap_plots(myplots, ncol = 1)

p25

p26 <- plot_grid(p25, legend, ncol = 1, rel_heights = c(1, .03))

p26

x.axis <- ggdraw() + draw_label("Expression Level", fontface='bold', angle = 90, size = 48)

p27 <- plot_grid(x.axis, p26, ncol = 2, rel_widths = c(.1, 1))

p27

cairo_pdf("Plots/Hypoxia/Hypoxia_VlnPlot.Ctl.Dll4.Rbpj.Notch1.pdf", width = 9, height = 56, family = "Arial")
p27
dev.off()

jpeg("Plots/Hypoxia/Hypoxia_VlnPlot.Ctl.Dll4.Rbpj.Notch1.jpeg", width = 8, height = 56, units = 'in', res = 800)
p27
dev.off()

##############################################

Liver <- readRDS("S:/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Liver_MacarenaFC/Figures_Paper/rds/Liver_20Jan21.Figures.All.Conditions.ppt_Rui.rds")

table(Liver@meta.data$Condition)

Liver <- subset(Liver, subset = Condition == "Control" | Condition == "Dll4KO" | Condition == "Dll4/MycKO" | Condition == "D4KO_aVEGF" | Condition == "Dll4KO(4d)+SL327")

palette_Maca_Vln <- c("#606060", "#F94040", "#C19EDD", "#BB005E", "#CCF9E8")

levels(Liver) <- c("Control", "Dll4KO", "Dll4/MycKO", "D4KO_aVEGF", "Dll4KO(4d)+SL327")

Liver@active.ident -> Liver@meta.data$Condition

# Trying Hypoxia genes

genes <- c("AK4", "P4HA2", "ADM", "VEGFA", "PFKL", "ERO1L", "PFKP", "P4HA1", "PGM2", "HK1", "SELENBP1", "CCNG2", "GYS1", "XPNPEP1", "KDM3A", "SLC2A1")

gene.names <- genes %>% str_to_title()

genes <- as.data.frame(genes)

gene.names <- as.data.frame(gene.names)

myplots <- vector('list', nrow(genes))


for (i in 1:nrow(genes)) {
  
  p21 <- VlnPlot(Liver, features = genes[i, 1], group.by = "Condition", cols = palette_Maca_Vln)+ NoLegend() + theme(text = element_text(family="Arial"), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.title.y=element_blank(), plot.title = element_text(face = "bold.italic"))+ggtitle(gene.names[i, 1])
  
  myplots[[i]] <- local({
    i <- i
    p21
  })
  
}

plegend <- VlnPlot(Liver, features = genes[i, 1], group.by = "Condition") +
  scale_fill_manual(values= palette_Maca_Vln, labels=c("Control", expression(italic("Dll4")^"iDEC"), expression(italic(Dll4/Myc)^"iDEC"), expression(italic(Dll4)^"iDEC"+"aVEGF"), expression(italic(Dll4)^"iDEC"+"SL327")))
plegend
legend <- get_legend(plegend + theme(legend.position = "bottom", legend.justification = "center", legend.direction = "horizontal", legend.text = element_text(size = 14)))

p25 <- patchwork::wrap_plots(myplots, ncol = 1)

p25

p26 <- plot_grid(p25, legend, ncol = 1, rel_heights = c(1, .03))

p26

x.axis <- ggdraw() + draw_label("Expression Level", fontface='bold', angle = 90, size = 48)

p27 <- plot_grid(x.axis, p26, ncol = 2, rel_widths = c(.1, 1))

p27

cairo_pdf("Plots/Hypoxia/Hypoxia_VlnPlot.Ctl.Dll4.Dll4_Myc.Dll4+aVEGF.Dll4+SL327.pdf", width = 9, height = 56, family = "Arial")
p27
dev.off()

jpeg("Plots/Hypoxia/Hypoxia_VlnPlot.Ctl.Dll4.Dll4_Myc.Dll4+aVEGF.Dll4+SL327.jpeg", width = 8, height = 56, units = 'in', res = 800)
p27
dev.off()


#############################################


# Feature Plot

Bestholtz_palette <- c("#DEDAD6", "#FEE392", "#FEC44E", "#FE9929", "#ED6F11", "#CC4C17", "#993411", "#65260C")




genes <- as.data.frame(genes)

#  Script

myplots <- vector('list', nrow(genes))


for (i in 1:nrow(genes)) {
  
  p21 <- FeaturePlot(Liver, pt.size = 1.5, features = genes[i, 1], order = T, combine = F, slot = "data",
                     # cols = Bestholtz_palette
  )
  
  p22 <- FeaturePlot(Liver, split.by = "Condition", order = T,features = genes[i, 1], combine = F, slot = "data",
                     # cols = Bestholtz_palette
  )
  
  p21 <- lapply(X = p21, FUN = function(p) p + 
                  scale_colour_gradientn(colors =  Bestholtz_palette) + 
                  NoAxes())
  
  p22 <- lapply(X = p22, FUN = function(p) p + NoLegend() +
                  scale_colour_gradientn(colors = Bestholtz_palette) +
                  NoAxes())
  
  p21 <- Reduce( `+`, p21 )+patchwork::plot_layout( ncol = 1 )
  
  p22 <- Reduce( `+`, p22 )+patchwork::plot_layout( ncol = 4 )
  
  p23 <- list(p21, p22)
  
  
  design <- c(patchwork::area(1, 1, 1, 1), patchwork::area(1, 2, 1, 5))
  
  p24 <- Reduce( `+`,  p23)+patchwork::plot_layout(design = design)
  
  
  myplots[[i]] <- local({
    i <- i
    p24
  })
  
  
}

p20 <- patchwork::wrap_plots(myplots, ncol = 1)

p20

#myplots[1:nrow(genes)] + patchwork::plot_layout(byrow = T, widths = 25, heights = 5)
cairo_pdf("Plots/Figure_4/Suplementary/Hypoxia_FeaturePlots.pdf",  width = 25, height = 90, family = "Arial")
p20
dev.off()


###############################################################

#  GSEA Heatmaps. HallMark Analysis

Liver <- readRDS("S:/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Liver_MacarenaFC/Figures_Paper/rds/Liver_20Jan21.Figures.All.Conditions.ppt_Rui.rds")

table(Liver@meta.data$Condition)

########################################################

# Dll4/aVEGF vs Dll4/Myc


Dll4aVEGFvDll4Myc <- subset(Liver, subset = Condition == "D4KO_aVEGF" | Condition == "Dll4/MycKO")

cond.genes <- wilcoxauc(Dll4aVEGFvDll4Myc, 'Condition')

# Dll4 vs Control

Dll4vCtl <- subset(Liver, subset = Condition == "Dll4KO" | Condition == "Control")

cond.genes <- wilcoxauc(Dll4vCtl, 'Condition')

# Notch1vs Control 

Notch1vCtl <- subset(Liver, subset = Condition == "NOTCH1KO" | Condition == "Control")

cond.genes <- wilcoxauc(Notch1vCtl, 'Condition')

# RbpjvsControl

RbpjvCtl <- subset(Liver, subset = Condition == "RBPJKO" | Condition == "Control")

cond.genes <- wilcoxauc(RbpjvCtl, 'Condition')

# we have all the genes for each cluster
dplyr::count(cond.genes, group)

msigdbr_species()

m_df<- msigdbr(species = "Mus musculus", category = 13)

h_df<- msigdbr(species = "Homo sapiens", category = 7)


table(m_df$gs_subcat)

table(h_gene_sets$gs_cat)

# Using Hallmark genes only

h_gene_sets = msigdbr(species = "mouse", category = "H")

head(m_df)

fgsea_sets<- h_gene_sets %>% split(x = .$gene_symbol, f = .$gs_name)

fgsea_sets$HALLMARK_MYC_TARGETS_V1


# Myc Targets


cond.genes %>%
  dplyr::filter(group == "NOTCH1KO") %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 20)

# Good to see ODC1, KCNE3, PEG10, TPI1, MIF

# select only the feature and auc columns for fgsea, which statistics to use is an open question

cond.dll4.genes<- cond.genes %>%
  dplyr::filter(group == "NOTCH1KO") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)


cond.dll4.genes$feature <- tolower(cond.dll4.genes$feature) %>% str_to_title()

genesTables <- cond.dll4.genes %>%
  +   mutate(rank = rank(cond.dll4.genes$padj,  ties.method = "random")) 

ranks<- deframe(cond.dll4.genes)

head(ranks)

ranks

fgseaRes<- fgsea(fgsea_sets, stats = ranks)

# Tidy Data with auc

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))


fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  head()

write.xlsx(fgseaResTidy, "Tables/Hypoxia/fgseaResTidy.NOTCH1KOvCtl.auc.xlsx")

capture.output(summary(fgsea_sets), file = "My New File.txt")

# BarPlot

# only plot the top 20 pathways

ggplot(fgseaResTidy %>% filter(padj < 0.05) %>% head(n= 20), aes(reorder(pathway, NES), NES))

ggplot(fgseaResTidy %>% filter(padj < 0.05) %>% head(n= 20), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES < 7.5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()

# GSEA Style Plot

plotEnrichment(fgsea_sets[["HALLMARK_COMPLEMENT"]],
               ranks, gseaParam = 1) + labs(title="HALLMARK_COMPLEMENT")

fgsea_sets[["RBM34_TARGET_GENES"]]
fgsea_sets[["MYC_UP.V1_UP"]]


p1 <- plotEnrichment(fgsea_sets[["HALLMARK_HYPOXIA"]],
                     ranks, ticksSize = 1) + labs(title="HALLMARK_HYPOXIA")

cairo_pdf("Plots/Figure_4/Suplementary/GSEA_Hypoxia_RbpjvsCtl.pdf",  width = 7, height = 6, family = "Arial")
p1
dev.off()


#  Heatmap

cond.fgsea <- read.xlsx("./Tables/Hypoxia/fgseaResTidy.Dll4KO_Notchq_RbpjvCtl.auc.xlsx", sheet = 4, sep.names = " ", rowNames = T)

liver.gs <- subset(Liver, subset = Condition == "Dll4KO" | Condition == "NOTCH1KO" | Condition == "RBPJKO")

Idents(liver.gs) <- "Condition"

avgexp <- AverageExpression(liver.gs, return.seurat = T)

avgexp@meta.data$Condition <- avgexp@active.ident

avgexp <- Seurat::AddMetaData(avgexp, cond.fgsea)

avgexp@meta.data$Condition <- avgexp@active.ident

levels(avgexp) <- c("Dll4KO", "NOTCH1KO", "RBPJKO")

avgexp@meta.data$Condition <- avgexp@active.ident

GS.hallmark.HS <- getGeneSets(library = "H")

ES.seurat <- enrichIt(obj = avgexp, gene.sets = GS.hallmark.HS, groups = 3, cores = 2)

names(ES.seurat)

cairo_pdf("Plots/Hypoxia/Figure_GSEA_Heatmap_Dll4vsControl_Notch1vsControl_RbpjvControl_unscaled.pdf", width = 5, height = 12, family = "Arial")
dittoHeatmap(avgexp, genes = NULL, metas = names(cond.fgsea), 
             annot.by = "Condition",
             # order.by = c(2,1, 3),
             annot.colors = palette_Maca_GSEA,
             fontsize = 7, 
             cluster_cols = F,
             cluster_rows = F,
             heatmap.colors = BuRd,
             # heatmap.colors.max.scaled = BuRd,
             scaled.to.max = F,
             scale = "none",
             complex = T, data.out = F)
# legend_breaks = c(1, 0.8, 0.6, 0.4, 0.2, 0),
# main = expression(Control~italic(Dll4)^iDEC~italic(Notch1)^iDEC~italic(Rbpj)^iDEC))
dev.off()



