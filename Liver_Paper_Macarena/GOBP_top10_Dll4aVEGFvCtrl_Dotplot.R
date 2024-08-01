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

#  Comparing ranks Hypoxia genes

Liver <- readRDS("//Tierra/SC/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Liver_MacarenaFC/Figures_Paper/rds/Liver_20Jan21.Figures.All.Conditions.OldDll4Myc.rds")

table(Liver@meta.data$Condition)

Liver <- subset(Liver, subset = Condition == "Control" | Condition == "Dll4KO" | Condition == "Dll4/MycKO" | Condition == "D4KO_aVEGF")

Idents(Liver) <- "Condition"

levels(Liver) <- c("Control", "Dll4KO", "Dll4/MycKO", "D4KO_aVEGF")

Liver@active.ident -> Liver@meta.data$Condition

table(Liver@meta.data$Condition)

ranks <- read.csv2("Tables/GSEA_Fig7/ranks.Dll4aVEGFvsCtrl.BP.csv", sep = ",")

ranksdll4 <- read.csv2("Tables/GSEA_Fig7/ranks.Dll4vsCtrl.BP.csv", sep = ",")

Bestholtz_palette <- c("#DEDAD6", "#FEE392", "#FEC44E", "#FE9929", "#ED6F11", "#CC4C17", "#993411", "#65260C")

# Load fgsea dataset

msigdbr_species()

m_df<- msigdbr(species = "Mus musculus", category = 13)

h_df<- msigdbr(species = "Homo sapiens", category = 7)


table(m_df$gs_subcat)

table(h_gene_sets$gs_cat)

# Using Hallmark genes only

h_gene_sets = msigdbr(species = "mouse", subcategory = "GO:BP")

fgsea_sets<- h_gene_sets %>% split(x = .$gene_symbol, f = .$gs_name) %>% list.filter("Type" > 50)

fgsea_sets$GOBP_AGING
head(m_df)

Cellular_Macromolecule_BS_Process_genes <- fgsea_sets$GOBP_CELLULAR_MACROMOLECULE_BIOSYNTHETIC_PROCESS

RNA_Processing_genes <- fgsea_sets$GOBP_RNA_PROCESSING

Cell_Cycle_genes <- fgsea_sets$GOBP_CELL_CYCLE

Cellular_Macromolecule_BS_Process_genes <- as.data.frame(Cellular_Macromolecule_BS_Process_genes)

RNA_Processing_genes <- as.data.frame(RNA_Processing_genes)

Cell_Cycle_genes <- as.data.frame(Cell_Cycle_genes)

library(tidyverse)

# Getting the ranking order for the 3 GOBP

Cellular_Macromolecule_BS_Process_genes <- t(Cellular_Macromolecule_BS_Process_genes)

Cellular_Macromolecule_BS_Process_genes <- as.vector(Cellular_Macromolecule_BS_Process_genes)

Cellular_Macromolecule_rank <- ranks[ ranks$X %in% Cellular_Macromolecule_BS_Process_genes, ]

RNA_Processing_genes <- t(RNA_Processing_genes)

RNA_Processing_genes <- as.vector(RNA_Processing_genes)

RNA_Processing_rank <- ranks[ ranks$X %in% RNA_Processing_genes, ]

Cell_Cycle_genes <- t(Cell_Cycle_genes)

Cell_Cycle_genes <- as.vector(Cell_Cycle_genes)

Cell_Cycle_rank <- ranks[ ranks$X %in% Cell_Cycle_genes, ]

Cell_Cycle_rank <- ranks[ ranksdll4$X %in% Cell_Cycle_genes, ]

# Top 10 Hypoxia genes

Idents(Liver) <- "Condition"

levels(Liver) <- c("D4KO_aVEGF", "Dll4/MycKO", "Dll4KO", "Control")

Top10_Cellular_Macromolecule <- Cellular_Macromolecule_rank$X[c(1:10)] %>% toupper()

Top10_Cellular_Macromolecule_H <- Cellular_Macromolecule_rank$X[c(1:10)]

Top10_RNA_Processing <-RNA_Processing_rank$X[c(1:10)] %>% toupper()

Top10_RNA_Processing_H <- RNA_Processing_rank$X[c(1:10)]

Top10_Cell_Cycle <-head(Cell_Cycle_rank$X, 10) %>% toupper()

Top10_Cell_Cycle_H <- head(Cell_Cycle_rank$X, 10)

class(Top10_Cell_Cycle)

Top10_GOBPs <- c(Top10_Cellular_Macromolecule, Top10_RNA_Processing, Top10_Cell_Cycle) %>% unique()

Top10_GOBPs_H <- c(Top10_Cellular_Macromolecule_H, Top10_RNA_Processing_H, Top10_Cell_Cycle_H) %>% unique()


cairo_pdf("Plots/GOBP/Figure_DotPlot_Maca_GOBP_genes_Dll4aVEGFvCtrl.Top10_raw.pdf", width = 9, height = 2, family = "Arial")
DotPlot(Liver, features = Top10_GOBPs, col.min = 0, dot.scale = 5, scale = T)+
  theme(strip.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"),
        legend.title = element_text(size = 10), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"), 
        plot.margin = unit(c(0, 0, 0, 0), "cm")
  ) +
  scale_y_discrete(labels=c(bquote(italic(Dll4)^iDEC+aVEGF), bquote(italic(Dll4/Myc)^iDEC), bquote(italic(Dll4)^iDEC), "Control"))+
  scale_x_discrete(label = Top10_GOBPs_H)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  geom_vline(xintercept = c(10.5, 17.5),linetype = 2 )+
  theme(legend.box = "horizontal")+
  theme(legend.box = "horizontal", title = element_text(size = 10))+
  ggtitle("GOBP: Cellular Macromolecule Biosynthetic Process; RNA Processing; Cell Cycle")
dev.off()
#

# Comparing repeated entries

Top10_GOBPs_check <- c(Top10_Cellular_Macromolecule, Top10_RNA_Processing, Top10_Cell_Cycle)

max_length <- max(c(length(Top10_GOBPs_check), length(Top10_GOBPs)))    # Find out maximum length
max_length 

Top10_GOBPs_check_1 <- data.frame(col1 = c(Top10_GOBPs_check),                 # Create data frame with unequal vectors
                   col2 = c(Top10_GOBPs,
                            rep(NA, max_length - length(Top10_GOBPs))))
Top10_GOBPs_check_1                                            # Print final da

# Repeated values are: 
# GOBP_RNA_Processing -> RPS8, RPL35A, RPS7
#  I will add to get those 10 per category

Top10_Cellular_Macromolecule <- Cellular_Macromolecule_rank$X[c(1:10)] %>% toupper()

Top10_Cellular_Macromolecule_H <- Cellular_Macromolecule_rank$X[c(1:10)]

Top10_RNA_Processing <-RNA_Processing_rank$X[c(1:13)] %>% toupper()

Top10_RNA_Processing_H <- RNA_Processing_rank$X[c(1:13)]

Top10_Cell_Cycle <-head(Cell_Cycle_rank$X, 15) %>% toupper()

Top10_Cell_Cycle_H <- head(Cell_Cycle_rank$X, 15)

Top10_GOBPs <- c(Top10_Cellular_Macromolecule, Top10_RNA_Processing, Top10_Cell_Cycle) %>% unique()

Top10_GOBPs_H <- c(Top10_Cellular_Macromolecule_H, Top10_RNA_Processing_H, Top10_Cell_Cycle_H) %>% unique()

# Do the check

Top10_GOBPs_check <- c(Top10_Cellular_Macromolecule, Top10_RNA_Processing, Top10_Cell_Cycle)

max_length <- max(c(length(Top10_GOBPs_check), length(Top10_GOBPs)))    # Find out maximum length
max_length 

Top10_GOBPs_check_1 <- data.frame(col1 = c(Top10_GOBPs_check),                 # Create data frame with unequal vectors
                                  col2 = c(Top10_GOBPs,
                                           rep(NA, max_length - length(Top10_GOBPs))))
Top10_GOBPs_check_1                                            # Print final da

cairo_pdf("Plots/GOBP/Figure_DotPlot_Maca_GOBP_genes_Dll4aVEGFvCtrl.Top10_each_category.pdf", width = 10, height = 2, family = "Arial")
DotPlot(Liver, features = Top10_GOBPs, col.min = 0, dot.scale = 5, scale = T)+
  theme(strip.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"),
        legend.title = element_text(size = 10), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"), 
        plot.margin = unit(c(0, 0, 0, 0), "cm")
  ) +
  scale_y_discrete(labels=c(bquote(italic(Dll4)^iDEC+aVEGF), bquote(italic(Dll4/Myc)^iDEC), bquote(italic(Dll4)^iDEC), "Control"))+
  scale_x_discrete(label = Top10_GOBPs_H)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  geom_vline(xintercept = c(10.5, 20.5),linetype = 2 )+
  theme(legend.box = "horizontal", title = element_text(size = 10))+
  ggtitle("GOBP: Cellular Macromolecule Biosynthetic Process; RNA Processing; Cell Cycle")
dev.off()
#

###############################################################################################################################

# Top10 Cellular Macromolecule Biosynthetic Process + RNA Processing


Top10_GOBPs <- c(Top10_Cellular_Macromolecule, Top10_RNA_Processing) %>% unique()

Top10_GOBPs_H <- c(Top10_Cellular_Macromolecule_H, Top10_RNA_Processing_H) %>% unique()

# Do the check

Top10_GOBPs_check <- c(Top10_Cellular_Macromolecule, Top10_RNA_Processing)

max_length <- max(c(length(Top10_GOBPs_check), length(Top10_GOBPs)))    # Find out maximum length
max_length 

Top10_GOBPs_check_1 <- data.frame(col1 = c(Top10_GOBPs_check),                 # Create data frame with unequal vectors
                                  col2 = c(Top10_GOBPs,
                                           rep(NA, max_length - length(Top10_GOBPs))))
Top10_GOBPs_check_1                                            # Print final da

cairo_pdf("Plots/GOBP/Figure_DotPlot_Maca_GOBP_genes_Dll4aVEGFvCtrl.Top10_Biosynthsis_&_RNA_Processing.pdf", width = 8, height = 2, family = "Arial")
DotPlot(Liver, features = Top10_GOBPs, col.min = 0, dot.scale = 5, scale = T)+
  theme(strip.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"),
        legend.title = element_text(size = 10), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"), 
        plot.margin = unit(c(0, 0, 0, 0), "cm")
  ) +
  scale_y_discrete(labels=c(bquote(italic(Dll4)^iDEC+aVEGF), bquote(italic(Dll4/Myc)^iDEC), bquote(italic(Dll4)^iDEC), "Control"))+
  scale_x_discrete(label = Top10_GOBPs_H)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  geom_vline(xintercept = c(10.5),linetype = 2 )+
  theme(legend.box = "horizontal", title = element_text(size = 10))+
  ggtitle("GOBP: Cellular Macromolecule Biosynthetic Process; RNA Processing")
dev.off()


###############################################################################################################################

# Standalone Cell cycle plot

Top10_Cell_Cycle_H <- c("Mcm6", "Mcm5", "Lig1", "Mcm2", "Mcm3", "Top2a", "Stmn1", "Mki67", "Hmgb2", "Cenpf")

Top10_Cell_Cycle <- Top10_Cell_Cycle_H %>% toupper()

cairo_pdf("Plots/GOBP/Figure_DotPlot_Maca_GOBP_genes_Dll4aVEGFvCtrl.Top10_Cell_Cycle.pdf", width = 7, height = 2, family = "Arial")
DotPlot(Liver, features = Top10_Cell_Cycle, col.min = 0, dot.scale = 10, scale = T)+
  theme(strip.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"),
        legend.title = element_text(size = 10), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"), 
        plot.margin = unit(c(0, 0, 0, 0), "cm")
  ) +
  scale_y_discrete(labels=c(bquote(italic(Dll4)^iDEC+aVEGF), bquote(italic(Dll4/Myc)^iDEC), bquote(italic(Dll4)^iDEC), "Control"))+
  scale_x_discrete(label = Top10_Cell_Cycle_H)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  theme(legend.box = "horizontal", title = element_text(size = 10))+
  ggtitle("GOBP: Cell Cycle")
dev.off()
