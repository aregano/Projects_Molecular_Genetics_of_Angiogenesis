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

#############################################################

colors <- colorRampPalette(c("#0348A6", "#7AC5FF", "#C6FDEC", "#FFB433", "#FF4B20"))
colors <- c("#0348A6", "#7AC5FF", "#C6FDEC", "#FFB433", "#FF4B20")
my_palette_Rui_colors <- c("#E95A74", "#50B6EF", "#45FF8E", "#F4A753", "#A80519", "#880088", "#E28CF4", "#C1B80C","#FC0808",  "#0E47D8")
my_palette_Rui_colors_B <- c("#F4A753", "#E95A74", "#50B6EF", "#FC0808",  "#0E47D8", "#45FF8E", "#A80519", "#E28CF4", "#880088", "#C1B80C")
BuRd <- c("#052F60", "#2066AC", "#4392C3", "#92C5DE", "#D1E4EF", "#F7F7F7", "#FCDBC6", "#F4A482", "#D55F4D", "#B1172B", "#67001F")

#############################################################

# Trying GSEA with a for loop

Liver <- readRDS("S:/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Liver_MacarenaFC/Figures_Paper/rds/Liver_20Jan21.Figures.rds")

table(Liver@meta.data$Condition)

m_df<- msigdbr(species = "Mus musculus", category = 13)

h_df<- msigdbr(species = "Homo sapiens", category = 7)

h_gene_sets = msigdbr(species = "mouse", category = "H")

head(m_df)

fgsea_sets<- h_gene_sets %>% split(x = .$gene_symbol, f = .$gs_name)

fgsea_sets$HALLMARK_MYC_TARGETS_V1


liver.genes <- wilcoxauc(Liver, 'FigClustering')

dplyr::count(cluster.genes, group)

figclustering <- c("C0 - Unspecified quiescent capillaries", "C1a - Arterial capillaries", "C1v - Venous capillaries", "C2a - Large arteries", "C2v - Large veins",  "C3 - Activated capillaries", 
                    "C4 - Endothelial tip cells", "C5ip - Proliferating S-phase", "C5p - G2-M", "C6 - X")


liver.genes %>%
  dplyr::filter(group == figclustering[5]) %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 20)


myfgsea <- vector('list', 10)


for (i in 1:10) {
  
  cluster.genes<- liver.genes %>%
    dplyr::filter(group == figclustering[i]) %>%
    arrange(desc(auc)) %>% 
    dplyr::select(feature, auc)
  
  cluster.genes$feature <- tolower(cluster.genes$feature) %>% str_to_title()
  
  ranks<- deframe(cluster.genes)
  
 myfgsea[[i]] <- local({
    i <- i
    fgseaRes<- fgsea(fgsea_sets, stats = ranks)
  })
  
}

myfgsea[1]

write.xlsx(myfgsea, file = "./Tables/GSEA_Fig4/fgsea_Dll4KO_Clusters.xlsx")

write.xlsx(myfgsea, file = "./Tables/fgsea_Clusters.xlsx")

# select only the feature and auc columns for fgsea, which statistics to use is an open question

cluster.fgsea <- read.xlsx("./Tables/fgsea_Clusters.xlsx ", sheet = 12, sep.names = " ", rowNames = T)

selected_HallmarkGS <- read.csv2("./Tables/Hallmark_GeneSets.csv", header = F)

selected_HallmarkGS <- as.list(selected_HallmarkGS)

View(cluster.fgsea)

GS.hallmark.HS <- getGeneSets(library = "H")

# GS.hallmark.HS.selected <- GS.hallmark.HS[(colnames(GS.hallmark.HS) %in% selected.HallmarkGS), ]

avgexp <- AverageExpression(Liver, return.seurat = T)

ES.seurat <- enrichIt(obj = avgexp, gene.sets = GS.hallmark.HS, groups = 10, cores = 2)

ES.seurat <- ES.seurat[, (colnames(ES.seurat) %in% selected.HallmarkGS)]

selected.HallmarkGS <- c("HALLMARK_E2F_TARGETS",
                         "HALLMARK_G2M_CHECKPOINT",
                         "HALLMARK_MITOTIC_SPINDLE",
                         "HALLMARK_MYC_TARGETS_V1",
                         "HALLMARK_MYC_TARGETS_V2",
                         "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
                         "HALLMARK_MTORC1_SIGNALING",
                         "HALLMARK_GLYCOLYSIS",
                         "HALLMARK_FATTY_ACID_METABOLISM",
                         "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
                         "HALLMARK_PROTEIN_SECRETION",
                         "HALLMARK_HYPOXIA",
                         "HALLMARK_P53_PATHWAY",
                         "HALLMARK_DNA_REPAIR",
                         "HALLMARK_APOPTOSIS",
                         "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
                         "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
                         "HALLMARK_INFLAMMATORY_RESPONSE"
)


ES.seurat.selected <- ES.seurat[, (colnames(ES.seurat) %in% selected.HallmarkGS)]

ES.seurat.selected <- ES.seurat.selected[,selected.HallmarkGS]



cluster.fgsea.selected <- cluster.fgsea[, (colnames(cluster.fgsea) %in% selected.HallmarkGS)]


gene.set <- sample(x = rownames(x = object@data), size = 100, replace = FALSE)


avgexp <- Seurat::AddMetaData(avgexp, cluster.fgsea)

levels(avgexp) <- c("C0 - Unspecified quiescent capillaries", "C1a - Arterial capillaries", "C1v - Venous capillaries", "C2a - Large arteries", "C2v - Large veins",  "C3 - Activated capillaries", 
                    "C4 - Endothelial tip cells", "C5ip - Proliferating S-phase", "C5p - G2-M", "C6 - X","C1")

avgexp@meta.data$FigClustering <- avgexp@active.ident

metas.fgsea <- getMetas(avgexp, names.only = T)

metas.fgsea <- metas.fgsea[10:27]

selected.HallmarkGS <- as.vector(selected.HallmarkGS)

cairo_pdf("Plots/GSEA/Heatmaps/GSEA_Heatmap_Liver_Ctrl_Dll4_Notch1_Rbpj_unscaled.pdf", width = 10, height = 12, family = "Arial")
dittoHeatmap(avgexp, genes = NULL, metas = names(ES.seurat), 
             annot.by = "FigClustering",
             # order.by = c(),
             annot.colors = my_palette_Rui_colors_B,
             fontsize = 7, 
             cluster_cols = F,
             cluster_rows = F,
             heatmap.colors = BuRd,
             # heatmap.colors.max.scaled = BuRd,
             scaled.to.max = F,
             scale = "none",
             complex = F, data.out = F,
             # legend_breaks = c(1, 0.8, 0.6, 0.4, 0.2, 0),
             main = expression(Control~italic(Dll4)^iDEC~italic(Notch1)^iDEC~italic(Rbpj)^iDEC))
dev.off()


cairo_pdf("Plots/GSEA/Heatmaps/GSEA_Heatmap_selected_Liver_Ctrl_Dll4_Notch1_Rbpj_unscaled.pdf", width = 10, height = 6, family = "Arial")
dittoHeatmap(avgexp, genes = NULL, metas = names(ES.seurat.selected), 
             annot.by = "FigClustering",
             annot.colors = my_palette_Rui_colors_B,
             fontsize = 7, 
             cluster_cols = F,
             cluster_rows = F,
             heatmap.colors = BuRd,
             # heatmap.colors.max.scaled = BuRd,
             scaled.to.max = F,
             scale = "none",
             complex = F, data.out = F,
             # legend_breaks = c(1, 0.8, 0.6, 0.4, 0.2, 0),
             main = expression(Control~italic(Dll4)^iDEC~italic(Notch1)^iDEC~italic(Rbpj)^iDEC))
dev.off()

cairo_pdf("Plots/GSEA/Heatmaps/GSEA_Heatmap_Liver_Ctrl_Dll4_Notch1_Rbpj_scaled.pdf", width = 10, height = 12, family = "Arial")
dittoHeatmap(avgexp, genes = NULL, metas = names(ES.seurat), 
             annot.by = "FigClustering",
             annot.colors = my_palette_Rui_colors_B,
             fontsize = 7, 
             cluster_cols = F,
             cluster_rows = F,
             heatmap.colors = BuRd,
             # heatmap.colors.max.scaled = BuRd,
             scaled.to.max = F,
             scale = "row",
             complex = F, data.out = F,
             # legend_breaks = c(1, 0.8, 0.6, 0.4, 0.2, 0),
             main = expression(Control~italic(Dll4)^iDEC~italic(Notch1)^iDEC~italic(Rbpj)^iDEC))
dev.off()

cairo_pdf("Plots/GSEA/Heatmaps/GSEA_Heatmap_selected_Liver_Ctrl_Dll4_Notch1_Rbpj_scaled.pdf", width = 10, height = 6, family = "Arial")
dittoHeatmap(avgexp, genes = NULL, metas = names(ES.seurat.selected), 
             annot.by = "FigClustering",
             annot.colors = my_palette_Rui_colors_B,
             fontsize = 7, 
             cluster_cols = F,
             cluster_rows = F,
             heatmap.colors = BuRd,
             # heatmap.colors.max.scaled = BuRd,
             scaled.to.max = F,
             scale = "row",
             complex = F, data.out = F,
             # legend_breaks = c(1, 0.8, 0.6, 0.4, 0.2, 0),
             main = expression(Control~italic(Dll4)^iDEC~italic(Notch1)^iDEC~italic(Rbpj)^iDEC))
dev.off()

data_heatmap <- dittoHeatmap(avgexp, genes = NULL, metas = names(ES.seurat), 
                             annot.by = "FigClustering",
                             annot.colors = my_palette_Rui_colors_B,
                             fontsize = 7, 
                             cluster_cols = F,
                             cluster_rows = F,
                             heatmap.colors = BuRd,
                             complex = T, data.out = T)

list_mat_heatmap <- data_heatmap[["mat"]]

list_mat_heatmap

list_mat_heatmap.m = melt(list_mat_heatmap)
head(list_mat_heatmap.m)

list_mat_heatmap.m$value <- round(list_mat_heatmap.m$value, 2)

head(list_mat_heatmap.m)


cairo_pdf("Plots/GSEA/Heatmaps/GSEA_Heatmap_Liver_Ctrl_Dll4_Notch1_Rbpj_scores.pdf", width = 14, height = 12, family = "Arial")
ggplot(data = data.frame(list_mat_heatmap.m), aes(x = Var2, y = fct_rev(as_factor(Var1)), fill= value)) + geom_tile(color = "white")+ geom_text(aes(label = value), size = 3)+
  theme_minimal() + theme(text = element_text(size = 10), axis.title = element_blank(), axis.text.x.bottom = element_text(angle = 90, hjust = 1))+scale_fill_gradientn(colours = BuRd)
dev.off()


dittoHeatmap(avgexp, genes = NULL, metas = names(), 
             annot.by = "FigClustering",
             annot.colors = my_palette_Rui_colors_B,
             fontsize = 7, 
             cluster_cols = F,
             heatmap.colors = BuRd,
             complex = T, data.out = F)


jpeg("Plots/GSEA/GSEA_Heatmap_Liver_Ctrl_Dll4_Notch1_Rbpj.jpeg", width = 14, height = 12, units = 'in', res = 800)
dittoHeatmap(avgexp, genes = NULL, metas = names(ES.seurat), 
             annot.by = "FigClustering",
             annot.colors = my_palette_Rui_colors_B,
             fontsize = 7, 
             cluster_cols = F,
             heatmap.colors = BuRd)
dev.off()


