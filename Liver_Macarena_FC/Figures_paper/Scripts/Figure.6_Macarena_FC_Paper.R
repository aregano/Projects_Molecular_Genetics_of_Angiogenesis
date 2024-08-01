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
library(fgsea)
library(SingleCellExperiment)
library(escape)
library(reshape2)
library(clusterProfiler)
library(pathview)
library(enrichplot)
library(presto)
library(msigdbr)
library("tidyselect")
library(magrittr)

# rds

Liver <- readRDS("//Tierra/SC/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Liver_MacarenaFC/Figures_Paper/rds/Liver_20Jan21.Figures.All.Conditions.OldDll4Myc.rds")
table(Liver@meta.data$Condition)

Liver <- subset(Liver, subset = Condition == "Control" | Condition == "Dll4KO" | Condition == "Dll4/MycKO")

LiverDll4 <- subset(Liver, Condition == "Dll4KO")

Idents(Liver) <- Liver@meta.data$Condition

levels(Liver) <- c("Control", "Dll4KO", "Dll4/MycLOF")

Liver@active.ident -> Liver@meta.data$Condition

# Color Palettes

my_palette_Rui_colors <- c("#E95A74", "#50B6EF", "#45FF8E", "#F4A753", "#A80519", "#880088", "#E28CF4", "#C1B80C","#FC0808",  "#0E47D8")

my_palette_Rui_colors_B <- c("#F4A753", "#E95A74", "#50B6EF", "#FC0808",  "#0E47D8", "#45FF8E", "#A80519", "#E28CF4", "#880088", "#C1B80C")

custom_rdylblu_palette <- c("#005AD9", "#0070FF", "#FFE135", "#FFDA13", "#E60026", "#C60012")

macarena_palette <- c("#FDDC85", "#6E2914")

Bestholtz_palette <- c("#DEDAD6", "#FEE392", "#FEC44E", "#FE9929", "#ED6F11", "#CC4C17", "#993411", "#65260C")

palette_Maca_Vln <- c("#606060", "#F94040", "#05BE78", "#FFA040")

palette_Maca_GSEA <- c("#F94040", "mediumpurple1")

BuRd <- c("#052F60", "#2066AC", "#4392C3", "#92C5DE", "#D1E4EF", "#F7F7F7", "#FCDBC6", "#F4A482", "#D55F4D", "#B1172B", "#67001F")

###############################################################################################################################

## add the Arial font
font_add("Arial", regular = "arial.ttf",
         bold = "arialbd.ttf", italic = "ariali.ttf", bolditalic = "arialbi.ttf")

###############################################################################################################################

# Figure 6C. Pie Chart

# PieChart

library(ggplot2)
library(ggrepel)
library(openxlsx)
library(lessR)
library(Seurat)
library(reshape2)



my_palette_Rui_colors_B <- c("#F4A753", "#E95A74", "#50B6EF", "#FC0808",  "#0E47D8", "#45FF8E", "#A80519", "#E28CF4", "#880088", "#C1B80C")


cluster.fgsea <- read.xlsx("./Tables/fgsea_Clusters_Ctrl_Dll4KO_RbpjKO_Notch1KO.xlsx ", sheet = 13, sep.names = " ", rowNames = T)

cluster.fgsea <- round(cluster.fgsea, digits = 2)

gene.sets <- c(Abs_HALLMARK_E2F_TARGETS, Abs_HALLMARK_MYC_TARGETS_V1, Abs_HALLMARK_MYC_TARGETS_V2, Abs_HALLMARK_G2M_CHECKPOINT, Abs_HALLMARK_OXIDATIVE_PHOSPHORYLATION)

gene.sets.names <- c("E2F TARGETS", "MYC TARGETS V1", "MYC TARGETS V2", "G2M CHECKPOINT","OXIDATIVE PHOSPHORYLATION")

# YOu need to change the y label to the Abs value of the NES

class(cluster.fgsea)

bp<- ggplot(cluster.fgsea, aes(x="", y=Abs_HALLMARK_E2F_TARGETS, fill=lbls))+
  geom_bar(width = 1, stat = "identity")

bp

pie <- bp + coord_polar("y", start=0)
pie

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )

# Need to change the label in aes to the real NES score and set up the ggtitle

p1 <- pie + scale_fill_manual(values=my_palette_Rui_colors_B)+ 
  blank_theme +
  geom_text(aes(label = HALLMARK_E2F_TARGETS),
            position = position_stack(vjust = 0.5), show.legend = F)+
  theme_void()+NoLegend()+ggtitle(gene.sets.names[1])


cairo_pdf("Plots/Figure_6/Pie_Charts/Figure_6C_Pie_Chart_E2F_TARGETS.pdf", width = 5, height = 5, family = "Arial")

p1

dev.off()  

###############################################################################################################################

# Figure 6H.1 UMAP

DimPlot(Liver, label = T)

Idents(Liver) <- Liver@meta.data$FigClustering

p1 <- DimPlot(Liver, cols = my_palette_Rui_colors_B, pt.size = 1.2, group.by = "FigClustering", split.by = "Condition", combine = T)

p1

new_labels <- c("Control" = "Control", "Dll4KO" = "italic(Dll4)^iDEC", "Dll4/MycLOF" = "italic(Dll4/Myc)^iDEC")

p21 <- p1+facet_grid(~Condition, labeller = as_labeller(new_labels, default = label_parsed))+
  theme(strip.text.x = element_text(size = 24), plot.title = element_blank())+NoLegend()+xlab("UMAP_1")+ylab("UMAP_2")

p22 <- DimPlot(Liver, group.by = "FigClustering", cols = my_palette_Rui_colors_B, pt.size = 1.2)+theme(legend.text = element_text(size = 14), axis.title.x = element_blank(), axis.line.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks = element_blank(), plot.title = element_blank())

p21
p22

p23 <- list(p21, p22)

design <- c(patchwork::area(1, 1, 1, 3), patchwork::area(1, 4, 1, 4.5))

p24 <- Reduce( `+`,  p23)+patchwork::plot_layout(design = design)


cairo_pdf("Plots/Figure_6/Figure_6H.1_UMAP.pdf",  width = 29, height = 6, family = "Arial")
p24
dev.off()


jpeg("Plots/Figure_6/Figure_6H.1_UMAP.jpeg", width = 29, height = 6, units = 'in', res = 800)
p24
dev.off()



###############################################################################################################################


# Figure 6H.2 Bar Plot

table(Liver@meta.data$Condition)

p1 <- dittoBarPlot(Liver, var = "FigClustering", group.by = "Condition", x.reorder = c(1,3,2), color.panel = my_palette_Rui_colors_B, xlab = NULL, main = NULL, y.breaks = c(0,0.25,0.5, 0.75, 1), min = 0, max = 1)+
  theme(axis.text=element_text(size=15), axis.title = element_text(size = 18), legend.text = element_text(size = 16), axis.text.x = element_text(angle = 25))+
  scale_x_discrete(labels=c("Control", bquote(italic(Dll4)^iDEC), bquote(italic(Dll4/Myc)^iDEC)))+
  geom_col(width = 0.1)+
  # scale_y_discrete(labels=c("0", "25", "50", "75", "100"))
  scale_y_continuous(labels = c("0", "25", "50", "75", "100"))
p1

cairo_pdf("Plots/Figure_6/Figure_6H.2_BarPlot.pdf",  width = 6.5, height = 6, family = "Arial")
p1
dev.off()


jpeg("Plots/Figure_6/Figure_6H.2_BarPlot.jpeg", width = 6.5, height = 6, units = 'in', res = 800)
p1
dev.off()

###############################################################################################################################

# Figure 6I GSEA Heatmap

LiverD4KO <- subset(Liver, subset = Condition == "Control" | Condition == "Dll4KO")
LiverD4MycKO <- subset(Liver, subset = Condition == "Control" | Condition == "Dll4/MycKO")

LiverD4KO@meta.data$Condition <- LiverD4KO@active.ident

Idents(LiverD4MycKO) <- "Condition"
LiverD4MycKO@meta.data$Condition <- LiverD4MycKO@active.ident


GS.hallmark.HS <- getGeneSets(library = "H")

liver.genes <- wilcoxauc(Liver, 'Condition')
liverd4.genes <- wilcoxauc(LiverD4KO, 'Condition')
liverd4myc.genes <- wilcoxauc(LiverD4MycKO, 'Condition')

# we have all the genes for each cluster
dplyr::count(liverd4.genes, group)
dplyr::count(liverd4myc.genes, group)

msigdbr_show_species()

m_df<- msigdbr(species = "Mus musculus", category = 13)

h_df<- msigdbr(species = "Homo sapiens", category = 7)

# Using Hallmark genes only

h_gene_sets = msigdbr(species = "mouse", category = "H")

head(m_df)

fgsea_sets<- h_gene_sets %>% split(x = .$gene_symbol, f = .$gs_name)

fgsea_sets$HALLMARK_MYC_TARGETS_V1

table(Liver@meta.data$Condition)

# Myc Targets


liverd4.genes %>% 
  dplyr::filter(group == "Dll4KO") %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 20)

liverd4myc.genes %>%
  dplyr::filter(group == "Dll4/MycKO") %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 20)

liverd4.sel.genes<- liverd4.genes %>%
  dplyr::filter(group == "Dll4KO") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)

liverd4myc.sel.genes<- liverd4myc.genes %>%
  dplyr::filter(group == "Dll4/MycKO") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)

# Good to see ODC1, KCNE3, PEG10, TPI1, MIF

# select only the feature and auc columns for fgsea, which statistics to use is an open question

liverd4.genes$feature <- tolower(liverd4.genes$feature) %>% str_to_title()

liverd4myc.genes$feature <- tolower(liverd4myc.genes$feature) %>% str_to_title()


cond.d4.genes <- liverd4.genes %>%
  dplyr::filter(group == "Dll4KO") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)

cond.d4myc.genes <- liverd4myc.genes %>%
  dplyr::filter(group == "Dll4/MycKO") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc) 

genesTablesd4 <- cond.d4.genes %>% mutate(rank = rank(cond.d4.genes$auc,  ties.method = "random")) 
genesTablesd4myc <- cond.d4myc.genes %>% mutate(rank = rank(cond.d4myc.genes$auc,  ties.method = "random"))

ranksd4 <- deframe(cond.d4.genes)
ranksd4myc <- deframe(cond.d4myc.genes)

head(ranksd4)
head(ranksd4myc)

write.csv(ranksd4, "./Tables/GSEA_Fig6/ranks.Dll4vsCtrl.csv")

write.csv(ranksd4myc, "./Tables/GSEA_Fig6/ranks.Dll4MycvsCtrl.csv")

# Full one
E2F <- as.list(fgsea_sets$HALLMARK_E2F_TARGETS)

fgseaResD4<- fgsea(fgsea_sets, stats = ranksd4)

sort(ranksd4)
sort(unlist(ranksd4)) 
sort.int(ranksd4, decreasing = T)

ranksd4df <- as.data.frame(ranksd4)

ranksd4stat <- rownames(ranksd4df)

ranksd4stat1 <- ranksd4df$ranksd4

ranksd4stat2 <- ranksd4df$ranksd4[c(1:20, 30:40)]

class(ranksd4stat2)

View(ranksd4stat)
calcGseaStat(ranksd4stat1, 
             selectedStats = ranksd4stat2,
             returnLeadingEdge = T)

fgseaResD4myc<- fgsea(fgsea_sets, stats = ranksd4myc)
 
# Tidy Data with auc

fgseaResD4Tidy <- fgseaResD4 %>%
  as_tibble() %>%
  arrange(desc(NES))

fgseaResD4MycTidy <- fgseaResD4myc %>%
  as_tibble() %>%
  arrange(desc(NES))

fgseaResD4Tidy %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  head()

fgseaResD4MycTidy %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  head()

write.xlsx(fgseaResD4Tidy, "Tables/GSEA_Fig6/fgseaResD4Tidy.auc.xlsx")
write.xlsx(fgseaResD4MycTidy, "Tables/GSEA_Fig6/fgseaResD4MycTidy.auc.xlsx")

cond.fgsea <- read.xlsx("./Tables/GSEA_Fig6/fgseaResD4&D4MycTidy.auc.xlsx", sheet = 3, sep.names = " ", rowNames = T)

liver.gs <- subset(Liver, subset = Condition == "Dll4KO" | Condition == "Dll4/MycKO")

Idents(liver.gs) <- "Condition"

avgexp <- AverageExpression(liver.gs, return.seurat = T)

avgexp@meta.data$Condition <- avgexp@active.ident

ES.seurat <- enrichIt(obj = avgexp, gene.sets = GS.hallmark.HS, groups = 2, cores = 2)

avgexp <- Seurat::AddMetaData(avgexp, cond.fgsea)

avgexp@meta.data$Condition <- avgexp@active.ident

levels(avgexp) <- c("Dll4KO", "Dll4/MycKO")

avgexp@meta.data$Condition <- avgexp@active.ident

avgexp@meta.data$Condition <- as.factor(avgexp@meta.data$Condition)

# levels(avgexp) <- c("Dll4/MycKO", "Dll4KO")

avgexp@meta.data$Condition <- avgexp@active.ident

cairo_pdf("Plots/Figure_6/Figure_6I_GSEA_Heatmap_Dll4vsControl_Dll4MyciDECvsControl_unscaled.pdf", width = 5, height = 12, family = "Arial")
dittoHeatmap(avgexp, genes = NULL, metas = names(ES.seurat), 
             annot.by = "Condition",
             # order.by = c(2,1),
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


cairo_pdf("Plots/Figure_6/Figure_6I_GSEA_Heatmap_selected_Dll4vsControl_Dll4MyciDECvsControl_unscaled.pdf", width = 5, height = 6, family = "Arial")
dittoHeatmap(avgexp, genes = NULL, metas = names(ES.seurat.selected), 
             annot.by = "Condition",
             order.by = c(2,1),
             annot.colors = palette_Maca_GSEA,
             fontsize = 7, 
             cluster_cols = F,
             cluster_rows = F,
             heatmap.colors = BuRd,
             # heatmap.colors.max.scaled = BuRd,
             scaled.to.max = F,
             scale = "none",
             complex = T, data.out = F
             # legend_breaks = c(1, 0.8, 0.6, 0.4, 0.2, 0),
             )
dev.off()

cairo_pdf("Plots/Figure_6/Figure_6I_GSEA_Heatmap_Dll4vsControl_Dll4MyciDECvsControl_scaled.pdf", width = 5, height = 12, family = "Arial")
dittoHeatmap(avgexp, genes = NULL, metas = names(ES.seurat), 
             annot.by = "Condition",
             order.by = c(2,1),
             annot.colors = palette_Maca_GSEA,
             fontsize = 7, 
             cluster_cols = F,
             cluster_rows = F,
             heatmap.colors = BuRd,
             # heatmap.colors.max.scaled = BuRd,
             scaled.to.max = F,
             scale = "row",
             complex = T, data.out = F,
             # legend_breaks = c(1, 0.8, 0.6, 0.4, 0.2, 0),
             main = expression(Control~italic(Dll4)^iDEC~italic(Notch1)^iDEC~italic(Rbpj)^iDEC))
dev.off()

cairo_pdf("Plots/Figure_6/Figure_6I_GSEA_Heatmap_selected_Dll4vsControl_Dll4MyciDECvsControl_unscaled.pdf", width = 5, height = 6, family = "Arial")
dittoHeatmap(avgexp, genes = NULL, metas = names(ES.seurat.selected), 
             annot.by = "Condition",
             annot.colors = palette_Maca_GSEA,
             fontsize = 7, 
             cluster_cols = F,
             cluster_rows = F,
             heatmap.colors = BuRd,
             # heatmap.colors.max.scaled = BuRd,
             scaled.to.max = F,
             scale = "none",
             complex = T, data.out = F,
             # legend_breaks = c(1, 0.8, 0.6, 0.4, 0.2, 0),
             )
dev.off()

###########################################################################

# Fig6I.2 Dotplot of key Gene Set genes

# load ranks and Gene Sets

ranksd4 <- read.xlsx("./Tables/GSEA_Fig6/ranks.Dll4vsCtrl.xlsx", check.names = T, rowNames = F, rows = 1:32263)

ranksd4

h_gene_sets = msigdbr(species = "mouse", category = "H")

head(m_df)

fgsea_sets<- h_gene_sets %>% split(x = .$gene_symbol, f = .$gs_name)

class(fgsea_sets)

E2F_genes <- fgsea_sets$HALLMARK_E2F_TARGETS
Myc1_genes <- fgsea_sets$HALLMARK_MYC_TARGETS_V1
G2M_genes <- fgsea_sets$HALLMARK_G2M_CHECKPOINT
OxPhos_genes <- fgsea_sets$HALLMARK_OXIDATIVE_PHOSPHORYLATION

GeneSet_genes <- list(E2F = E2F_genes, Myc = Myc1_genes, G2M = G2M_genes, OxPhos = OxPhos_genes)

ranksd4 %>% select(matches(E2F_genes))

GeneSet_genes[1]

E2F_top10 <- ranksd4[ranksd4$gene %in% GeneSet_genes$E2F, ] %>% top_n(10) %>% extract2(1) %>% as.character()

E2F_top10

Myc_top10 <- ranksd4[ranksd4$gene %in% GeneSet_genes$Myc, ] %>% top_n(10) %>% extract2(1) %>% as.character()

Myc_top10

G2M_top10 <- ranksd4[ranksd4$gene %in% GeneSet_genes$G2M, ] %>% top_n(10) %>% extract2(1) %>% as.character()

G2M_top10

OxPhos_top10 <- ranksd4[ranksd4$gene %in% GeneSet_genes$OxPhos, ] %>% top_n(10) %>% extract2(1) %>% as.character()

C4_top10 <- c("Myc", "Odc1", "Vegfa", "Kcne3", "Esm1", "Apln", "Angpt2", "Cd34", "Fabp5", "Mif")

C5_top10 <- c("Top2a", "Stmn1", "Mki67", "Hmgb2", "Cenpf", "Pclaf", "Ube2c", "Hist1h2ap", "Hist1h2ae", "Prc1")

DotPlot_list <- c(E2F_top10, Myc_top10, G2M_top10, OxPhos_top10, C4_top10, C5_top10)

class(DotPlot_list)

DotPlot_list <- unique(DotPlot_list)

DotPlot_list_caps <- toupper(DotPlot_list) %>% unique()

DotPlot_list_caps

Idents(Liver) <- "Condition"

levels(Liver) <- c("Dll4/MycKO", "Dll4KO", "Control")

cairo_pdf("Plots/Figure_6/Figure_6I.2_DotPlot_GeneSets.pdf", width = 13, height = 2, family = "Arial")
DotPlot(Liver, features = DotPlot_list_caps, col.min = 0, dot.scale = 5, scale = T)+
  theme(strip.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"),
        legend.title = element_text(size = 10), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"), 
        plot.margin = unit(c(0, 0, 0, 0), "cm")
  ) +
  scale_y_discrete(labels=c(bquote(italic(Dll4/Myc)^iDEC), bquote(italic(Dll4)^iDEC), "Control"))+
  scale_x_discrete(label = DotPlot_list)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  geom_vline(xintercept = c(10.5, 21.5, 29.5, 39.5, 47.5),linetype = 2 )+
  theme(legend.box = "horizontal")
dev.off()

