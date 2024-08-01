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


# rds

Liver <- readRDS("//Tierra/SC/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Liver_MacarenaFC/Figures_Paper/rds/Liver_20Jan21.Figures.All.Conditions.OldDll4Myc.rds")

table(Liver@meta.data$Condition)

Liver <- subset(Liver, subset = Condition == "Control" | Condition == "Dll4KO" | Condition == "Dll4/MycKO" | Condition == "D4KO_aVEGF")

Idents(Liver) <- "Condition"

levels(Liver) <- c("Control", "Dll4KO", "Dll4/MycKO", "D4KO_aVEGF")

Liver@active.ident -> Liver@meta.data$Condition

Liver4d <- readRDS("S:/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Liver_MacarenaFC/Figures_Paper/rds/Liver4d.Mapping.rds")



# Color Palettes

my_palette_Rui_colors <- c("#E95A74", "#50B6EF", "#45FF8E", "#F4A753", "#A80519", "#880088", "#E28CF4", "#C1B80C","#FC0808",  "#0E47D8")

my_palette_Rui_colors_B <- c("#F4A753", "#E95A74", "#50B6EF", "#FC0808",  "#0E47D8", "#45FF8E", "#A80519", "#E28CF4", "#880088", "#C1B80C")

custom_rdylblu_palette <- c("#005AD9", "#0070FF", "#FFE135", "#FFDA13", "#E60026", "#C60012")

macarena_palette <- c("#FDDC85", "#6E2914")

Bestholtz_palette <- c("#DEDAD6", "#FEE392", "#FEC44E", "#FE9929", "#ED6F11", "#CC4C17", "#993411", "#65260C")

palette_Maca_Vln <- c("#606060", "#F94040", "#BB005E")

BuRd <- c("#052F60", "#2066AC", "#4392C3", "#92C5DE", "#D1E4EF", "#F7F7F7", "#FCDBC6", "#F4A482", "#D55F4D", "#B1172B", "#67001F")


palette_Maca_GSEA <- c("#F94040", "mediumpurple1", "#BB005E")

###############################################################################################################################

## add the Arial font
font_add("Arial", regular = "arial.ttf",
         bold = "arialbd.ttf", italic = "ariali.ttf", bolditalic = "arialbi.ttf")


###############################################################################################################################

# Figure 7D.1 UMAP

DimPlot(Liver, label = T)

Idents(Liver) <- Liver@meta.data$FigClustering

p1 <- DimPlot(Liver, cols = my_palette_Rui_colors_B, pt.size = 1.2, group.by = "FigClustering", split.by = "Condition", combine = T)

p1

new_labels <- c("Control" = "Control", "Dll4KO" = "italic(Dll4)^iDEC", "D4KO_aVEGF" = "italic(Dll4)^iDEC+aVEGF")

p21 <- p1+facet_grid(~Condition, labeller = as_labeller(new_labels, default = label_parsed))+
  theme(strip.text.x = element_text(size = 24), plot.title = element_blank())+NoLegend()+xlab("UMAP_1")+ylab("UMAP_2")

p22 <- DimPlot(Liver, group.by = "FigClustering", cols = my_palette_Rui_colors_B, pt.size = 1.2)+theme(legend.text = element_text(size = 14), axis.title.x = element_blank(), axis.line.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks = element_blank(), plot.title = element_blank())

p21
p22

p23 <- list(p21, p22)

design <- c(patchwork::area(1, 1, 1, 3), patchwork::area(1, 4, 1, 4.5))

p24 <- Reduce( `+`,  p23)+patchwork::plot_layout(design = design)


cairo_pdf("Plots/Figure_7/Figure_7D.1_UMAP.pdf",  width = 29, height = 6, family = "Arial")
p24
dev.off()


jpeg("Plots/Figure_7/Figure_7D.1_UMAP.jpeg", width = 29, height = 6, units = 'in', res = 800)
p24
dev.off()



###############################################################################################################################


# Figure 7D.2 Bar Plot

table(Liver@meta.data$Condition)

p1 <- dittoBarPlot(Liver, var = "FigClustering", group.by = "Condition", x.reorder = c(1,3,2), color.panel = my_palette_Rui_colors_B, xlab = NULL, main = NULL, y.breaks = c(0,0.25,0.5, 0.75, 1), min = 0, max = 1)+
  theme(axis.text=element_text(size=12), axis.title = element_text(size = 18), legend.text = element_text(size = 16), axis.text.x = element_text(angle = 25))+
  scale_x_discrete(labels=c("Control", bquote(italic(Dll4)^iDEC), bquote(italic(Dll4)^iDEC+aVEGF)))+
  geom_col(width = 0.1)+
  # scale_y_discrete(labels=c("0", "25", "50", "75", "100"))
  scale_y_continuous(labels = c("0", "25", "50", "75", "100"))
p1

cairo_pdf("Plots/Figure_7/Figure_7D.2_BarPlot.pdf",  width = 6.5, height = 6, family = "Arial")
p1
dev.off()


jpeg("Plots/Figure_7/Figure_7D.2_BarPlot.jpeg", width = 6.5, height = 6, units = 'in', res = 800)
p1
dev.off()


###############################################################################################################################

# Figura 7E. GSEA Heatmap.


LiverD4KOVEGF <- subset(Liver, subset = Condition == "Control" | Condition == "D4KO_aVEGF")

LiverD4KOVEGF@meta.data$Condition <- LiverD4KOVEGF@active.ident


GS.hallmark.HS <- getGeneSets(library = "H")

liver.genes <- wilcoxauc(Liver, 'Condition')
liverd4avegf.genes <- wilcoxauc(LiverD4KOVEGF, 'Condition')

# we have all the genes for each cluster
dplyr::count(liverd4avegf.genes, group)

msigdbr_show_species()

m_df<- msigdbr(species = "Mus musculus", category = 13)

h_df<- msigdbr(species = "Homo sapiens", category = 7)

# Using Hallmark genes only

h_gene_sets = msigdbr(species = "mouse", category = "H")

head(m_df)

fgsea_sets<- h_gene_sets %>% split(x = .$gene_symbol, f = .$gs_name)

fgsea_sets$HALLMARK_MYC_TARGETS_V1


# Myc Targets


liverd4avegf.genes %>% 
  dplyr::filter(group == "D4KO_aVEGF") %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 20)

liverd4avegf.sel.genes<- liverd4avegf.genes %>%
  dplyr::filter(group == "D4KO_aVEGF") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)


# Good to see ODC1, KCNE3, PEG10, TPI1, MIF

# select only the feature and auc columns for fgsea, which statistics to use is an open question

liverd4avegf.genes$feature <- tolower(liverd4avegf.genes$feature) %>% str_to_title()

cond.d4avegf.genes <- liverd4avegf.genes %>%
  dplyr::filter(group == "D4KO_aVEGF") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)

genesTablesd4avegf <- cond.d4avegf.genes %>% mutate(rank = rank(cond.d4avegf.genes$auc,  ties.method = "random")) 

ranksd4avegf <- deframe(cond.d4avegf.genes)

head(ranksd4avegf)

ranksd4avegf.d <- as.data.frame(ranksd4avegf)

write.xlsx(ranksd4avegf.d, "Tables/GSEA_Fig7/ranksDll4KOaVEGF.auc.xlsx", rowNames = T)

# Full one

fgseaResD4aVEGF<- fgsea(fgsea_sets, stats = ranksd4avegf)

# Tidy Data with auc

fgseaResD4aVEGFTidy <- fgseaResD4aVEGF %>%
  as_tibble() %>%
  arrange(desc(NES))

fgseaResD4aVEGFTidy %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  head()


write.xlsx(fgseaResD4aVEGFTidy, "Tables/GSEA_Fig7/fgseaResD4aVEGFTidy.auc.xlsx")

cond.fgsea <- read.xlsx("./Tables/GSEA_Fig7/fgseaResD4&D4Myc&D4aVEGFTidy.auc.xlsx", sheet = 4, sep.names = " ", rowNames = T)

liver.gs <- subset(Liver, subset = Condition == "Dll4KO" | Condition == "Dll4/MycKO" | Condition == "D4KO_aVEGF")

Idents(liver.gs) <- "Condition"

avgexp <- AverageExpression(liver.gs, return.seurat = T)

avgexp@meta.data$Condition <- avgexp@active.ident

ES.seurat <- enrichIt(obj = avgexp, gene.sets = GS.hallmark.HS, groups = 3, cores = 2)


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
                         "HALLMARK_INFLAMMATORY_RESPONSE",
                         "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
)

ES.seurat.selected <- ES.seurat[, (colnames(ES.seurat) %in% selected.HallmarkGS)]

ES.seurat.selected <- ES.seurat.selected[,selected.HallmarkGS]

avgexp <- Seurat::AddMetaData(avgexp, cond.fgsea)

avgexp@meta.data$Condition <- avgexp@active.ident

levels(avgexp) <- c("Dll4KO", "Dll4/MycKO", "D4KO_aVEGF")

avgexp@meta.data$Condition <- avgexp@active.ident

cairo_pdf("Plots/Figure_7/Figure_7E_GSEA_Heatmap_Dll4vsControl_Dll4MyciDECvsControl_Dll4KOaVEGFvControl_unscaled.pdf", width = 5, height = 12, family = "Arial")
dittoHeatmap(avgexp, genes = NULL, metas = names(ES.seurat), 
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


cairo_pdf("Plots/Figure_7/Figure_7E_GSEA_Heatmap_selected_Dll4vsControl_Dll4MyciDECvsControl_Dll4KOaVEGFvControl_unscaled.pdf", width = 5, height = 6, family = "Arial")
dittoHeatmap(avgexp, genes = NULL, metas = names(ES.seurat.selected), 
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
             complex = T, data.out = F,
             # legend_breaks = c(1, 0.8, 0.6, 0.4, 0.2, 0),
             )
dev.off()

cairo_pdf("Plots/Figure_7/Figure_7E_GSEA_Heatmap_Dll4vsControl_Dll4MyciDECvsControl_Dll4KOaVEGFvControl_scaled.pdf", width = 5, height = 12, family = "Arial")
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
             scale = "row",
             complex = T, data.out = F)
             # legend_breaks = c(1, 0.8, 0.6, 0.4, 0.2, 0),
             # main = expression(Control~italic(Dll4)^iDEC~italic(Notch1)^iDEC~italic(Rbpj)^iDEC))
dev.off()

cairo_pdf("Plots/Figure_7/Figure_7E_GSEA_Heatmap_selected_Dll4vsControl_Dll4MyciDECvsControl_Dll4KOaVEGFvControl_scaled.pdf", width = 5, height = 6, family = "Arial")
dittoHeatmap(avgexp, genes = NULL, metas = names(ES.seurat.selected), 
             annot.by = "Condition",
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
             )
dev.off()

###############################################################################################################################

# Figure 7E.2 DotPlot GSEA

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

C2a_top10 <- c("Ly6a", "Ly6c1", "Atp13a3", "Ednrb", "Clu", "Adgrg6", "Cst3", "Efnb1", "Tm4sf1", "Fbln2")

C1a_top10 <- c("Ntn4", "Ltbp4", "Msr1", "Hes1", "Ehd3", "Serinc3", "Il6st", "Dll4", "Efnb2",  "Epas1")

DotPlot_list <- c(E2F_top10, Myc_top10, G2M_top10, OxPhos_top10, C4_top10, C5_top10, C2a_top10)

class(DotPlot_list)

DotPlot_list <- unique(DotPlot_list)

DotPlot_list_caps <- toupper(DotPlot_list)

DotPlot_list_caps

Idents(Liver) <- "Condition"

levels(Liver) <- c("D4KO_aVEGF", "Dll4/MycKO", "Dll4KO", "Control")

cairo_pdf("Plots/Figure_7/Figure_7E.2_DotPlot_GeneSets.pdf", width = 15, height = 2.2, family = "Arial")
DotPlot(Liver, features = DotPlot_list_caps, col.min = 0, dot.scale = 5, scale = T)+
  theme(strip.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"),
        legend.title = element_text(size = 10), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"), 
        plot.margin = unit(c(0, 0, 0, 0), "cm")
  ) +
  scale_y_discrete(labels=c(bquote(italic(Dll4)^iDEC+aVEGF), bquote(italic(Dll4/Myc)^iDEC), bquote(italic(Dll4)^iDEC), "Control"))+
  scale_x_discrete(label = DotPlot_list)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  geom_vline(xintercept = c(10.5, 21.5, 29.5, 39.5, 47.5, 56.5),linetype = 2 )+
  theme(legend.box = "horizontal")
dev.off()

cairo_pdf("Plots/Figure_7/Figure_7E.2a_DotPlot_GeneSets.pdf", width = 15, height = 2.2, family = "Arial")
DotPlot(Liver, features = DotPlot_list_caps, col.min = 0, dot.scale = 5, scale = T)+
  theme(strip.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"),
        legend.title = element_text(size = 10), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"), 
        plot.margin = unit(c(0, 0, 0, 0), "cm")
  ) +
  scale_y_discrete(labels=c(bquote(italic(Dll4)^iDEC+aVEGF), bquote(italic(Dll4/Myc)^iDEC), bquote(italic(Dll4)^iDEC), "Control"))+
  scale_x_discrete(label = DotPlot_list)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  geom_vline(xintercept = c(10.5, 21.5, 29.5, 39.5, 47.5, 56.5),linetype = 2 )+
  theme(legend.box = "horizontal")
dev.off()

#################################################################################################################################


# Figure 7F. Violin Plot

genes <- list(c("MYC", "ODC1"))
genes <- as.data.frame(genes)
gene.names <- list(c("Myc", "Odc1"))
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
  scale_fill_manual(values= palette_Maca_Vln, labels=c("Control", expression(italic(Dll4)^"iDEC"), expression(italic(Dll4)^"iDEC"+aVEGF)))
plegend
legend <- get_legend(plegend + theme(legend.position = "bottom", legend.justification = "center", legend.direction = "horizontal", legend.text = element_text(size = 10), legend.key.size = unit(0.3, "cm")))

p25 <- patchwork::wrap_plots(myplots, ncol = 1)

p26 <- plot_grid(p25, legend, ncol = 1, rel_heights = c(1, .07))

x.axis <- ggdraw() + draw_label("Expression Level", fontface='bold', angle = 90, size = 24)

p27 <- plot_grid(x.axis, p26, ncol = 2, rel_widths = c(.1, 1))

cairo_pdf("Plots/Figure_7/Figure_7F_VlnPlot.pdf", width = 3.5, height = 6, family = "Arial")
p27
dev.off()

jpeg("Plots/Figure_7/Figure_7F_VlnPlot.jpeg", width = 3.5, height = 6, units = 'in', res = 800)
p27
dev.off()



###############################################################################################################################

# Figure 7J.1 UMAP

DimPlot(Liver4d, label = T, reduction = "ref.umap")

Idents(Liver4d) <- Liver@meta.data$FigClustering

p1 <- DimPlot(Liver4d, cols = my_palette_Rui_colors_B, pt.size = 1.2, reduction = "ref.umap", group.by = "FigClustering", split.by = "Condition", combine = T)

p1

new_labels <- c("Control(4d)" = "Control", "Dll4KO(4d)+Vehicle" = "italic(Dll4)^iDEC+Veh", "Dll4KO(4d)+SL327" = "italic(Dll4)^iDEC+SL327")

p21 <- p1+facet_grid(~Condition, labeller = as_labeller(new_labels, default = label_parsed))+
  theme(strip.text.x = element_text(size = 24), plot.title = element_blank())+NoLegend()+xlab("UMAP_1")+ylab("UMAP_2")

p22 <- DimPlot(Liver4d, group.by = "FigClustering", reduction = "ref.umap", cols = my_palette_Rui_colors_B, pt.size = 1.2)+theme(legend.text = element_text(size = 14), axis.title.x = element_blank(), axis.line.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks = element_blank(), plot.title = element_blank())

p21
p22

p23 <- list(p21, p22)

design <- c(patchwork::area(1, 1, 1, 3), patchwork::area(1, 4, 1, 4.5))

p24 <- Reduce( `+`,  p23)+patchwork::plot_layout(design = design)


cairo_pdf("Plots/Figure_7/Figure_7J.1_UMAP.pdf",  width = 29, height = 6, family = "Arial")
p24
dev.off()


jpeg("Plots/Figure_7/Figure_7J.1_UMAP.jpeg", width = 29, height = 6, units = 'in', res = 800)
p24
dev.off()



###############################################################################################################################


# Figure 7J.2 Bar Plot

table(Liver@meta.data$Condition)

p1 <- dittoBarPlot(Liver4d, var = "FigClustering", group.by = "Condition", x.reorder = c(1,3,2), color.panel = my_palette_Rui_colors_B, xlab = NULL, main = NULL, y.breaks = c(0,0.25,0.5, 0.75, 1), min = 0, max = 1)+
  theme(axis.text=element_text(size=12), axis.title = element_text(size = 18), legend.text = element_text(size = 16), axis.text.x = element_text(angle = 25))+
  scale_x_discrete(labels=c("Control", bquote(italic(Dll4)^iDEC+Veh), bquote(italic(Dll4)^iDEC+SL327)))+
  geom_col(width = 0.1)+
  # scale_y_discrete(labels=c("0", "25", "50", "75", "100"))
  scale_y_continuous(labels = c("0", "25", "50", "75", "100"))
p1

cairo_pdf("Plots/Figure_7/Figure_7J.2_BarPlot.pdf",  width = 6.5, height = 6, family = "Arial")
p1
dev.off()


jpeg("Plots/Figure_7/Figure_7J.2_BarPlot.jpeg", width = 6.5, height = 6, units = 'in', res = 800)
p1
dev.off()


###############################################################################################################################
