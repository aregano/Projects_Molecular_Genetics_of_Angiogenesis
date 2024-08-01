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
library(msigdbr)

# rds

Liver <- readRDS("//Tierra/SC/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Liver_MacarenaFC/Figures_Paper/rds/Liver_20Jan21.Figures.All.Conditions.OldDll4Myc.rds")

table(Liver@meta.data$Condition)

Liver <- subset(Liver, subset = Condition == "Control" | Condition == "Dll4KO" | Condition == "Dll4/MycKO" | Condition == "D4KO_aVEGF")

Idents(Liver) <- "Condition"

levels(Liver) <- c("Control", "Dll4KO", "Dll4/MycKO", "D4KO_aVEGF")

Liver@active.ident -> Liver@meta.data$Condition

table(Liver@meta.data$Condition)

Liver4d <- readRDS("S:/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Liver_MacarenaFC/Figures_Paper/rds/Liver4d.Mapping.rds")

table(Liver4d@meta.data$Condition)

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

# Aternative Dotplot for Specific Genes 7E.3

# Get list of all top3 genes of each Condition, 27 total

E2F_genes_Condition <- read_excel_allsheets("Tables/GSEA_Fig7/E2F_ranks.Dll4KO.Dll4MycKO.Dll4KOaVEGF.xlsx")
Myc_genes_Condition <- read_excel_allsheets("Tables/GSEA_Fig7/Myc_ranks.Dll4KO.Dll4MycKO.Dll4KOaVEGF.xlsx")
OxPhos_genes_Condition <- read_excel_allsheets("Tables/GSEA_Fig7/OxPhos_ranks.Dll4KO.Dll4MycKO.Dll4KOaVEGF.xlsx")

Top3genes_GeneSets_Conditions <- c(E2F_genes_Condition[[1]][[1]][1:3], E2F_genes_Condition[[2]][[1]][c(2:4)], E2F_genes_Condition[[3]][[1]][c(2,3,5)],
                                   Myc_genes_Condition[[1]][[1]][1:3], Myc_genes_Condition[[2]][[1]][c(2, 4, 5)], Myc_genes_Condition[[3]][[1]][c(3,8,10)],
                                   OxPhos_genes_Condition[[1]][[1]][1:3], OxPhos_genes_Condition[[2]][[1]][c(1,5,8)], OxPhos_genes_Condition[[3]][[1]][c(2,3,5)])

unique(Top3genes_GeneSets_Conditions)
Top3genes_GeneSets_Conditions

C4_Tip_Cells <- c("Vegfa", "Kcne3", "Esm1", "Apln")
Cp_G2M <- c("Top2a", "Mki67", "Hmgb2", "Cenpf")
C1a_Arterial_capillaries <- c("Ntn4", "Ltbp4", "Msr1", "Ehd3", "Efnb2")

DotPlot_list <- c(Top3genes_GeneSets_Conditions, C4_Tip_Cells, Cp_G2M, C1a_Arterial_capillaries)

class(DotPlot_list)

DotPlot_list <- unique(DotPlot_list)

DotPlot_list_caps <- toupper(DotPlot_list)

DotPlot_list_caps

Idents(Liver) <- "Condition"

levels(Liver) <- c("D4KO_aVEGF", "Dll4/MycKO", "Dll4KO", "Control")

cairo_pdf("Plots/Figure_7/Figure_7E.3_DotPlot_GeneSets.pdf", width = 11, height = 2.2, family = "Arial")
DotPlot(Liver, features = DotPlot_list_caps, col.min = 0, dot.scale = 5, scale = T)+
  theme(strip.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"),
        legend.title = element_text(size = 10), legend.text  = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"), 
        plot.margin = unit(c(0, 0, 0, 0), "cm")
  ) +
  scale_y_discrete(labels=c(bquote(italic(Dll4)^iDEC+aVEGF), bquote(italic(Dll4/Myc)^iDEC), bquote(italic(Dll4)^iDEC), "Control"))+
  scale_x_discrete(label = DotPlot_list)+
  scale_colour_gradientn(colours = Bestholtz_palette)+
  geom_vline(xintercept = c(9.5, 18.5, 27.5, 31.5, 35.5),linetype = 2 )+
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

# Supplementary Featureplots

# Tip Cells + Proliferating

TipCell_Proliferating <- c("KCNE3", "ESM1")
TipCell_Proliferating <- as.data.frame(TipCell_Proliferating)
TipCell_Proliferating_text <- c("Kcne3", "Esm1")

myplots <- vector('list', nrow(TipCell_Proliferating))

for (i in 1:nrow(TipCell_Proliferating)) {
  
  p21 <- FeaturePlot(Liver4d, reduction = "ref.umap", pt.size = 1.3, features = TipCell_Proliferating[i, 1], order = T, combine = F)
  
  p22 <- FeaturePlot(Liver4d, reduction = "ref.umap", split.by = "Condition", order = T,features = TipCell_Proliferating[i, 1], combine = F)
  
  p21 <- lapply(X = p21, FUN = function(p) p + NoAxes()+ scale_colour_gradientn(colors = Bestholtz_palette))
  
  p22 <- lapply(X = p22, FUN = function(p) p + NoLegend() + NoAxes()+ scale_colour_gradientn(colors = Bestholtz_palette))
  
  p21 <- Reduce( `+`, p21 )+patchwork::plot_layout( ncol = 1 )
  
  p22 <- Reduce( `+`, p22 )+patchwork::plot_layout( ncol = 3 )
  
  p23 <- list(p21, p22)
  
  
  design <- c(patchwork::area(1, 1, 1, 1), patchwork::area(1, 2, 1, 4))
  
  p24 <- Reduce( `+`,  p23)+patchwork::plot_layout(design = design)
  
  
  myplots[[i]] <- local({
    i <- i
    p24
  })
  
  
}

p25 <- patchwork::wrap_plots(myplots, ncol = 1)

#myplots[1:nrow(genes)] + patchwork::plot_layout(byrow = T, widths = 25, heights = 5)
p25


cairo_pdf("Plots/Figure_7/Suplementary/Figure_S7.1_FeaturePlot.TipCell.pdf", width = 16, height = 6, family = "Arial")
p25
dev.off()

jpeg("Plots/Figure_7/Suplementary/Figure_S7.1_FeaturePlot.TipCell.jpeg", width = 16, height = 6, units = 'in', res = 800)
p25
dev.off()


##################################################################################

# Figure Supplementary 7.2. Violin Plot

genes <- list(c("KCNE3", "ESM1"))
genes <- as.data.frame(genes)
gene.names <- list(c("Kcne3", "Esm1"))
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
legend <- get_legend(plegend + theme(legend.position = "bottom", legend.justification = "center", legend.direction = "horizontal", legend.text = element_text(size = 18)))

p25 <- patchwork::wrap_plots(myplots, ncol = 1)

p25

p26 <- plot_grid(p25, legend, ncol = 1, rel_heights = c(1, .1))

p26

x.axis <- ggdraw() + draw_label("Expression Level", fontface='bold', angle = 90, size = 24)

p27 <- plot_grid(x.axis, p26, ncol = 2, rel_widths = c(.05, 1))

p27

cairo_pdf("Plots/Figure_7/Supplementary/Figure_S7.2_VlnPlot.pdf", width = 6, height = 6, family = "Arial")
p27
dev.off()

jpeg("Plots/Figure_7/Supplementary/Figure_S7.2_VlnPlot.jpeg", width = 6, height = 6, units = 'in', res = 800)
p27
dev.off()


####################################################################

# Supplementary 7.3 Ctrl 2w, Dll4KO4d, Dll4KO2w

Liver <- readRDS("//Tierra/SC/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Liver_MacarenaFC/Figures_Paper/rds/Liver_20Jan21.Figures.All.Conditions.ppt_Rui.rds")

levels(Liver@meta.data$Condition)

Liver <- subset(Liver, subset = Condition == "Control" | Condition == "Dll4KO(4d)+Vehicle" | Condition == "Dll4KO")

levels(Liver) <- c("Control", "Dll4KO(4d)+Vehicle", "Dll4KO")

Liver@meta.data$Condition <- Liver@active.ident

DimPlot(Liver)

p1 <- DimPlot(Liver, cols = my_palette_Rui_colors_B, pt.size = 1.2, split.by = "Condition", group.by = "FigClustering", combine = T)


new_labels <- c("Control" = "Control", "Dll4KO(4d)+Vehicle" = "italic(Dll4)^iDEC4d+Vehicle", "Dll4KO" = "italic(Dll4)^iDEC2w")
p1

# jpeg("Plots/Figure_4/Figure_4C_Axes_Conditions.jpeg", width = 28, height = 6, units = 'in', res = 800)
# p1+facet_grid(~Condition, labeller = as_labeller(new_labels, default = label_parsed))+theme(strip.text.x = element_text(size = 24), legend.text = element_text(size = 14))
# dev.off()

p21 <- p1+facet_grid(~Condition, labeller = as_labeller(new_labels, default = label_parsed))+
  theme(strip.text.x = element_text(size = 24), plot.title = element_blank())+NoLegend()+xlab("UMAP_1")+ylab("UMAP_2")

p22 <- DimPlot(Liver, cols = my_palette_Rui_colors_B, group.by = "FigClustering", pt.size = 1.2)+theme(legend.text = element_text(size = 14), axis.title.x = element_blank(), axis.line.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks = element_blank(), plot.title = element_blank())

p21
p22

p23 <- list(p21, p22)

design <- c(patchwork::area(1, 1, 1, 3), patchwork::area(1, 4, 1, 4.5))

p24 <- Reduce( `+`,  p23)+patchwork::plot_layout(design = design)


cairo_pdf("Plots/Figure_7/Supplementary/Figure_S7.3_UMAP.pdf",  width = 28, height = 6, family = "Arial")
p24
dev.off()


jpeg("Plots/Figure_7/Supplementary/Figure_S7.3_UMAP.jpeg", width = 28, height = 6, units = 'in', res = 800)
p24
dev.off()

###################################################################

# Supplementary S7.4 BarPlot

p1 <- dittoBarPlot(Liver, var = "FigClustering", group.by = "Condition", color.panel = my_palette_Rui_colors_B, xlab = NULL, main = NULL, x.reorder = c(1,3,2), y.breaks = c(0,0.25,0.5, 0.75, 1), min = 0, max = 1)+
  theme(axis.text=element_text(size=15), axis.title = element_text(size = 18), legend.text = element_text(size = 16), axis.text.x = element_text(angle = 25))+
  scale_x_discrete(labels=c("Control", bquote(italic(Dll4)^iDEC4d+Veh), bquote(italic(Dll4)^iDEC2w)))+
  geom_col(width = 0.1)+
  # scale_y_discrete(labels=c("0", "25", "50", "75", "100"))
  scale_y_continuous(labels = c("0", "25", "50", "75", "100"))
p1

cairo_pdf("Plots/Figure_7/Supplementary/Figure_S7.4_BarPlot.pdf",  width = 7, height = 6, family = "Arial")
p1
dev.off()


jpeg("Plots/Figure_7/Supplementary/Figure_S7.4_BarPlot.jpeg", width = 7, height = 6, units = 'in', res = 800)
p1
dev.off()

###################################################################

# Supplementary S7.5 ViolinPlot

palette_Vln_Dll4 <- c("#606060", "#C19EDD", "#F94040")

# Figure Supplementary 7.2. Violin Plot

genes <- list(c("KCNE3", "ESM1", "VEGFA", "APLN", "ODC1", "STMN1"))
genes <- as.data.frame(genes)
gene.names <- list(c("Kcne3", "Esm1", "Vegfa", "Apln", "Odc1", "Stmn1"))
gene.names <- as.data.frame(gene.names)

#  Script

myplots <- vector('list', nrow(genes))


for (i in 1:nrow(genes)) {
  
  p21 <- VlnPlot(Liver, features = genes[i, 1], group.by = "Condition", cols = palette_Vln_Dll4)+ NoLegend() + theme(text = element_text(family="Arial"), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.title.y=element_blank(), plot.title = element_text(face = "bold.italic"))+ggtitle(gene.names[i, 1])
  
  myplots[[i]] <- local({
    i <- i
    p21
  })
  
}

plegend <- VlnPlot(Liver, features = genes[i, 1], group.by = "Condition") +
  scale_fill_manual(values= palette_Vln_Dll4, labels=c("Control", expression(italic(Dll4)^"iDEC4d"), expression(italic(Dll4)^"iDEC2w")))
plegend
legend <- get_legend(plegend + theme(legend.position = "bottom", legend.justification = "center", legend.direction = "horizontal", legend.text = element_text(size = 18)))

p25 <- patchwork::wrap_plots(myplots, ncol = 1)

p25

p26 <- plot_grid(p25, legend, ncol = 1, rel_heights = c(1, .1))

p26

x.axis <- ggdraw() + draw_label("Expression Level", fontface='bold', angle = 90, size = 24)

p27 <- plot_grid(x.axis, p26, ncol = 2, rel_widths = c(.05, 1))

p27

cairo_pdf("Plots/Figure_7/Supplementary/Figure_S7.5_VlnPlot.pdf", width = 6, height = 18, family = "Arial")
p27
dev.off()

jpeg("Plots/Figure_7/Supplementary/Figure_S7.5_VlnPlot.jpeg", width = 6, height = 18, units = 'in', res = 800)
p27
dev.off()


#  Supplementary S7A Vln Plot

Liver <- subset(Liver, subset = Condition == "Control" | Condition == "Dll4KO" | Condition == "D4KO_aVEGF")

Idents(Liver) <- "Condition"

levels(Liver) <- c("Control", "Dll4KO", "D4KO_aVEGF")

Liver@active.ident -> Liver@meta.data$Condition

table(Liver@meta.data$Condition)

palette_Maca_Vln <- c("#606060", "#F94040", "#BB005E")

genes <- list(c("DLL4"))
genes <- as.data.frame(genes)
gene.names <- list(c("Dll4"))
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
  scale_fill_manual(values= palette_Maca_Vln, labels=c("Control", expression(italic("Dll4"^"iDEC")), expression(italic(Dll4)^"iDEC"+"aVEGF")))
plegend
legend <- get_legend(plegend + theme(legend.position = "bottom", legend.justification = "center", legend.direction = "horizontal", legend.text = element_text(size = 14)))

p25 <- patchwork::wrap_plots(myplots, ncol = 1)

p26 <- plot_grid(p25, legend, ncol = 1, rel_heights = c(1, .1))

x.axis <- ggdraw() + draw_label("Expression Level", fontface='bold', angle = 90, size = 12)

p27 <- plot_grid(x.axis, p26, ncol = 2, rel_widths = c(.07, 1))

p27

cairo_pdf("Plots/Figure_7/Supplementary/Figure_S7A_VlnPlot.pdf", width = 5, height = 3, family = "Arial")
p27
dev.off()

jpeg("Plots/Figure_7/Supplementary/Figure_S7A_VlnPlot.jpeg", width = 5, height = 3, units = 'in', res = 800)
p27
dev.off()

###############################################################################################################################


#  Supplementary S7C GSEA Biological Processes

#  Supplementary S7D Vln Plot



palette_Maca_Vln <- c("#606060", "#F94040", "#CCF9E8")

genes <- list(c("DLL4"))
genes <- as.data.frame(genes)
gene.names <- list(c("Dll4"))
gene.names <- as.data.frame(gene.names)

#  Script

myplots <- vector('list', nrow(genes))

table(Liver4d@meta.data$Condition)

for (i in 1:nrow(genes)) {
  
  p21 <- VlnPlot(Liver4d, features = genes[i, 1], group.by = "Condition", cols = palette_Maca_Vln)+ NoLegend() + theme(text = element_text(family="Arial"), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.title.y=element_blank(), plot.title = element_text(face = "bold.italic"))+ggtitle(gene.names[i, 1])
  
  myplots[[i]] <- local({
    i <- i
    p21
  })
  
}

plegend <- VlnPlot(Liver4d, features = genes[i, 1], group.by = "Condition") +
  scale_fill_manual(values= palette_Maca_Vln, labels=c("Control", expression(italic("Dll4")^"iDEC(4d)"), expression(italic(Dll4)^"iDEC(4d)"+"SL327")))
plegend
legend <- get_legend(plegend + theme(legend.position = "bottom", legend.justification = "center", legend.direction = "horizontal", legend.text = element_text(size = 14)))

p25 <- patchwork::wrap_plots(myplots, ncol = 1)

p26 <- plot_grid(p25, legend, ncol = 1, rel_heights = c(1, .1))

x.axis <- ggdraw() + draw_label("Expression Level", fontface='bold', angle = 90, size = 12)

p27 <- plot_grid(x.axis, p26, ncol = 2, rel_widths = c(.07, 1))

p27

cairo_pdf("Plots/Figure_7/Supplementary/Figure_S7D_VlnPlot.pdf", width = 5, height = 3, family = "Arial")
p27
dev.off()

jpeg("Plots/Figure_7/Supplementary/Figure_S7D_VlnPlot.jpeg", width = 5, height = 3, units = 'in', res = 800)
p27
dev.off()

#########################################################################


#  Suplementary S7C GSEA Biological Processes

Liver <- readRDS("//Tierra/SC/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Liver_MacarenaFC/Figures_Paper/rds/Liver_20Jan21.Figures.All.Conditions.OldDll4Myc.rds")
table(Liver@meta.data$Condition)

LiverD4KO <- subset(Liver, subset = Condition == "Control" | Condition == "Dll4KO")
LiverD4aVEGF <- subset(Liver, subset = Condition == "Control" | Condition == "D4KO_aVEGF")

LiverD4KO@meta.data$Condition <- LiverD4KO@active.ident

Idents(LiverD4aVEGF) <- "Condition"
LiverD4aVEGF@meta.data$Condition <- LiverD4aVEGF@active.ident




# we have all the genes for each cluster
dplyr::count(liverd4.genes, group)
dplyr::count(liverd4avegf.genes, group)

msigdbr_show_species()

m_df<- msigdbr(species = "Mus musculus", category = 13)

h_df<- msigdbr(species = "Homo sapiens", category = 7)

table(m_df$gs_subcat)

table(h_gene_sets$gs_cat)

# Using GO:BP subcategory genes only

h_gene_sets = msigdbr(species = "mouse", subcategory = "GO:BP")

fgsea_sets<- h_gene_sets %>% split(x = .$gene_symbol, f = .$gs_name) %>% list.filter("Type" > 50)

fgsea_sets$GOBP_AGING

liver.genes <- wilcoxauc(Liver, 'Condition')
liverd4.genes <- wilcoxauc(LiverD4KO, 'Condition')
liverd4avegf.genes <- wilcoxauc(LiverD4aVEGF, 'Condition')


# Myc Targets


liverd4.genes %>% 
  dplyr::filter(group == "Dll4KO") %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 20)

liverd4avegf.genes %>%
  dplyr::filter(group == "D4KO_aVEGF") %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 20)

liverd4.sel.genes<- liverd4.genes %>%
  dplyr::filter(group == "Dll4KO") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)

liverd4avegf.sel.genes<- liverd4avegf.genes %>%
  dplyr::filter(group == "D4KO_aVEGF") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)

# select only the feature and auc columns for fgsea, which statistics to use is an open question

liverd4.genes$feature <- tolower(liverd4.genes$feature) %>% str_to_title()

liverd4avegf.genes$feature <- tolower(liverd4avegf.genes$feature) %>% str_to_title()


cond.d4.genes <- liverd4.genes %>%
  dplyr::filter(group == "Dll4KO") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)

cond.d4avegf.genes <- liverd4avegf.genes %>%
  dplyr::filter(group == "D4KO_aVEGF") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc) 

genesTablesd4 <- cond.d4.genes %>% mutate(rank = rank(cond.d4.genes$auc,  ties.method = "random")) 
genesTablesd4avegf <- cond.d4avegf.genes %>% mutate(rank = rank(cond.d4avegf.genes$auc,  ties.method = "random"))

ranksd4 <- deframe(cond.d4.genes)
ranksd4avegf <- deframe(cond.d4avegf.genes)

head(ranksd4)
head(ranksd4avegf)

write.csv(ranksd4, "./Tables/GSEA_Fig7/ranks.Dll4vsCtrl.BP.csv")

write.csv(ranksd4avegf, "./Tables/GSEA_Fig7/ranks.Dll4aVEGFvsCtrl.BP.csv")

# Full one

fgseaResD4<- fgsea(fgsea_sets, stats = ranksd4)

fgseaResD4aVEGF<- fgsea(fgsea_sets, stats = ranksd4avegf)

myfgsea <- list(fgseaResD4, fgseaResD4aVEGF)

condition_short <- c("Dll4KO", "D4KO_aVEGF")

names(myfgsea) <- condition_short

write.xlsx(myfgsea, file = "./Tables/fgsea_Sup7C_Condition_Dll4_v_Control_Dll4aVEGF_v_Control_GO_BP.xlsx")

myfgsea.filter <- NULL

for (i in 1:2) {
  myfgsea.filter[[i]] <- local({
    i <- i
    myfgsea[[condition_short[i]]][size > 50 & padj < 0.05]
  })
}



names(myfgsea.filter) <- condition_short

write.xlsx(myfgsea.filter, file = "./Tables/GSEA_Fig7/fgsea_Sup7C_Condition_Dll4_v_Control_Dll4aVEGF_v_Control_GO_BP.filter.xlsx")

########################Heatmap##############################

# Traspose this data!

cond.fgsea <- read.xlsx("./Tables/GSEA_Fig7/fgsea_Sup7C_Condition_Dll4_v_Control_Dll4Myc_Control_Dll4aVEGF_v_Control_GO_BP.filter.xlsx", sheet = 4, sep.names = " ", rowNames = F)

liver.gs <- subset(Liver, subset = Condition == "Dll4KO" | Condition == "Dll4/MycKO" | Condition == "D4KO_aVEGF")

########################Dll4##########################

dll4 <- gsea_liver[["Dll4KO"]]

dll4t <- setNames(data.frame(t(dll4[ , - 1])), dll4[ , 1])

rownames(dll4t)

colnames(dll4t)

selected.GOBP <- read.xlsx("Tables/GSEA_Fig7/GOBP_genes.xlsx", sheet = 2)

colnames(selected.GOBP)

dll4.selected <- dll4t[, (colnames(dll4t) %in% colnames(selected.GOBP))]

#####################Dll4/Myc#######################

dll4myc <- gsea_liver[["Dll4MycKO"]]

dll4myct <- setNames(data.frame(t(dll4myc[ , - 1])), dll4myc[ , 1])

rownames(dll4myct)

colnames(dll4myct)

dll4myc.selected <- dll4myct[, (colnames(dll4myct) %in% colnames(selected.GOBP))]

#####################Dll4_aVEGF########################

dll4aVEGF <- gsea_liver[["D4KO_aVEGF"]]

dll4aVEGFt <- setNames(data.frame(t(dll4aVEGF[ , - 1])), dll4aVEGF[ , 1])

rownames(dll4aVEGFt)

colnames(dll4aVEGFt)

dll4aVEGF.selected <- dll4aVEGFt[, (colnames(dll4aVEGFt) %in% colnames(selected.GOBP))]


#################List################

class(dll4.selected)

liver.gs <- subset(Liver, subset = Condition == "Dll4KO" | Condition == "Dll4/MycKO" | Condition == "D4KO_aVEGF")

GOBP <- list(dll4.selected, dll4myc.selected, dll4aVEGF.selected)

condition_short <- c("Dll4KO", "Dll4/MycKO", "D4KO_aVEGF")

names(GOBP) <- condition_short

write.xlsx(GOBP, "Tables/GSEA_Fig7/GOBP_Macarena_selected.xlsx")

Idents(liver.gs) <- "Condition"

avgexp <- AverageExpression(liver.gs, return.seurat = T)

avgexp@meta.data$Condition <- avgexp@active.ident

selected.GOBP <- read.xlsx("Tables/GSEA_Fig7/GOBP_Macarena_selected.xlsx", sheet = 4, rowNames = T)

avgexp <- Seurat::AddMetaData(avgexp, selected.GOBP)

avgexp@meta.data$Condition <- avgexp@active.ident

levels(avgexp) <- c("Dll4KO", "Dll4/MycKO", "D4KO_aVEGF")

avgexp@meta.data$Condition <- avgexp@active.ident

cairo_pdf("Plots/Figure_7/Supplementary/Figure_S7C_GSEA_GOBP_Dll4vsControl_Dll4MyciDECvsControl_Dll4KOaVEGFvControl_unscaled.pdf", width = 5, height = 12, family = "Arial")
dittoHeatmap(avgexp, genes = NULL, metas = names(selected.GOBP), 
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


#################################################################################################################################


# Figure. Violin Plot Adm Odc1

genes <- list(c("ADM", "ODC1"))
genes <- as.data.frame(genes)
gene.names <- list(c("Adm", "Odc1"))
gene.names <- as.data.frame(gene.names)

palette_Maca_Vln <- c("#606060", "#F94040", "#C19EDD", "#BB005E")

levels(Liver) <- c("Control", "Dll4KO", "Dll4/MycKO", "D4KO_aVEGF")

Liver@active.ident -> Liver@meta.data$Condition

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
  scale_fill_manual(values= palette_Maca_Vln, labels=c("Control", expression(italic(Dll4)^"iDEC"), expression(italic(Dll4/Myc)^"iDEC"), expression(italic(Dll4)^"iDEC"+aVEGF)))
plegend
legend <- get_legend(plegend + theme(legend.position = "bottom", legend.justification = "center", legend.direction = "horizontal", legend.text = element_text(size = 7), legend.key.size = unit(0.2, "cm")))

p25 <- patchwork::wrap_plots(myplots, ncol = 1)

p26 <- plot_grid(p25, legend, ncol = 1, rel_heights = c(1, .05))

x.axis <- ggdraw() + draw_label("Expression Level", fontface='bold', angle = 90, size = 20)

p27 <- plot_grid(x.axis, p26, ncol = 2, rel_widths = c(.1, 1))

p27

cairo_pdf("Plots/Figure_7/Figure_7_VlnPlot_odc1_Adm.pdf", width = 3.5, height = 6, family = "Arial")
p27
dev.off()

jpeg("Plots/Figure_7/Figure_7_VlnPlot_odc1_Adm.jpeg", width = 3.5, height = 6, units = 'in', res = 800)
p27
dev.off()

########################################################

# Figure S7D. Violin Plot

genes <- list(c("HES1", "MSR1", "LTBP4", "EFNB2"))
genes <- as.data.frame(genes)
gene.names <- list(c("Hes1", "Msr1", "Ltbp4", "Efnb2"))
gene.names <- as.data.frame(gene.names)

Liver <- subset(Liver, subset = Condition == "Control" | Condition == "Dll4KO" | Condition == "D4KO_aVEGF")

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

cairo_pdf("Plots/Figure_7/Supplementary/Figure_S7D_VlnPlot.pdf", width = 3.5, height = 12, family = "Arial")
p27
dev.off()

jpeg("Plots/Figure_7/Supplementary/Figure_S7D_VlnPlot.jpeg", width = 3.5, height = 12, units = 'in', res = 800)
p27
dev.off()
