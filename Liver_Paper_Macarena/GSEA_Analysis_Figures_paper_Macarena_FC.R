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

# Gene Set Database
# https://cran.r-project.org/web/packages/msigdbr/vignettes/msigdbr-intro.html

# Seurat
# https://crazyhottommy.github.io/scRNA-seq-workshop-Fall-2019/scRNAseq_workshop_3.html#downstream_analysis_of_scrnaseq_data 

# Escape
# http://www.bioconductor.org/packages/release/bioc/vignettes/escape/inst/doc/vignette.html

# ClusterProfiler
# https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/

# https://uclouvain-cbio.github.io/WSBIM2122/sec-gsea.html

##########################################################################################################################

# Seurat Tutorial Regular Plotting

# for GSEA, we need the information of all genes, Seurat is just too slow if we test
# all 20,000 genes. instead let's try presto which performs a fast Wilcoxon rank sum test 

# library(devtools)
# install_github('immunogenomics/presto')

Liver <- readRDS("S:/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Liver_MacarenaFC/Figures_Paper/rds/Liver_20Jan21.Figures.All.Conditions.OldDll4Myc.rds")

table(Liver@meta.data$Condition)

Dll4 <- subset(Liver, subset = Condition == "Dll4KO")

dll4.genes <- wilcoxauc(Dll4, 'FigClustering')

# we have all the genes for each cluster
dplyr::count(dll4.genes, group)

msigdbr_show_species()

m_df<- msigdbr(species = "Mus musculus", category = 13)

h_df<- msigdbr(species = "Homo sapiens", category = 7)


table(m_df$gs_subcat)

table(h_gene_sets$gs_cat)

# Using Hallmark genes only

h_gene_sets = msigdbr(species = "mouse", category = "H")

head(m_df)

fgsea_sets<- h_gene_sets %>% split(x = .$gene_symbol, f = .$gs_name)

hypoxia_genes <- fgsea_sets$HALLMARK_HYPOXIA

table(Dll4@meta.data$FigClustering)

# Myc Targets


dll4.genes %>%
  dplyr::filter(group == "C4 - Endothelial tip cells") %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 20)

# Good to see ODC1, KCNE3, PEG10, TPI1, MIF

# select only the feature and auc columns for fgsea, which statistics to use is an open question

cluster4.genes<- dll4.genes %>%
  dplyr::filter(group == "C4 - Endothelial tip cells") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)

# padj

cluster4.genes<- dll4.genes %>%
  dplyr::filter(group == "C4 - Endothelial tip cells") %>%
  arrange(padj) %>% 
  dplyr::select(feature, padj)

# logFC. Does not give good results

cluster4.genes<- dll4.genes %>%
  dplyr::filter(group == "C4 - Endothelial tip cells") %>%
  arrange(desc(logFC)) %>% 
  dplyr::select(feature, logFC)


cluster4.genes$feature <- tolower(cluster4.genes$feature) %>% str_to_title()

genesTables <- cluster4.genes %>%
  +   mutate(rank = rank(cluster4.genes$padj,  ties.method = "random")) 

ranks<- deframe(cluster4.genes)

head(ranks)

rank = rank(cluster4.genes$padj,  ties.method = "random")

ranks

# Trial

fgseaResShort<- fgsea(fgsea_sets, stats = ranks, nperm = 1000, scoreType = "pos")

# Full one

fgseaRes<- fgsea(fgsea_sets, stats = ranks)

# Tidy Data with auc

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))


fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  head()

write.xlsx(fgseaResTidy, "Tables/fgseaResTidy.auc.xlsx")

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


plotEnrichment(fgsea_sets[["HALLMARK_MYC_TARGETS_V1"]],
               ranks) + labs(title="HALLMARK_MYC_TARGETS_V1")


######################################################################################################

# Using escape to produce other interesting plots

Idents(Dll4) <- "FigClustering"

# Palette

colors <- colorRampPalette(c("#0348A6", "#7AC5FF", "#C6FDEC", "#FFB433", "#FF4B20"))
colors <- c("#0348A6", "#7AC5FF", "#C6FDEC", "#FFB433", "#FF4B20")
my_palette_Rui_colors <- c("#E95A74", "#50B6EF", "#45FF8E", "#F4A753", "#A80519", "#880088", "#E28CF4", "#C1B80C","#FC0808",  "#0E47D8")
my_palette_Rui_colors_B <- c("#F4A753", "#E95A74", "#50B6EF", "#FC0808",  "#0E47D8", "#45FF8E", "#A80519", "#E28CF4", "#880088", "#C1B80C")
BuRd <- c("#052F60", "#2066AC", "#4392C3", "#92C5DE", "#D1E4EF", "#F7F7F7", "#FCDBC6", "#F4A482", "#D55F4D", "#B1172B", "#67001F")

###############################################################################

# Script

GS.hallmark <- getGeneSets(library = "H", species = "Mus musculus")

# Use the Homo sapiens Gene Set because Carlos has stored the gene names in Homo sapiens nomenclature

GS.hallmark.HS <- getGeneSets(library = "H")

ES.seurat <- enrichIt(obj = Dll4, gene.sets = GS.hallmark.HS, groups = 923, cores = 2)

sce <- as.SingleCellExperiment(Dll4, assay = "RNA")

ES.sce <- enrichIt(obj = sce, gene.sets = GS.hallmark.HS, groups = 923, cores = 2)

## if working with a Seurat object
Dll4 <- Seurat::AddMetaData(Dll4, ES.seurat)

## if working with an SCE object
met.data <- merge(colData(sce), ES.sce, by = "row.names", all=TRUE)
row.names(met.data) <- met.data$Row.names
met.data$Row.names <- NULL
colData(sce) <- met.data

########################################################################################

# Heatmaps

dittoHeatmap(Dll4, genes = NULL, metas = names(ES.seurat), 
             annot.by = "FigClustering",
             annot.colors = my_palette_Rui_colors_B,
             fontsize = 7, 
             cluster_cols = F,
             heatmap.colors = BuRd)

dittoHeatmap(Dll4, genes = NULL, 
             metas = c("HALLMARK_MYC_TARGETS_V1", "HALLMARK_DNA_REPAIR", "HALLMARK_OXIDATIVE_PHOSPHORYLATION"), 
             annot.by = "FigClustering",
             annot.colors = my_palette_Rui_colors_B,
             fontsize = 7,
             heatmap.colors = BuRd)

######################################################################################

# Violin Plots

dittoPlot(Dll4, var = c("HALLMARK_MYC_TARGETS_V1"),
                group.by = "FigClustering", plots = c("jitter", "vlnplot", "boxplot"), legend.show = T, color.panel = my_palette_Rui_colors_B,
                ylab = "Enrichment Scores", 
                theme = theme_classic() + theme(axis.ticks.x = element_blank(), axis.text.x.bottom = element_blank(), plot.title = element_text(size = 10, hjust = 0.5)))

cairo_pdf("Plots/GSEA/GSEA_VlnPlot.pdf", width = 14, height = 8, family = "Arial")
dittoPlot(Dll4, var = c("HALLMARK_MYC_TARGETS_V1"),
          group.by = "FigClustering", plots = c("jitter", "vlnplot", "boxplot"), legend.show = T, color.panel = my_palette_Rui_colors_B,
          ylab = "Enrichment Scores", 
          theme = theme_classic() + theme(axis.ticks.x = element_blank(), axis.text.x.bottom = element_blank(), plot.title = element_text(size = 10, hjust = 0.5)))

dev.off()

jpeg("Plots/GSEA/GSEA_VlnPlot.jpeg", width = 14, height = 8, units = 'in', res = 800)
dittoPlot(Dll4, var = c("HALLMARK_MYC_TARGETS_V1"),
          group.by = "FigClustering", plots = c("jitter", "vlnplot", "boxplot"), legend.show = T, color.panel = my_palette_Rui_colors_B,
          ylab = "Enrichment Scores", 
          theme = theme_classic() + theme(axis.ticks.x = element_blank(), axis.text.x.bottom = element_blank(), plot.title = element_text(size = 10, hjust = 0.5)))

dev.off()

# Example with multi_dittoPlot

multi_dittoPlot(Dll4, vars = c("HALLMARK_MYC_TARGETS_V1"), ncol = 1,
                group.by = "FigClustering", plots = c("jitter", "vlnplot", "boxplot"), legend.show = T,
                ylab = "Enrichment Scores",
                theme = theme_classic() + theme(axis.ticks.x = element_blank(), axis.text.x.bottom = element_blank(), plot.title = element_text(size = 10, hjust = 0.5)))

###################################################################################

# HexDensity plots

dittoScatterHex(sce, x.var = "HALLMARK_MYC_TARGETS_V1", 
                y.var = "HALLMARK_KRAS_SIGNALING_DN", 
                colors = my_palette_Rui_colors_B,
                do.label = T,
                do.contour = TRUE) + 
  scale_fill_gradientn(colors = BuRd) 

###################################################################################

# RidgePlots

## Seurat object example
ES2 <- data.frame(Dll4[[]], Idents(Dll4))
colnames(ES2)[ncol(ES2)] <- "cluster"

## plot
ridgeEnrichment(ES2, gene.set = "HALLMARK_DNA_REPAIR", group = "cluster", add.rug = TRUE)


cairo_pdf("Plots/GSEA/GSEA_RidgePlot.pdf", width = 14, height = 8, family = "Arial")
ridgeEnrichment(ES2, gene.set = "HALLMARK_MYC_TARGETS_V1", group = "FigClustering", 
                facet = "Condition", add.rug = TRUE, colors = my_palette_Rui_colors_B)
dev.off()

jpeg("Plots/GSEA/GSEA_RidgePlot.jpeg", width = 14, height = 8, units = 'in', res = 800)
ridgeEnrichment(ES2, gene.set = "HALLMARK_MYC_TARGETS_V1", group = "FigClustering", 
                facet = "Condition", add.rug = TRUE, colors = my_palette_Rui_colors_B)
dev.off()

###################################################################################

# Split Violin Plots

splitEnrichment(ES2, split = "FigClustering", gene.set = "HALLMARK_DNA_REPAIR", colors = my_palette_Rui_colors_B)

###################################################################################

# Expanded Analysis

PCA <- performPCA(ES2, groups = "HALLMARK_DNA_REPAIR")

ES2[["HALLMARK_DNA_REPAIR"]]



pcaEnrichment(PCA, PCx = "PC1", PCy = "PC2", contours = TRUE)
###################################################################################

#  Get Significance

output <- getSignificance(ES2, group = "FigClustering", fit = "linear.model")

#####################################################################################

# Quick look in the avgexpDll4 data

avgexpDll4 <- AverageExpression(Dll4, return.seurat = T)

ES.seurat.Dll4 <- enrichIt(obj = avgexpDll4, gene.sets = GS.hallmark.HS, groups = 10, cores = 2)

avgexpDll4 <- Seurat::AddMetaData(avgexpDll4, ES.seurat.Dll4)

levels(avgexpDll4) <- c("C0 - Unspecified quiescent capillaries", "C1a - Arterial capillaries", "C1v - Venous capillaries", "C2a - Large arteries", "C2v - Large veins",  "C3 - Activated capillaries", 
                    "C4 - Endothelial tip cells", "C5ip - Proliferating S-phase", "C5p - G2-M", "C6 - X","C1")

avgexpDll4@meta.data$FigClustering <- avgexpDll4@active.ident

dittoHeatmap(avgexpDll4, genes = NULL, metas = names(ES.seurat), 
                             annot.by = "FigClustering",
                             annot.colors = my_palette_Rui_colors_B,
                             fontsize = 7, 
                             cluster_cols = F,
                             heatmap.colors.max.scaled = BuRd,
                             heatmap.colors = BuRd,
                             complex = T, data.out = F,
                              scaled.to.max = T,
                              main = expression(italic(Dll4)^iDEC))


################################################################################

# Trying out the Carlos Liver sample (Ctrl, Dll4KO, Notch1KO, RbpjKO)

Liver <- readRDS("S:/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Liver_MacarenaFC/Figures_Paper/rds/Liver_20Jan21.Figures.rds")

GS.hallmark.HS <- getGeneSets(library = "H")

table(Liver@meta.data$Condition)

ES.seurat <- enrichIt(obj = Liver, gene.sets = GS.hallmark.HS, groups = 3780, cores = 2)

sce <- as.SingleCellExperiment(Liver, assay = "RNA")

ES.sce <- enrichIt(obj = sce, gene.sets = GS.hallmark.HS, groups = 3780, cores = 2)

## if working with a Seurat object
Liver <- Seurat::AddMetaData(Liver, ES.seurat)

## if working with an SCE object
met.data <- merge(colData(sce), ES.sce, by = "row.names", all=TRUE)
row.names(met.data) <- met.data$Row.names
met.data$Row.names <- NULL
colData(sce) <- met.data

saveRDS(Liver, "./rds/Liver.Figures.GSEAHallmark.Ctrl.Dll4KO.Notch1KO.RbpjKO.rds")

########################################################################################

# Doing the GSEA with the Seurat tutorial, to get the proper NES values


liver.genes <- wilcoxauc(Liver, 'FigClustering')

# we have all the genes for each cluster
dplyr::count(liver.genes, group)

msigdbr_show_species()

m_df<- msigdbr(species = "Mus musculus", category = 13)

h_df<- msigdbr(species = "Homo sapiens", category = 7)

# Using Hallmark genes only

h_gene_sets = msigdbr(species = "mouse", category = "H")

head(m_df)

fgsea_sets<- h_gene_sets %>% split(x = .$gene_symbol, f = .$gs_name)

fgsea_sets$HALLMARK_MYC_TARGETS_V1

table(Liver@meta.data$FigClustering)

# Myc Targets


liver.genes %>%
  dplyr::filter(group == "C4 - Endothelial tip cells") %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 20)

cluster4.genes<- liver.genes %>%
  dplyr::filter(group == "C4 - Endothelial tip cells") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)

# Good to see ODC1, KCNE3, PEG10, TPI1, MIF

# select only the feature and auc columns for fgsea, which statistics to use is an open question

liver.genes$feature <- tolower(liver.genes$feature) %>% str_to_title()

cluster.genes <- liver.genes %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc) 

genesTables <- cluster.genes %>% mutate(rank = rank(cluster.genes$auc,  ties.method = "random")) 

ranks<- deframe(cluster.genes)

head(ranks)

# Full one

fgseaRes<- fgsea(fgsea_sets, stats = ranks)

# Tidy Data with auc

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))


fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  head()

write.xlsx(fgseaResTidy, "Tables/fgseaResTidy.auc.xlsx")

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


plotEnrichment(fgsea_sets[["HALLMARK_MYC_TARGETS_V1"]],
               ranks) + labs(title="HALLMARK_MYC_TARGETS_V1")

plotEnrichment(fgsea_sets[["HALLMARK_HYPOXIA"]],
               ranks) + labs(title="HALLMARK_HYPOXIA")



#######################################################################################

# Heatmap per cluster with AverageExpression

Liver <- readRDS("S:/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Liver_MacarenaFC/Figures_Paper/rds/Liver.Figures.GSEAHallmark.Ctrl.Dll4KO.Notch1KO.RbpjKO.rds")

gene.set <- sample(x = rownames(x = object@data), size = 100, replace = FALSE)

avgexp <- AverageExpression(Liver, return.seurat = T)


ES.seurat <- enrichIt(obj = avgexp, gene.sets = GS.hallmark.HS, groups = 10, cores = 2)

avgexp <- Seurat::AddMetaData(avgexp, ES.seurat)

levels(avgexp) <- c("C0 - Unspecified quiescent capillaries", "C1a - Arterial capillaries", "C1v - Venous capillaries", "C2a - Large arteries", "C2v - Large veins",  "C3 - Activated capillaries", 
                    "C4 - Endothelial tip cells", "C5ip - Proliferating S-phase", "C5p - G2-M", "C6 - X","C1")

avgexp@meta.data$FigClustering <- avgexp@active.ident

cairo_pdf("Plots/GSEA/Heatmaps/GSEA_Heatmap_Liver_Ctrl_Dll4_Notch1_Rbpj_unscaled.pdf", width = 14, height = 12, family = "Arial")
dittoHeatmap(avgexp, genes = NULL, metas = names(ES.seurat), 
             annot.by = "FigClustering",
             annot.colors = my_palette_Rui_colors_B,
             fontsize = 7, 
             cluster_cols = F,
             heatmap.colors = BuRd,
             # heatmap.colors.max.scaled = BuRd,
             scaled.to.max = F,
             scale = "none",
             complex = F, data.out = F,
             # legend_breaks = c(1, 0.8, 0.6, 0.4, 0.2, 0),
             main = expression(Control~italic(Dll4)^iDEC~italic(Notch1)^iDEC~italic(Rbpj)^iDEC))
dev.off()

cairo_pdf("Plots/GSEA/Heatmaps/GSEA_Heatmap_Liver_Ctrl_Dll4_Notch1_Rbpj_scaled.pdf", width = 14, height = 12, family = "Arial")
dittoHeatmap(avgexp, genes = NULL, metas = names(ES.seurat), 
             annot.by = "FigClustering",
             annot.colors = my_palette_Rui_colors_B,
             fontsize = 7, 
             cluster_cols = F,
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
             heatmap.colors = BuRd,
             complex = T, data.out = T)

list_mat_heatmap <- data_heatmap[["mat"]]

list_mat_heatmap.m = melt(list_mat_heatmap)
head(list_mat_heatmap.m)

list_mat_heatmap.m$value <- round(list_mat_heatmap.m$value, 2)

head(list_mat_heatmap.m)

ggplot(data = data.frame(list_mat_heatmap.m), aes(x = Var2, y = Var1, fill= value)) + geom_tile(color = "white")+ geom_text(aes(label = value), size = 3)+
  theme_minimal() + theme(text = element_text(size = 10), axis.title = element_blank(), axis.text.x.bottom = element_text(angle = 90, hjust = 1))+scale_fill_gradientn(colours = BuRd)



cairo_pdf("Plots/GSEA/GSEA_Heatmap_Liver_Ctrl_Dll4_Notch1_Rbpj.pdf", width = 14, height = 12, family = "Arial")
dittoHeatmap(avgexp, genes = NULL, metas = names(ES.seurat), 
             annot.by = "FigClustering",
             annot.colors = my_palette_Rui_colors_B,
             fontsize = 7, 
             cluster_cols = F,
             heatmap.colors = BuRd,
             complex = T, data.out = F)
dev.off()

jpeg("Plots/GSEA/GSEA_Heatmap_Liver_Ctrl_Dll4_Notch1_Rbpj.jpeg", width = 14, height = 12, units = 'in', res = 800)
dittoHeatmap(avgexp, genes = NULL, metas = names(ES.seurat), 
             annot.by = "FigClustering",
             annot.colors = my_palette_Rui_colors_B,
             fontsize = 7, 
             cluster_cols = F,
             heatmap.colors = BuRd)
dev.off()


##################################################

# Trying with ClusterProfiler

# SET THE DESIRED ORGANISM HERE
organism = "org.Mm.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

# reading in data from deseq2

liver.genes <- wilcoxauc(Liver, 'FigClustering')
df = liver.genes

# we want the log2 fold change 
original_gene_list <- df$logFC

# name the vector
names(original_gene_list) <- df$feature %>% tolower() %>% str_to_title()

original_gene_list

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

View(gene_list)

# Analysis with gseGO

gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "ALIAS", 
             minGSSize = 10, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none",
             by = "fgsea")


# Using MSigDb sets

msig_m <- msigdbr(species = "Mus musculus", category = "H")

msig_m <- msigdbr(species = "Mus musculus", category = "H") %>%
  dplyr::select(gs_name, gene_symbol) %>%
  dplyr::rename(ont = gs_name, gene = gene_symbol)
msig_m

liver.genes.d <- as.data.frame(liver.genes)

msig_ora <- enricher(gene = liver.genes.d$feature,
                     TERM2GENE = msig_m) %>%
  as_tibble
msig_ora

###################################GO:BP GSEA#########################################

m_df<- msigdbr(species = "Mus musculus", category = 13)

h_df<- msigdbr(species = "Homo sapiens", category = 7)


table(m_df$gs_subcat)

table(h_gene_sets$gs_cat)

# Using GO:BP subcategory genes only

h_gene_sets = msigdbr(species = "mouse", subcategory = "GO:BP")

fgsea_sets<- h_gene_sets %>% split(x = .$gene_symbol, f = .$gs_name) %>% list.filter("Type" > 50)

fgsea_sets$GOBP_AGING

#  After this use the GSEA for loop to get the fgsea table

liver.genes <- wilcoxauc(Liver, 'FigClustering')

dplyr::count(cluster.genes, group)

figclustering <- c("C0 - Unspecified quiescent capillaries", "C1a - Arterial capillaries", "C1v - Venous capillaries", "C2a - Large arteries", "C2v - Large veins",  "C3 - Activated capillaries", 
                   "C4 - Endothelial tip cells", "C5ip - Proliferating S-phase", "C5p - G2-M", "C6 - X")


liver.genes %>%
  dplyr::filter(group == figclustering[2]) %>%
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

class(myfgsea)

figclustering_short <- c("C0", "C1a", "C1v", "C2a", "C2v",  "C3", "C4", "C5ip", "C5p", "C6")

names(myfgsea) <- figclustering_short

write.xlsx(myfgsea, file = "./Tables/fgsea_Clusters_GO_BP.xlsx")

myfgsea.filter <- vector('list', 10)

for (i in 1:10) {
  myfgsea.filter[[i]] <- local({
    i <- i
    myfgsea[[figclustering_short[i]]][size > 50 & padj < 0.05]
  })
}

names(myfgsea.filter) <- figclustering_short


write.xlsx(myfgsea.filter, file = "./Tables/fgsea_Clusters_GO_BP.size.filter.xlsx")

read_excel_allsheets <- function(filename, tibble = FALSE) {
  # I prefer straight data.frames
  # but if you like tidyverse tibbles (the default with read_excel)
  # then just pass tibble = TRUE
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

fgsea.GOBP <- read_excel_allsheets("./Tables/fgsea_Clusters_GO_BP.size.filter.xlsx")


###########################################################Supplementary###################

#  Suplementary 6D GSEA:BP Dll4 v Dll4/Myc

Liver <- readRDS("//Tierra/SC/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Liver_MacarenaFC/Figures_Paper/rds/Liver_20Jan21.Figures.All.Conditions.OldDll4Myc.rds")
table(Liver@meta.data$Condition)

# Dll4 v Dll4/MycKO

Liver <- subset(Liver, subset = Condition == "Dll4KO" | Condition == "Dll4/MycKO")

# Dll4KO v Dll4KO+aVEGF

Liver <- subset(Liver, subset = Condition == "Dll4KO" | Condition == "D4KO_aVEGF")

# GSEA Analysis

liver.genes <- wilcoxauc(Liver, 'Condition')

dplyr::count(liver.genes, group)

condition <- c("Dll4KO", "D4KO_aVEGF")

liver.genes %>%
  dplyr::filter(group == condition[2]) %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 20)


myfgsea <- vector('list', 2)


for (i in 1:2) {
  
  dll4_v_dll4myc.genes<- liver.genes %>%
    dplyr::filter(group == condition[i]) %>%
    arrange(desc(auc)) %>% 
    dplyr::select(feature, auc)
  
  dll4_v_dll4myc.genes$feature <- tolower(dll4_v_dll4myc.genes$feature) %>% str_to_title()
  
  ranks<- deframe(dll4_v_dll4myc.genes)
  
  myfgsea[[i]] <- local({
    i <- i
    fgseaRes<- fgsea(fgsea_sets, stats = ranks)
  })
  
}

myfgsea[1]

class(myfgsea)

condition_short <- c("Dll4KO", "D4KO_aVEGF")

names(myfgsea) <- condition_short

write.xlsx(myfgsea, file = "./Tables/fgsea_Sup7C_Condition_Dll4_v_Dll4aVEGF_GO_BP.xlsx")

write.xlsx(myfgsea, file = "./Tables/fgsea_Sup6D_Condition_Dll4_v_Dll4Myc_GO_BP.xlsx")

myfgsea.filter <- vector('list', 2)

for (i in 1:2) {
  myfgsea.filter[[i]] <- local({
    i <- i
    myfgsea[[condition_short[i]]][size > 50 & padj < 0.05]
  })
}

names(myfgsea.filter) <- condition_short

write.xlsx(myfgsea, file = "./Tables/fgsea_Sup7C_Condition_Dll4_v_Dll4aVEGF_GO_BP.size.filter.xlsx")

write.xlsx(myfgsea.filter, file = "./Tables/fgsea_Sup6D_Condition_Dll4_v_Dll4Myc_GO_BP.size.filter.xlsx")

read_excel_allsheets <- function(filename, tibble = FALSE) {
  # I prefer straight data.frames
  # but if you like tidyverse tibbles (the default with read_excel)
  # then just pass tibble = TRUE
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

fgsea.GOBP <- read_excel_allsheets("./Tables/fgsea_Sup6D_Condition_Dll4_v_Dll4Myc_GO_BP.size.filter.xlsx")


########################################################


cond.genes <- wilcoxauc(Liver, 'Condition')

Dll4vCtl <- subset(Liver, subset = Condition == "RBPJKO" | Condition == "Control")

cond.genes <- wilcoxauc(Dll4vCtl, 'Condition')

# we have all the genes for each cluster
dplyr::count(cond.genes, group)

msigdbr_show_species()

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
  dplyr::filter(group == "RBPJKO") %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 20)

# Good to see ODC1, KCNE3, PEG10, TPI1, MIF

# select only the feature and auc columns for fgsea, which statistics to use is an open question

cond.dll4.genes<- cond.genes %>%
  dplyr::filter(group == "RBPJKO") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)


cond.dll4.genes$feature <- tolower(cond.dll4.genes$feature) %>% str_to_title()

genesTables <- cond.dll4.genes %>%
  +   mutate(rank = rank(cond.dll4.genes$padj,  ties.method = "random")) 

ranks<- deframe(cond.dll4.genes)

head(ranks)

ranks

# Trial

fgseaResShort<- fgsea(fgsea_sets, stats = ranks, nperm = 1000, scoreType = "pos")

# Full one

fgseaRes<- fgsea(fgsea_sets, stats = ranks)

# Tidy Data with auc

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))


fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  head()

write.xlsx(fgseaResTidy, "Tables/fgseaResTidy.RbpjvCtl.auc.xlsx")

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

