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
library(magrittr)


ranksd4 <- read.xlsx("./Tables/GSEA_Fig6/ranks.Dll4vsCtrl.xlsx", check.names = T, rowNames = F, rows = 1:32263)

ranksd4myc <- read.xlsx("./Tables/GSEA_Fig6/ranks.Dll4MycvsCtrl.xlsx", check.names = T, rowNames = F, rows = 1:32263)

ranksd4avegdf <- read.xlsx("./Tables/GSEA_Fig7/ranksDll4KOaVEGF.auc.xlsx", check.names = T, rowNames = F, rows = 1:32263)

ranks <- list(ranksd4, ranksd4myc, ranksd4avegdf)

# ranks <- read_excel_allsheets("Tables/GSEA_Fig7/ranks.Dll4KO.Dll4MycKO.Dll4KOaVEGF.auc.xlsx", tibble = T)
  
ranks

h_gene_sets = msigdbr(species = "mouse", category = "H")

head(m_df)

fgsea_sets<- h_gene_sets %>% split(x = .$gene_symbol, f = .$gs_name)

class(fgsea_sets)

E2F_genes <- fgsea_sets$HALLMARK_E2F_TARGETS
Myc1_genes <- fgsea_sets$HALLMARK_MYC_TARGETS_V1
OxPhos_genes <- fgsea_sets$HALLMARK_OXIDATIVE_PHOSPHORYLATION

GeneSet_genes <- list(E2F = E2F_genes, Myc = Myc1_genes, OxPhos = OxPhos_genes)

GeneSet_genes$E2F


# E2F_rank <- for (i in length(ranks)) {ranks[ranks[[i]][["gene"]] %in% E2F_genes] %>% extract2(1) %>% as.character()
# }

# E2F ranks

E2F_ranksd4 <- ranksd4[ranksd4$gene %in% GeneSet_genes$E2F, ] %>% extract2(1) %>% as.character()

E2F_ranksd4myc <- ranksd4myc[ranksd4myc$X1 %in% GeneSet_genes$E2F, ] %>% extract2(1) %>% as.character()

E2F_ranksd4avegf <- ranksd4avegdf[ranksd4avegdf$X %in% GeneSet_genes$E2F, ] %>% extract2(1) %>% as.character()

E2F_ranks <- list(E2F_ranksd4 = E2F_ranksd4, E2F_ranksd4myc = E2F_ranksd4myc, E2F_ranksd4avegf = E2F_ranksd4avegf)

write.xlsx(E2F_ranks, "Tables/GSEA_Fig7/E2F_ranks.Dll4KO.Dll4MycKO.Dll4KOaVEGF.xlsx", rowNames = T)

# Myc ranks

Myc_ranksd4 <- ranksd4[ranksd4$gene %in% GeneSet_genes$Myc, ] %>% extract2(1) %>% as.character()

Myc_ranksd4myc <- ranksd4myc[ranksd4myc$X1 %in% GeneSet_genes$Myc, ] %>% extract2(1) %>% as.character()

Myc_ranksd4avegf <- ranksd4avegdf[ranksd4avegdf$X %in% GeneSet_genes$Myc, ] %>% extract2(1) %>% as.character()

Myc_ranks <- list(Myc_ranksd4 = Myc_ranksd4, Myc_ranksd4myc = Myc_ranksd4myc, Myc_ranksd4avegf = Myc_ranksd4avegf)

write.xlsx(Myc_ranks, "Tables/GSEA_Fig7/Myc_ranks.Dll4KO.Dll4MycKO.Dll4KOaVEGF.xlsx", rowNames = T)

# OxPhos ranks

OxPhos_ranksd4 <- ranksd4[ranksd4$gene %in% GeneSet_genes$OxPhos, ] %>% extract2(1) %>% as.character()

OxPhos_ranksd4myc <- ranksd4myc[ranksd4myc$X1 %in% GeneSet_genes$OxPhos, ] %>% extract2(1) %>% as.character()

OxPhos_ranksd4avegf <- ranksd4avegdf[ranksd4avegdf$X %in% GeneSet_genes$OxPhos, ] %>% extract2(1) %>% as.character()

OxPhos_ranks <- list(OxPhos_ranksd4 = OxPhos_ranksd4, OxPhos_ranksd4myc = OxPhos_ranksd4myc, OxPhos_ranksd4avegf = OxPhos_ranksd4avegf)

write.xlsx(OxPhos_ranks, "Tables/GSEA_Fig7/OxPhos_ranks.Dll4KO.Dll4MycKO.Dll4KOaVEGF.xlsx", rowNames = T)

# Get list of all top3 genes of each Condition, 27 total

E2F_genes_Condition <- read_excel_allsheets("Tables/GSEA_Fig7/E2F_ranks.Dll4KO.Dll4MycKO.Dll4KOaVEGF.xlsx")
Myc_genes_Condition <- read_excel_allsheets("Tables/GSEA_Fig7/Myc_ranks.Dll4KO.Dll4MycKO.Dll4KOaVEGF.xlsx")
OxPhos_genes_Condition <- read_excel_allsheets("Tables/GSEA_Fig7/OxPhos_ranks.Dll4KO.Dll4MycKO.Dll4KOaVEGF.xlsx")

Top3genes_GeneSets_Conditions <- c(E2F_genes_Condition[[1]][[1]][1:3], E2F_genes_Condition[[2]][[1]][1:3], E2F_genes_Condition[[3]][[1]][1:3],
                                   Myc_genes_Condition[[1]][[1]][1:3], Myc_genes_Condition[[2]][[1]][1:3], Myc_genes_Condition[[3]][[1]][1:3],
                                   OxPhos_genes_Condition[[1]][[1]][1:3], OxPhos_genes_Condition[[2]][[1]][1:3], OxPhos_genes_Condition[[3]][[1]][1:3])

Top3genes_GeneSets_Conditions

# C4 Tip Cells: Vegfa, Kcne3, Esm1, Apln
# Cp G2-M: Top2a, Mki67, Hmgb2,Cenpf
# C1a Arterial capillaries: Ntn4, Ltpb4,Msr1,Ehd3

