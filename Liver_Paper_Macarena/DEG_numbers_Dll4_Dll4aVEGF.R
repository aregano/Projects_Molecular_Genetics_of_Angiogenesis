library(Seurat)
library(ggplot2)
library(patchwork)
library("dittoSeq")
library("magick")
library("RColorBrewer")
library(cowplot)
library(readxl)
library(pheatmap)

BuRd <- c("#052F60", "#2066AC", "#4392C3", "#92C5DE", "#D1E4EF", "#F7F7F7", "#FCDBC6", "#F4A482", "#D55F4D", "#B1172B", "#67001F")

palette_Maca_HM <- c("#606060", "#F94040", "#BB005E")


Liver <- readRDS("S:/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Liver_MacarenaFC/Figures_Paper/rds/Liver_20Jan21.Figures.All.Conditions.OldDll4Myc.rds")

table(Liver@meta.data$Condition)

LiverDEG <- subset(Liver, subset = Condition == "Control" | Condition == "Dll4KO")

# Dll4KO+aVEGF v Ctrl

LiverDEG <- subset(Liver, subset = Condition == "Control" | Condition == "D4KO_aVEGF")

Idents(LiverDEG) <- "Condition"

table(LiverDEG@meta.data$Condition)

DEG <- FindMarkers(LiverDEG, ident.1 = "Dll4KO", test.use = "wilcox")

DEG <- FindMarkers(LiverDEG, ident.1 = "D4KO_aVEGF", test.use = "wilcox")

filter_DEG <- filter(DEG, p_val_adj < 0.05)

write.xlsx(filter_DEG, "Tables/DEG_Dll4vCtrl_padj0.05_logfc0.25.xlsx", colNames = T)

write.xlsx(filter_DEG, "Tables/DEG_Dll4_aVEGFvCtrl_padj0.05_logfc0.25.xlsx", colNames = T)


# Super Heatmap

genes <- list(c("STMN1", "MCM2", "CDKN1A", "MYC", "ODC1", "ESM1", "KCNE3", "ANGPT2", "APLN", "RPL32", "MSR1", "LTBP4", "WNT2"))
genes.names <- c("Stmn1", "Mcm2", "Cdkn1a", "Myc", "Odc1", "Esm1", "Kcne3", "Angpt2", "Apln", "Rpl32", "Msr1", "Ltbp4", "Wnt2")

LiverHM <- subset(Liver, subset = Condition == "Control" | Condition == "Dll4KO" | Condition == "D4KO_aVEGF")

avgexp <- AverageExpression(LiverHM, return.seurat = T)

levels(avgexp) <- c("Control", "Dll4KO", "D4KO_aVEGF")

avgexp@active.ident -> avgexp@meta.data$Condition


cairo_pdf("Plots/Figure_7/Figure_7D_Heatmap_Liver_Ctrl_Dll4_Dll4aVEGF_AllGenes_scaled.pdf", width = 4, height = 7, family = "Arial")
dittoHeatmap(avgexp, var = "Condition",
             annot.by = "Condition",
             slot = "data",
             annot.colors = palette_Maca_HM,
             fontsize = 7, 
             cluster_cols = F,
             highlight.features = c("STMN1", "MCM2", "CDKN1A", "MYC", "ODC1", "ESM1", "KCNE3", "ANGPT2", "APLN", "RPL32", "MSR1", "LTBP4", "WNT2"),
             heatmap.colors.max.scaled = BuRd,
             scale = "row",
             heatmap.colors = BuRd,
             # complex = T, 
             scaled.to.max = F,
             data.out = F,
             main = expression(Control~italic(Dll4)^iDEC~italic(Dll4)^iDEC+aVEGF))
dev.off()




# Trying with pheatmap

add.flag(p1, genes, 0.2)             



vector <- avgexp@assays[["RNA"]]@data

is.na(vector) %>% table()

is.nan(vector) %>% table()

is.infinite(vector) %>% table()

dim(vector)

head(vector)

pheatmap(vector, color = BuRd, scale = "row")

vector.wo.zero <-filter_if(vector, is.numeric, all_vars((.) != 0)) #remove all rows with n= 0
rnames <- vector.wo.zero$`Gene symbol`#select name
vector.wo.zero <- vector.wo.zero[-c(1:2)]# remove gene symbol
mouse.matrix <-(as.matrix(vector.wo.zero))
rownames(mouse.matrix) <- rnames # assign row names
mouse.matrix <- t(mouse.matrix) #transpose
mouseUT <- scale(mouse.matrix)