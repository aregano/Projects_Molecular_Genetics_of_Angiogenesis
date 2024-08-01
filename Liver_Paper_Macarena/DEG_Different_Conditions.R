# Libraries
library(Seurat)
library(openxlsx)
library(VennDiagram)
library("ggVennDiagram")
library("gridExtra")    


Liver <- readRDS("S:/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Liver_MacarenaFC/Figures_Paper/rds/Liver_20Jan21.Figures.All.Conditions.OldDll4Myc.rds")

table(Liver@meta.data$Condition)

Idents(Liver) <- "Condition"

# Ctrl vs Dll4, Ctrl vs Dll4/Myc and Ctrl vs Dll4 + Anti-VEGFA.

# Ctrl vs Dll4,

DEG_CtrlvDll4 <- FindMarkers(Liver, logfc.threshold = 0.25, ident.1 = "Dll4KO", ident.2 = "Control", test.use = "wilcox", only.pos = F, verbose = T, min.diff.pct = -Inf)
filter_DEG_CtrlvDll4 <- filter(DEG_CtrlvDll4, p_val_adj < 0.05)
#  Upregulated
count(filter_DEG_CtrlvDll4, avg_log2FC > 0.25)
#  Downregulated
count(filter_DEG_CtrlvDll4, avg_log2FC < -0.25)

# logFCthreshold > 0.5
count(filter_DEG_CtrlvDll4, avg_log2FC > 0.5 | avg_log2FC < -0.5 )
#  Upregulated
count(filter_DEG_CtrlvDll4, avg_log2FC > 0.5)
#  Downregulated
count(filter_DEG_CtrlvDll4, avg_log2FC < -0.5 )

# Ctrl vs Dll4/Myc

DEG_CtrlvDll4Myc <- FindMarkers(Liver, logfc.threshold = 0.25, ident.1 = "Dll4/MycKO", ident.2 = "Control", test.use = "wilcox", only.pos = F, verbose = T, min.diff.pct = -Inf)
filter_DEG_CtrlvDll4Myc <- filter(DEG_CtrlvDll4Myc, p_val_adj < 0.05)
#  Upregulated
count(filter_DEG_CtrlvDll4Myc, avg_log2FC > 0.25)
#  Downregulated
count(filter_DEG_CtrlvDll4Myc, avg_log2FC < -0.25)

# logFCthreshold > 0.5
count(filter_DEG_CtrlvDll4Myc, avg_log2FC > 0.5 | avg_log2FC < -0.5 )
#  Upregulated
count(filter_DEG_CtrlvDll4Myc, avg_log2FC > 0.5)
#  Downregulated
count(filter_DEG_CtrlvDll4Myc, avg_log2FC < -0.5 )


# Ctrl vs Dll4 + Anti-VEGFA

DEG_CtrlvDll4aVEGF <- FindMarkers(Liver, logfc.threshold = 0.25, ident.1 = "D4KO_aVEGF", ident.2 = "Control", test.use = "wilcox", only.pos = F, verbose = T, min.diff.pct = -Inf)
filter_DEG_CtrlvDll4aVEGF <- filter(DEG_CtrlvDll4aVEGF, p_val_adj < 0.05)
#  Upregulated
count(filter_DEG_CtrlvDll4aVEGF, avg_log2FC > 0.25)
#  Downregulated
count(filter_DEG_CtrlvDll4aVEGF, avg_log2FC < -0.25)

# logFCthreshold > 0.5
count(filter_DEG_CtrlvDll4aVEGF, avg_log2FC > 0.5 | avg_log2FC < -0.5 )
#  Upregulated
count(filter_DEG_CtrlvDll4aVEGF, avg_log2FC > 0.5)
#  Downregulated
count(filter_DEG_CtrlvDll4aVEGF, avg_log2FC < -0.5 )


# Dll4KO v Dll4/MycKO

DEG_Dll4vDll4Myc <- FindMarkers(Liver, logfc.threshold = 0.25, ident.1 = "Dll4/MycKO", ident.2 = "Dll4KO", test.use = "wilcox", only.pos = F, verbose = T, min.diff.pct = -Inf)
filter_DEG_Dll4vDll4Myc <- filter(DEG_Dll4vDll4Myc, p_val_adj < 0.05)
#  Upregulated
count(filter_DEG_Dll4vDll4Myc, avg_log2FC > 0.25)
#  Downregulated
count(filter_DEG_Dll4vDll4Myc, avg_log2FC < -0.25)

# logFCthreshold > 0.5
count(filter_DEG_Dll4vDll4Myc, avg_log2FC > 0.5 | avg_log2FC < -0.5 )
#  Upregulated
count(filter_DEG_Dll4vDll4Myc, avg_log2FC > 0.5)
#  Downregulated
count(filter_DEG_Dll4vDll4Myc, avg_log2FC < -0.5 )



# Dll4KO v Dll4KOaVEGF

DEG_Dll4vDll4aVEGF <- FindMarkers(Liver, logfc.threshold = 0.25, ident.1 = "D4KO_aVEGF", ident.2 = "Dll4KO", test.use = "wilcox", only.pos = F, verbose = T, min.diff.pct = -Inf)
filter_DEG_Dll4vDll4aVEGF <- filter(DEG_Dll4vDll4aVEGF, p_val_adj < 0.05)
#  Upregulated
count(filter_DEG_Dll4vDll4aVEGF, avg_log2FC > 0.25)
#  Downregulated
count(filter_DEG_Dll4vDll4aVEGF, avg_log2FC < -0.25)

# logFCthreshold > 0.5
count(filter_DEG_Dll4vDll4aVEGF, avg_log2FC > 0.5 | avg_log2FC < -0.5 )
#  Upregulated
count(filter_DEG_Dll4vDll4aVEGF, avg_log2FC > 0.5)
#  Downregulated
count(filter_DEG_Dll4vDll4aVEGF, avg_log2FC < -0.5 )

# Merge

DEG_total <- list(filter_DEG_CtrlvDll4, filter_DEG_CtrlvDll4Myc, filter_DEG_CtrlvDll4aVEGF, filter_DEG_Dll4vDll4Myc, filter_DEG_Dll4vDll4aVEGF)

names(DEG_total) <- c("CtrlvDll4", "CtrlvDll4Myc", "CtrlvDll4aVEGF", "Dll4vDll4Myc", "Dll4vDll4aVEGF")

write.xlsx(DEG_total, "Tables/DEG/DEG_Numbers_Conditions.wilcox.logfc0.25.padj0.05.xlsx", rowNames = T)


# Venn Diagram of DEGs with the different conditions

DEGs_Upregulated <- read_excel("Tables/DEG/DEG_Numbers_Conditions.wilcox.logfc0.25.padj0.05.xlsx", sheet = 8)

DEGs_Downregulated <- read_excel("Tables/DEG/DEG_Numbers_Conditions.wilcox.logfc0.25.padj0.05.xlsx", sheet = 9)

display_venn(
  DEGs_Upregulated,
  category.names = c("Set 1" , "Set 2 " , "Set 3", "Set 4", "Set5"),
  fill = c("#999999", "#E69F00", "#56B4E9", "#009E73", "black")
)

ggVennDiagram(DEGs_Upregulated, label_alpha = 0)

venn.diagram(DEGs_Upregulated, filename = NULL)

# Venn Diagram function

display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}

Venn_DEGs_Upregulated <- list(
  DEGs_Upregulated$CtrlvDll4 %>% na.omit(),
  DEGs_Upregulated$CtrlvDll4Myc %>% na.omit(),
  DEGs_Upregulated$CtrlvDll4aVEGF %>% na.omit(),
  DEGs_Upregulated$Dll4vDll4Myc %>% na.omit(),
  DEGs_Upregulated$Dll4vDll4aVEGF %>% na.omit()
  )

Venn_DEGs_Upregulated_3groups <- list(
  DEGs_Upregulated$CtrlvDll4 %>% na.omit(),
  DEGs_Upregulated$CtrlvDll4Myc %>% na.omit(),
  DEGs_Upregulated$CtrlvDll4aVEGF %>% na.omit()
)

Venn_DEGs_Downregulated <- list(
  DEGs_Downregulated$CtrlvDll4 %>% na.omit(),
  DEGs_Downregulated$CtrlvDll4Myc %>% na.omit(),
  DEGs_Downregulated$CtrlvDll4aVEGF %>% na.omit(),
  DEGs_Downregulated$Dll4vDll4Myc %>% na.omit(),
  DEGs_Downregulated$Dll4vDll4aVEGF %>% na.omit()
)

Venn_DEGs_Downregulated_3groups <- list(
  DEGs_Downregulated$CtrlvDll4 %>% na.omit(),
  DEGs_Downregulated$CtrlvDll4Myc %>% na.omit(),
  DEGs_Downregulated$CtrlvDll4aVEGF %>% na.omit()
)
cairo_pdf("Plots/Figure_7/Supplementary/Venn_Diagram_DEGs_Upregulated_5groups.pdf", width = 5, height = 5, family = "Arial")
display_venn(Venn_DEGs_Upregulated,   category.names = c("CtrlvDll4" , "CtrlvDll4Myc" , "CtrlvDll4aVEGF", "Dll4vDll4Myc", "Dll4vDll4aVEGF"),
fill = c("#F94040", "#C19EDD", "#BB005E", "#009E73", "#56B4E9"), main = "DEGs Upregulated 5groups")
dev.off()

cairo_pdf("Plots/Figure_7/Supplementary/Venn_Diagram_DEGs_Upregulated_3groups.pdf", width = 5, height = 5, family = "Arial")
display_venn(Venn_DEGs_Upregulated_3groups,   category.names = c("CtrlvDll4" , "CtrlvDll4Myc" , "CtrlvDll4aVEGF"),
             fill = c("#F94040", "#C19EDD", "#BB005E"), main = "DEGs Upregulated 3groups")
dev.off()

cairo_pdf("Plots/Figure_7/Supplementary/Venn_Diagram_DEGs_Downregulated_3groups.pdf", width = 5, height = 5, family = "Arial")
display_venn(Venn_DEGs_Downregulated_3groups,   category.names = c("CtrlvDll4" , "CtrlvDll4Myc" , "CtrlvDll4aVEGF"),
             fill = c("#F94040", "#C19EDD", "#BB005E"), main = "DEGs Downregulated 3groups")
dev.off()