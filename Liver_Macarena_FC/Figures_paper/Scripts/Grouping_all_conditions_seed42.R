# Here I will gather all rds, with all conditions, and make 3 final objects:
# 1. All Conditions
# 2. Liver 2weeks deletion
# 3. Liver 4 days deletion
# In order to group them correctly, I will also downsample to a max # of cells of 1000, using a set seed, 42.

library(Seurat)

Liver_20Jan21 <- readRDS("S:/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Liver_MacarenaFC/Figures_Paper/rds/Groups/Liver_20Jan21.Figures.rds")
Liver4 <- readRDS("S:/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Liver_MacarenaFC/Figures_Paper/rds/Groups/Liver_Ctrl4d_Dll4+inh.V3.rds")
Liver_Heart <- readRDS("S:/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Liver_MacarenaFC/Figures_Paper/rds/Groups/Dll4MycLOF.Fresh.rds")
Liver_MultipleGroups <- readRDS("S:/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Liver_MacarenaFC/Figures_Paper/rds/Groups/Liver_MultipleGroups.Mapping.rds")
Liver3 <- readRDS("S:/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Liver_MacarenaFC/Figures_Paper/rds/Groups/Liver3.mkI.rds")


table(Liver_20Jan21@meta.data$Condition)
table(Liver4@meta.data$Condition)
table(Liver_Heart@meta.data$Condition)
table(Liver_MultipleGroups@meta.data$Condition)
table(Liver3@meta.data$Condition)

#  Liver_Heart

Idents(Liver_Heart) <- "Condition"

DimPlot(Liver_Heart)

Liver_Heart <- subset(Liver_Heart, subset = Condition == "Dll4/MycLOF", downsample = 1000, seed = 42)

condition <- c("Dll4/MycLOF")

names(condition) <- levels(Liver_Heart)

Liver_Heart@active.ident -> Liver_Heart@meta.data$Condition

# Liver_Multiple_Groups

Idents(Liver_MultipleGroups) <- "Condition"

Liver_MultipleGroups <- subset(Liver_MultipleGroups, downsample = 1000, seed = 42)

# Liver3

Idents(Liver3) <- "hash.ID"

Liver3 <- subset(Liver3, subset = hash.ID == "BC304", downsample = 1000, seed = 42)

condition <- c("Jag1/Jag2/Dll1KO")

names(condition) <- levels(Liver3)

Liver3 <- RenameIdents(Liver3, condition)

Liver3@active.ident -> Liver3@meta.data$Condition

# Liver4

Liver4 <- subset(Liver4, downsample = 1000, seed = 42)

# Liver_MultiplGroups

Idents(Liver_MultipleGroups) <- Liver_MultipleGroups@meta.data$Condition

Liver_MultipleGroups <- subset(Liver_MultipleGroups, downsample = 1000, seed = 42)

# Merging datasets
 
Liver1 <- merge(Liver_MultipleGroups, Liver4)
Liver2 <- merge(Liver1, Liver3)
Liver <- merge(Liver2, Liver_Heart)

table(Liver@meta.data$Condition)

saveRDS(Liver, "rds/Groups/Liver.Figures.merge.All.Conditions.rds")
