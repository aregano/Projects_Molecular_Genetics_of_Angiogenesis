# Starting analysis on the Liver sample

library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(dittoSeq)
library(patchwork)

###############################

# Quality Control Analysis performed by Carlos Torroja

# This is a QC check only

Liver <- readRDS("")

id <- c("Liver")
names(id) <- levels(Liver)
Liver <- RenameIdents(Liver, id)
Liver@meta.data$Sample <- Liver@active.ident
Liver@meta.data$orig.ident <- Liver@active.ident

VlnPlot(Liver, features = c("nFeature_RNA", "nCount_RNA", "subsets_MT_percent"), group.by = "orig.ident", ncol = 3, pt.size = 0)


plot1 <- FeatureScatter(Liver, feature1 = "nFeature_RNA",
                        group.by = "orig.ident", feature2 = "subsets_MT_percent") +
  geom_vline(xintercept = c(200,5000),linetype = 2 ) +
  geom_hline(yintercept = 15 ,linetype = 2)
plot2 <- FeatureScatter(Liver, feature1 = "nCount_RNA",
                        group.by = "orig.ident",feature2 = "nFeature_RNA")+
  geom_hline(yintercept = c(200,5000),linetype = 2 )
plot1 / plot2

# Carlos did not follow the exact same criteria as Irepan, the cutouts are more generous in Carlos' analysis.
# I will keep them the same since it is how he did the first Liver analysis

####################################

# EXAMPLE OF QC: used in all samples besides Ctl, Dll4KO, RbpjKO, Notch1KO

# Normalize data

Liver <- NormalizeData(Liver, 
                       normalization.method = "LogNormalize",
                       scale.factor = 10000)

# Identification of highly variable features (feature selection)

Liver <- FindVariableFeatures(Liver, 
                              selection.method = "vst",
                              nfeatures = 2000)
top10 <- head(VariableFeatures(Liver), 10)
top10

plot1 <- VariableFeaturePlot(Liver)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

# Scaling the data

all.genes <- rownames(Liver)
Liver <- ScaleData(Liver, features = all.genes)


## Perform linear dimensional reduction

Liver <- RunPCA(Liver,
                features = VariableFeatures(object = Liver))


# Examine and visualize PCA results a few different ways

print(Liver[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Liver, dims = 1:2, reduction = "pca")
DimPlot(Liver,group.by = "orig.ident", reduction = "pca")
DimHeatmap(Liver, dims = 1:9, cells = 500, balanced = TRUE)

ElbowPlot(Liver, ndims = 30)


Liver <- FindNeighbors(Liver, dims = 1:30)
Liver <- FindClusters(Liver, resolution = 0.1)
Liver <- FindClusters(Liver, resolution = 0.2)
Liver <- FindClusters(Liver, resolution = 0.35)
Liver <- FindClusters(Liver, resolution = 0.5)
Liver <- FindClusters(Liver, resolution = 0.5)


DimPlot(Liver)

#Run UMAP

Liver <- RunUMAP(Liver, dims = 1:30)

DimPlot(Liver, group.by = "RNA_snn_res.0.35")
DimPlot(Liver, group.by = "RNA_snn_res.0.35", split.by = "Condition")

table(Liver@meta.data$hash.ID)

Liver_singlets <- subset(x = Liver, subset = (hash.ID == "HTO302" | hash.ID == "HTO307" | hash.ID == "HTO306" | hash.ID == "HTO301" | hash.ID == "HTO305"))


#  Now to take out Blood cells

FeaturePlot(Liver_singlets,
            features = c("PTPRC", # Leukocyte lineage
                         "CD14", # Macrophage
                         "CD3E", # T cell                        
                         "CD79A", # B cell
                         "PECAM1", # ECs
                         "CDH5") # ECs
)

Liver_singlets <- FindNeighbors(Liver_singlets, dims = 1:30)
Liver_singlets <- FindClusters(Liver_singlets, resolution = 1.5)

DimPlot(Liver_singlets, label = T, label.box = T, group.by = "RNA_snn_res.1.5")

# As suspected, blood cells belong to the outlier clusters (22, 27, 21, 23, 20, 18). The main one in the middle is the liver cluster 


# So I will just take out the blood cells clusters (22, 27, 21, 23, 20, 18)

VlnPlot(Liver_singlets, features =  "PTPRC",
        group.by =  "RNA_snn_res.1.5", 
        pt.size = 0.05 ) + theme(legend.position="none")

Liver <- Liver_singlets@meta.data$RNA_snn_res.1.5  %in%  as.character(c(0:17, 19, 24:26, 28))

p1 <- DimPlot(Liver_singlets, reduction = "umap", cells = Liver,
              group.by = "RNA_snn_res.1.5",
              label = T, pt.size = 0.5)
p1

x <- Liver_singlets@meta.data$RNA_snn_res.1.5  %in%  as.character(c(0:17, 19, 24:26, 28))
Liver_nonblood <- subset(Liver_singlets, cells = colnames(Liver_singlets)[x])
DefaultAssay(Liver_nonblood) <- "RNA"

Liver_nonblood <- FindVariableFeatures(Liver_nonblood, 
                                       selection.method = "vst",
                                       nfeatures = 2000)
Liver_nonblood <- RunPCA(Liver_nonblood,
                         features = VariableFeatures(object = Liver_nonblood))
Liver_nonblood <- FindNeighbors(Liver_nonblood, dims = 1:20)
Liver_nonblood <- FindClusters(Liver_nonblood, resolution = 0.1)
Liver_nonblood <- FindClusters(Liver_nonblood, resolution = 0.3)
Liver_nonblood <- FindClusters(Liver_nonblood, resolution = 0.05)

Liver_nonblood <- RunUMAP(Liver_nonblood, dims = 1:20)

DimPlot(Liver_nonblood)

####################################################

# DATASET MERGING

# Here I will gather all rds, with all conditions, and make 3 final objects:
# 1. All Conditions
# 2. Liver 2weeks deletion
# 3. Liver 4 days deletion
# In order to group them correctly, I will also downsample to a max # of cells of 1000, using a set seed, 42.

library(Seurat)

Liver_20Jan21 <- readRDS("S:/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Liver_MacarenaFC/Figures_Paper/rds/Groups/Liver_20Jan21.Figures.rds")
Liver4 <- readRDS("S:/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Liver_MacarenaFC/Figures_Paper/rds/Groups/Liver_Ctrl4d_Dll4+inh.V3.rds")
Liver_MultipleGroups <- readRDS("S:/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Liver_MacarenaFC/Figures_Paper/rds/Groups/Liver_MultipleGroups.Mapping.rds")
Liver3 <- readRDS("S:/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Liver_MacarenaFC/Figures_Paper/rds/Groups/Liver3.mkI.rds")


table(Liver_20Jan21@meta.data$Condition)
table(Liver4@meta.data$Condition)
table(Liver_MultipleGroups@meta.data$Condition)
table(Liver3@meta.data$Condition)

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

Idents(Liver4) <- "Condition"

Liver4 <- subset(Liver4, downsample = 1000, seed = 42)

# Liver_MultiplGroups

Idents(Liver_MultipleGroups) <- Liver_MultipleGroups@meta.data$Condition

Liver_MultipleGroups <- subset(Liver_MultipleGroups, downsample = 1000, seed = 42)

# Merging datasets

Liver1 <- merge(Liver_MultipleGroups, Liver4)
Liver <- merge(Liver1, Liver3)


table(Liver@meta.data$Condition)

saveRDS(Liver, "rds/Groups/Liver.Figures.merge.All.Conditions.rds")

# Mapping them using the Starting scRNASeq dataset from Carlos Torroja

# Read rds object with Ctl, Dll4KO, Notch1KO and RbpjKO where initial DimRed analysis was performed

Carlos <- readRDS("")

# Read rds object where all conditions have been merged to a normalized count of 1000 cells/condition

Liver <- readRDS("")

Liver.query <- Liver
Liver.anchors <- FindTransferAnchors(reference = Carlos, query = Liver.query,
                                     dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = Liver.anchors, refdata = Carlos$NewClustering,
                            dims = 1:30)
Liver.query <- AddMetaData(Liver.query, metadata = predictions)

table(Liver.query@mata.data$predicted.id)

Liver.query@mata.data$predicted.id -> Liver.query@meta.data$FigClustering
