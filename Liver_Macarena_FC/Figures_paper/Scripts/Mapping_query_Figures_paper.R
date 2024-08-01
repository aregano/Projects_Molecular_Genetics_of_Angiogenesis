# Now onto making the mapping

library(Seurat)
library(ggplot2)

##############################################################################################################

# rds

Carlos <- readRDS("S:/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Liver_MacarenaFC/Figures_Paper/rds/Groups/Liver_20Jan21.Figures.rds")

Liver.query <- readRDS("./rds/Groups/Liver.Figures.merge.All.Conditions.rds")

Liver.query <- subset(Liver.query, downsample = 1000, seed = 42)

table(Liver.query@meta.data$Condition)

#  Mapping

Liver.anchors <- FindTransferAnchors(reference = Carlos, query = Liver.query,
                                     dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = Liver.anchors, refdata = Carlos$FigClustering,
                            dims = 1:30)
Liver.query <- AddMetaData(Liver.query, metadata = predictions)

# Doesn't work here
# Liver.query$prediction.match <- Liver.query$predicted.id == Liver.query$
# table(Liver.query$prediction.match)

table(Liver.query$predicted.id)

# Unimodal Map Projection

Carlos <- RunUMAP(Carlos, dims = 1:7, reduction = "pca", return.model = TRUE, seed.use = 123456)
Liver.query <- MapQuery(anchorset = Liver.anchors, reference = Carlos, query = Liver.query,
                        refdata = list(celltype = "FigClustering"), reference.reduction = "pca", reduction.model = "umap")

p1 <- DimPlot(Carlos, reduction = "umap", group.by = "FigClustering", label = TRUE, label.size = 3,
              repel = TRUE) + NoLegend() + ggtitle("Reference annotations")

p2 <- DimPlot(Liver.query, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
              label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
p1+p2

Idents(Liver.query) <- "Condition"

Liver.query@meta.data$Condition <- Liver.query@active.ident

Liver.query@meta.data$FigClustering <- Liver.query@meta.data$predicted.id

DimPlot(Liver.query, reduction = "ref.umap",split.by = "Condition", group.by = "FigClustering")

saveRDS(Liver.query, "./rds/Liver_20Jan21.Figures.All.Conditions.rds")

#  Trials

FeaturePlot(Liver.query, features = c("RSPO3", "GJA5", "KCNE3", "ODC1", "MKI67", "PEG10"))

table(Liver.query@meta.data$Condition)

##############################################################################################################

# Mapping of Liver 4 days deletion: Control + Dll4KO Veh + Dll4KO SL327

Liver4d <- subset(Liver4d, downsample = 1000, seed = 42)

Carlos <- readRDS("S:/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Liver_MacarenaFC/Figures_Paper/rds/Liver_20Jan21.Figures.rds")


table(Liver4d@meta.data$Condition)

Liver4d -> Liver.query

#  Mapping

Liver.anchors <- FindTransferAnchors(reference = Carlos, query = Liver.query,
                                     dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = Liver.anchors, refdata = Carlos$FigClustering,
                            dims = 1:30)
Liver.query <- AddMetaData(Liver.query, metadata = predictions)

# Doesn't work here
# Liver.query$prediction.match <- Liver.query$predicted.id == Liver.query$
# table(Liver.query$prediction.match)

table(Liver.query$predicted.id)

# Unimodal Map Projection

Carlos <- RunUMAP(Carlos, dims = 1:7, reduction = "pca", return.model = TRUE, seed.use = 123456)
Liver.query <- MapQuery(anchorset = Liver.anchors, reference = Carlos, query = Liver.query,
                        refdata = list(celltype = "FigClustering"), reference.reduction = "pca", reduction.model = "umap")

p1 <- DimPlot(Carlos, reduction = "umap", group.by = "FigClustering", label = TRUE, label.size = 3,
              repel = TRUE) + NoLegend() + ggtitle("Reference annotations")

p2 <- DimPlot(Liver.query, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
              label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
p1+p2

Liver.query@meta.data$Condition <- Liver.query@active.ident

Liver.query@meta.data$FigClustering <- Liver.query@meta.data$predicted.id

DimPlot(Liver.query, reduction = "ref.umap",split.by = "Condition", group.by = "FigClustering")

saveRDS(Liver.query, "./rds/Liver_20Jan21.Figures.All.Conditions.rds")

#  Trials

FeaturePlot(Liver.query, features = c("RSPO3", "GJA5", "KCNE3", "ODC1", "MKI67", "PEG10"))

table(Liver.query@meta.data$Condition)

saveRDS(Liver.query, "./rds/Liver4d.Mapping.rds")
