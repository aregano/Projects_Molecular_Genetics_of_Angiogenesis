# Now onto making the mapping

library(Seurat)
library(ggplot2)

##############################################################################################################

# rds

Carlos <- readRDS("S:/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Liver_MacarenaFC/Figures_Paper/rds/Groups/Liver_20Jan21.Figures.rds")

Liver <- readRDS("S:/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Liver_MacarenaFC/Figures_Paper/rds/Liver_20Jan21.Figures.All.Conditions.ppt_Rui.rds")

table(Liver@meta.data$Condition)

Liver.query <- readRDS("./rds/Groups/Liver.Figures.merge.All.Conditions.rds")

Liver.query <- subset(Liver.query, downsample = 1000, seed = 42)

Liver.query <- subset(Liver, subset = Condition == "Dll4KO" | Condition == "Control(4d)" | Condition == "Dll4KO(4d)+SL327" | Condition == "Dll4KO(4d)+Vehicle")


Liver.query@meta.data$Condition <- Liver.query@active.ident

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

saveRDS(Liver.query, "./rds/Liver_20Jan21.Figures.Fig7Sup.2.rds")

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


#######################################################

# Kalucka Mapping

kalucka_norm <- readRDS("S:/LAB_RB/LAB/Alvaro/Bioinformatics/Joanna_Kalucka_EC_atlas/kalucka_norm.rds")

kalucka_norm <- readRDS("/run/user/2041069230/gvfs/smb-share:server=tierra,share=sc/LAB_RB/LAB/Alvaro/Bioinformatics/Joanna_Kalucka_EC_atlas/kalucka_norm.rds")

kalucka_Liver <- kalucka_norm@meta.data$Tissue  %in%  c("liver")
kalucka_Liver<- subset(kalucka_norm, cells = colnames(kalucka_norm)[kalucka_Liver])
kalucka_Liver$Cluster <- as.factor(kalucka_Liver$Cluster)

Liver <- ScaleData(Liver)

kalucka_Liver <- subset(kalucka_norm, subset = Tissue == "liver")

# Rerunning uMAP on Kalucka_Heart

kalucka_Liver <- FindVariableFeatures(kalucka_Liver, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(kalucka_Liver)
VariableFeaturePlot(kalucka_Liver)
kalucka_Liver <- ScaleData(kalucka_Liver, features = all.genes)   
kalucka_Liver <- RunPCA(kalucka_Liver, features = VariableFeatures(object = kalucka_Liver))

kalucka_Liver <- RunUMAP(kalucka_Liver, dims = 1:30, return.model=TRUE)

DimPlot(kalucka_Liver, group.by = "Cluster")

kalucka_Liver <- RenameCells(kalucka_Liver, add.cell.id = 1)

Cells(kalucka_Liver)

Cells(Liver)

Liver <- FindVariableFeatures(Liver, selection.method = "vst", nfeatures = 2000)
VariableFeaturePlot(Liver)
all.genes <- rownames(Liver)
Liver <- ScaleData(Liver, features = all.genes)   
Liver <- RunPCA(Liver, features = VariableFeatures(object = Liver))

##################################### Doing the Mapping for Control sample v Kalucka ############################

Carlos <- readRDS("S:/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Liver_MacarenaFC/Figures_Paper/rds/Liver_20Jan21.Figures.rds")


all.genes <- rownames(Liver)

genes.human <- convertHumanGeneList(all.genes)

#  Mapping

kalucka.Anchors <- FindTransferAnchors(reference = kalucka_Liver, query = Liver,
                                       dims = 1:30, reference.reduction = "pca")

predictions <- TransferData(anchorset = kalucka.Anchors, refdata = kalucka_Liver$Cluster,
                            dims = 1:30)

Control_Heart <- AddMetaData(Control_Heart, metadata = predictions)

table(Control_Heart$predicted.id)

# Unimodal Map Projection

# kalucka_Heart <- RunUMAP(kalucka_Heart, dims = 1:30, reduction = "pca")
Control_Heart <- MapQuery(anchorset = kalucka.Anchors, reference = kalucka_Heart, query = Control_Heart,
                          refdata = list(celltype = "Cluster"), reference.reduction = "pca", reduction.model = "umap")

p1 <- DimPlot(kalucka_Heart, reduction = "umap", group.by = "Cluster", label = TRUE, label.size = 3, label.box = T,
              repel = TRUE) + NoLegend() + ggtitle("Kalucka Heart")

p2 <- DimPlot(Control_Heart, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE, label.box = T,
              label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Control Heart 20Oct21")
p1+p2

saveRDS(Control_Heart, "./rds/Control_Heart.20Oct21.mapped.rds")

Control_Heart@meta.data$Condition <- Control_Heart@active.ident

Control_Heart@meta.data$FigClustering <- Control_Heart@meta.data$predicted.id

DimPlot(Control_Heart, reduction = "ref.umap",split.by = "Condition", group.by = "FigClustering")

saveRDS(Control_Heart, "./rds/Liver_20Jan21.Figures.All.Conditions.rds")
