library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(BiocStyle)
library(dittoSeq)
library(dplyr)

Liver_Heart.data <- Read10X(data.dir = "raw_data/")

Liver_Heart <- CreateSeuratObject(counts = Liver_Heart.data$`Gene Expression`, 
                                          project = "Liver_Heart20Jan22", assay = "RNA",
                                          min.cells = 3 , min.features = 200)

Liver_Heart[["HTO"]] <- CreateAssayObject(counts = Liver_Heart.data$`Antibody Capture`[, colnames(x = Liver_Heart)])

Liver_Heart[["CMO"]] <- CreateAssayObject(counts = Liver_Heart.data$`Multiplexing Capture`[, colnames(x = Liver_Heart)])



rm(Liver_Heart.data)

saveRDS(Liver_Heart, "../Sample_Liver_Heart_20Jan22/rds/Liver_Heart_raw.rds")

DefaultAssay(Liver_Heart) <- "HTO"
DefaultAssay(Liver_Heart)
rownames(Liver_Heart)

# Cluster cells in terms of GEx count data

DefaultAssay(Liver_Heart) <- "RNA"
DefaultAssay(Liver_Heart)

Liver_Heart <- NormalizeData(Liver_Heart)
Liver_Heart <- FindVariableFeatures(Liver_Heart)
Liver_Heart <- ScaleData(Liver_Heart)
Liver_Heart <- RunPCA(Liver_Heart, verbose = FALSE)
Liver_Heart <- FindNeighbors(Liver_Heart, dims = 1:30)
Liver_Heart <- FindClusters(Liver_Heart, resolution = 0.8, verbose = FALSE)
Liver_Heart <- RunUMAP(Liver_Heart, dims = 1:30)
DimPlot(Liver_Heart, label = TRUE)

#  Adding HTO data as an independent assay

Liver_Heart <- NormalizeData(Liver_Heart, assay = "HTO", normalization.method = "CLR")

Liver_Heart <- HTODemux(Liver_Heart, assay = "HTO", positive.quantile = 0.99)

table(Liver_Heart$HTO_classification.global)

# Group cells based on the max HTO signal
Idents(Liver_Heart) <- "HTO_maxID"
RidgePlot(Liver_Heart, assay = "HTO", features = rownames(Liver_Heart[["HTO"]])[1:3], ncol = 1)

FeatureScatter(Liver_Heart, feature1 = "HTO305", feature2 = "HTO307")

Idents(Liver_Heart) <- "HTO_classification.global"
VlnPlot(Liver_Heart, features = "nCount_RNA", pt.size = 0.1, log = TRUE)

# First, we will remove negative cells from the object
Liver_Heart.subset <- subset(Liver_Heart, idents = "Negative", invert = TRUE)

# Calculate a tSNE embedding of the HTO data
DefaultAssay(Liver_Heart.subset) <- "HTO"
Liver_Heart.subset <- ScaleData(Liver_Heart.subset, features = rownames(Liver_Heart.subset),
                                 verbose = FALSE)
Liver_Heart.subset <- RunPCA(Liver_Heart.subset, features = rownames(Liver_Heart.subset), approx = FALSE)
Liver_Heart.subset <- RunTSNE(Liver_Heart.subset, dims = 1:3, perplexity = 100)
DimPlot(Liver_Heart.subset)

# To increase the efficiency of plotting, you can subsample cells using the num.cells argument
HTOHeatmap(Liver_Heart, assay = "HTO", ncells = 30000)


# Extract the singlets
Liver_Heart.singlet <- subset(Liver_Heart, idents = "Singlet")

# Select the top 1000 most variable features
Liver_Heart.singlet <- NormalizeData(Liver_Heart.singlet)
Liver_Heart.singlet <- FindVariableFeatures(Liver_Heart.singlet, selection.method = "mean.var.plot")

# Scaling RNA data, we only scale the variable features here for efficiency
Liver_Heart.singlet <- ScaleData(Liver_Heart.singlet, features = VariableFeatures(Liver_Heart.singlet))

# Run PCA
Liver_Heart.singlet <- RunPCA(Liver_Heart.singlet, features = VariableFeatures(Liver_Heart.singlet))

# We select the top 10 PCs for clustering and tSNE based on PCElbowPlot
Liver_Heart.singlet <- FindNeighbors(Liver_Heart.singlet, reduction = "pca", dims = 1:10)
Liver_Heart.singlet <- FindClusters(Liver_Heart.singlet, resolution = 0.6, verbose = FALSE)
Liver_Heart.singlet <- RunTSNE(Liver_Heart.singlet, reduction = "pca", dims = 1:10)

# Projecting singlet identities on TSNE visualization
DimPlot(Liver_Heart.singlet, group.by = "HTO_classification")



# QC Analysis
DefaultAssay(Liver_Heart.Singlets) <- "RNA"
Liver_Heart.Singlets[["percent.mt"]] <- PercentageFeatureSet(Liver_Heart.Singlets, pattern = "^mt-")
head(Liver_Heart.Singlets@meta.data, 50)

VlnPlot(Liver_Heart.Singlets, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        group.by = "orig.ident", ncol = 3, pt.size = 0)

plot1 <- FeatureScatter(Liver_Heart.Singlets, feature1 = "nFeature_RNA",
                        group.by = "orig.ident", feature2 = "percent.mt") +
  geom_vline(xintercept = c(200,5000),linetype = 2 ) +
  geom_hline(yintercept = 15 ,linetype = 2)
plot2 <- FeatureScatter(Liver_Heart.Singlets, feature1 = "nCount_RNA",
                        group.by = "orig.ident",feature2 = "nFeature_RNA")+
  geom_hline(yintercept = c(200,5000),linetype = 2 )
plot1 / plot2

density(Liver_Heart.Singlets@meta.data$nFeature_RNA)

DenPlot(Liver_Heart.Singlets, x, group.by, title = NULL, subtitle = NULL,
        show.mean = FALSE)

nFeature_RNA <- Liver_Heart.Singlets@meta.data$nFeature_RNA

nFeature_RNA <- as.data.frame(nFeature_RNA)

# Produce Density plot like Carlos

nFeature_RNA %>%
  ggplot(aes(x=nFeature_RNA)) +
  geom_density(fill='#C8DCE7', size = 1, alpha = 0.7) +
  scale_x_continuous() +
  ggtitle("Genes Detected Density Plot")+
  ylab("density") +
  xlab("detected")+
  geom_vline(xintercept = 600, color = "red", size = 1)

# Produce Cummulative plot like Carlos

ggplot(nFeature_RNA, aes(nFeature_RNA)) + stat_ecdf(geom = "point", color = "#F8756C")+
  scale_x_continuous() +
  ggtitle("Genes Detected Cummulative Plot")+
  ylab("CumSum_Features") +
  xlab("total_Features")+
  geom_vline(xintercept = 600, color = "red", size = 1)

# Going over mitochondrial genes

percent_mt <- Liver_Heart.Singlets@meta.data$percent.mt
percent_mt <- as.data.frame(percent_mt)

percent_mt %>%
  ggplot(aes(x=percent_mt)) +
  geom_density(fill='#C8DCE7', size = 1, alpha = 0.7) +
  scale_x_continuous() +
  scale_y_discrete()+
  ggtitle("MT Content Density Plot")+
  ylab("density") +
  xlab("subset_MT_percent")+
  ylim(-0.4, 0.4)+
  geom_vline(xintercept = 25, color = "red", size = 1)+
  geom_jitter(aes(x = percent_mt, y = 0), height = 0.1, alpha = 0.01)


# Name cells in the conditions category

Liver_Heart.singlet -> Liver_Heart.Singlets

Idents(Liver_Heart.Singlets) <- Liver_Heart.Singlets@meta.data$hash.ID 

Liver_Heart.Singlets <- RenameIdents(Liver_Heart.Singlets, 'HTO305' = 'Misc', 'HTO302' = 'Dll4/MycLOF_Frozen', 'HTO307' = 'Dll4/MycLOF')

Liver_Heart.Singlets@meta.data$Condition <- Liver_Heart.Singlets@active.ident

# Moving onto CMO Demultiplexing. 

FeaturePlot(Liver_Heart.Singlets, 
            c("CMO301","CMO302","CMO303","CMO304","CMO305"),ncol = 2,
            order = T)

#  Interestingly CMO309 shows a number of reads, which makes no sense
FeaturePlot(Liver_Heart.Singlets, 
            c("CMO306","CMO307","CMO308","CMO309","CMO310","CMO311","CMO312"),ncol = 2,
            order = T)

# It seems I need to go to the specific folders for each CMO section in the cellranger multi output

CMO301 <- Read10X(data.dir = "raw_data/per_sample_outs/Control_Non_ECs_Heart/count/sample_feature_bc_matrix/")

CMO_301_Seurat <- CreateSeuratObject(CMO301$`Gene Expression`, project = "Liver_Heart20Jan22", assay = "RNA", min.cells = 3 , min.features = 200 )

CMO_301_Seurat[["HTO"]] <- CreateAssayObject(counts = CMO301$`Antibody Capture`[, colnames(x = CMO_301_Seurat)])

CMO_301_Seurat[["CMO"]] <- CreateAssayObject(counts = CMO301$`Multiplexing Capture`[, colnames(x = CMO_301_Seurat)])

CMO_301_Seurat <- NormalizeData(CMO_301_Seurat)
CMO_301_Seurat <- FindVariableFeatures(CMO_301_Seurat)
CMO_301_Seurat <- ScaleData(CMO_301_Seurat)
CMO_301_Seurat <- RunPCA(CMO_301_Seurat, verbose = FALSE)
CMO_301_Seurat <- FindNeighbors(CMO_301_Seurat, dims = 1:30)
CMO_301_Seurat <- FindClusters(CMO_301_Seurat, resolution = 0.8, verbose = FALSE)
CMO_301_Seurat <- RunUMAP(CMO_301_Seurat, dims = 1:30)
DimPlot(CMO_301_Seurat, label = TRUE) 

DefaultAssay(Liver_Heart.Singlets) <- "CMO"
Liver_Heart.Singlets <- NormalizeData(Liver_Heart.Singlets, assay = "CMO", normalization.method = "CLR")

Liver_Heart.Singlets <- HTODemux(Liver_Heart.Singlets, assay = "CMO", positive.quantile = 0.99)



# Selecting specific subset of cells

Dll4_MycLOF_Fresh <- subset(Liver_Heart.Singlets, subset = hash.ID == "HTO307-TotalSeqB")

levels(Dll4_MycLOF_Fresh) <- "Dll4/MycLOF"
Idents(Dll4_MycLOF_Fresh) <- "Dll4/MycLOF"
Dll4_MycLOF_Fresh@meta.data$Condition <- Dll4_MycLOF_Fresh@active.ident


saveRDS(Liver_Heart.Singlets, "./rds/Liver_Heart.Singlets.rds")