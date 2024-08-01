library(Seurat)
DefaultAssay(Sample_28Jul21) <- "ADT"


DefaultAssay(Sample_28Jul21) <- "RNA"

# Select the top 1000 most variable features
Sample_28Jul21 <- FindVariableFeatures(Sample_28Jul21, selection.method = "vst", nfeatures = 2000)

variable <- FindVariableFeatures(Sample_28Jul21, selection.method = "vst", nfeatures = 2000)

# Scaling RNA data, we only scale the variable features here for efficiency
Sample_28Jul21 <- ScaleData(Sample_28Jul21, features = VariableFeatures(Sample_28Jul21))

# Run PCA
Sample_28Jul21 <- RunPCA(Sample_28Jul21, features = VariableFeatures(Sample_28Jul21))

ElbowPlot(Sample_28Jul21)


# We select the top 10 PCs for clustering and tSNE based on PCElbowPlot
Sample_28Jul21 <- FindNeighbors(Sample_28Jul21, reduction = "pca", dims = 1:14)
Sample_28Jul21 <- FindClusters(Sample_28Jul21, resolution = 0.6, verbose = FALSE)
Sample_28Jul21 <- RunUMAP(Sample_28Jul21, reduction = "pca", dims = 1:14)

# Projecting singlet identities on TSNE visualization
DimPlot(Sample_28Jul21)

#  Now to take out Blood cells

FeaturePlot(Sample_28Jul21,
            features = c("Ptprc", # Leukocyte lineage
                         "Cd14", # Macrophage
                         "Cd3e", # T cell                        
                         "Cd79a") # B cell
)

p1 <- DimPlot(Sample_28Jul21)

p1

Blood.cells <- CellSelector(p1, ident = "BloodCell")

# You can use a selector or clustering, or expression levels, to weed out blood cells 

Sample_28Jul21 <- subset(Sample_28Jul21, cells = Blood.cells, invert = TRUE)

Sample_28Jul21 <- subset(Sample_28Jul21, subset = Ptprc > 0 | Cd14 > 0 | Cd3e > 0 | Cd79a > 0, invert = TRUE)

FeaturePlot(Sample_28Jul21, features = c("Cdh5", "tdTomato-WPRE-sv40pA"))

Idents(Sample_28Jul21) <- "hash.ID"

DimPlot(Sample_28Jul21)

Condition <- c("Control", "Dll4LOF", "Dll4/MycLOF")

levels(Sample_28Jul21) -> names(Condition) 

Sample_28Jul21 <- RenameIdents(Sample_28Jul21, Condition)

Sample_28Jul21@meta.data$Condition <- Sample_28Jul21@active.ident

# Lets look again at different PCAs and decide the best one

# For Loop, going for 9 dim settings: 2,4,7,10,13,17,21,25,30

dim <- c(4, 7, 10, 13, 17, 20, 23, 27, 30)
n <- as.data.frame(list(c(1:9)))


myplots <- vector('list', nrow(n))

for (i in 1:length(dim)){
  FindNeighbors(Sample_28Jul21, dims = 1:dim[i])
  Sample_28Jul21 <- FindClusters(Sample_28Jul21, resolution = 0.3)
  Sample_28Jul21 <- RunUMAP(Sample_28Jul21, dims = 1:dim[i])
  myplots[[i]] <- local({
    i <- i
    p1 <- DimPlot(Sample_28Jul21, group.by = "Condition") + ggtitle(paste(dim[i],"PCAs"))+ theme(plot.title = element_text(hjust = 0.5))
  })
  
}

patchwork::wrap_plots(myplots, ncol = 3)

# I will go for 19 dimensions since that is the same as the previous liver analysis, and the plots look all the same after 17 dimensions

ElbowPlot(Sample_28Jul21, ndims = 30)

# THe ElbowPlot would tell us that 15 dimensions is a good threshold

rui_palette <- c( "#D0C9D7", "#71718D", "#AAAA55", "#E2E21C", "#FFE200", "#FFAA00", "#FF7100", "#FF3800", "#FF0000")


FeaturePlot(Sample_28Jul21, features = c("Dll4", "Myc"), split.by = "Condition", cols = rui_palette, order = T)


DimPlot(Sample_28Jul21, reduction = "umap", label = T, label.box = T, label.size = 5)
DimPlot(Sample_28Jul21, reduction = "umap", split.by = "Condition", label.size = 5)
DimPlot(Sample_28Jul21, split.by = "hash.ID")


Sample_28Jul21 <- FindVariableFeatures(Sample_28Jul21, 
                               selection.method = "vst",
                               nfeatures = 2000)
Sample_28Jul21 <- RunPCA(Sample_28Jul21,
                 features = VariableFeatures(object = Sample_28Jul21))
Sample_28Jul21 <- FindNeighbors(Sample_28Jul21, dims = 1:19)
Sample_28Jul21 <- FindClusters(Sample_28Jul21, resolution = 0.1)
Sample_28Jul21 <- FindClusters(Sample_28Jul21, resolution = 0.3)
Sample_28Jul21 <- FindClusters(Sample_28Jul21, resolution = 0.5)

Sample_28Jul21 <- RunUMAP(Sample_28Jul21, dims = 1:19)

# Now let's rename all HTO idents to each condition
# BC305: "Control", BC302: "Dll4KO+Vehicle", BC303: "Dll4KO+SL327", BC304: "Dll4/MycKO"


conditions_names <- c("Dll4/MycKO", "Dll4KO(4d)+Vehicle", "Control(4d)", "Dll4KO(4d)+SL327")

names(conditions_names) <- levels(Sample_28Jul21)
Sample_28Jul21 <- RenameIdents(Sample_28Jul21, conditions_names)

Sample_28Jul21@meta.data$Condition <- Sample_28Jul21@active.ident

saveRDS(Sample_28Jul21, "./rds/FrozenRetina.ECs.Alvaro.rds")
