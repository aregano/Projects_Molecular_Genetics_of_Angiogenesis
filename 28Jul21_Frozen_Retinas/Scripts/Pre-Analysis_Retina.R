# Starting analysis on the Liver sample

library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(dittoSeq)
library(patchwork)

Sample_28Jul21 <- readRDS("./rds/28Jul21.FrozenRetinas.ESCells.rds")

DimPlot(Sample_28Jul21)

# Loop to look at different Dimensions and plot a UMAP

# For Loop, going for 9 dim settings: 2,4,7,10,13,17,21,25,30

dim <- c(2, 4, 7, 10, 13, 17, 21, 25, 30)
n <- as.data.frame(list(c(1:9)))


myplots <- vector('list', nrow(n))

for (i in 1:length(dim)){
  FindNeighbors(Sample_28Jul21, dims = 1:dim[i])
  Sample_28Jul21 <- FindClusters(Sample_28Jul21, resolution = 0.3)
  Sample_28Jul21 <- RunUMAP(Sample_28Jul21, dims = 1:dim[i])
  myplots[[i]] <- local({
    i <- i
    p1 <- DimPlot(Sample_28Jul21, group.by = "Condition", label = TRUE, repel = T, label.box = TRUE) + ggtitle(paste(dim[i],"PCAs"))+ theme(plot.title = element_text(hjust = 0.5))+ NoAxes()+ NoLegend()
  })
  
}

p25 <- patchwork::wrap_plots(myplots, ncol = 3)

p25

for (i in 1:length(dim)){
  FindNeighbors(Sample_28Jul21, dims = 1:dim[i])
  Sample_28Jul21 <- FindClusters(Sample_28Jul21, resolution = 0.3)
  reduction <- paste("umap", dim[i], "dim")
  Sample_28Jul21 <- RunUMAP(Sample_28Jul21, dims = 1:dim[i], reduction.name = reduction)
  
}


#  We can see that below 7 dims is useless. The right dims must be between 16-21


# Now we should check the quality of the reads and features in the UMAP
FeaturePlot(Sample_28Jul21, reduction = "umap",
            features = c("nCount_RNA", "nFeature_RNA"))

# It does not look bad I think

# Get the HTOs in the meta.data so I can sort singlets by using that criteria

DefaultAssay(Sample_28Jul21) <- "ADT"

p1 <- FeaturePlot(Sample_28Jul21, c("B0305","B0306","B0307"),
                  order = T,ncol = 1,
                  cols = c("lightgrey", "darkgreen")) 

p1

# Demultiplexing HTOs

# If you have a very large dataset we suggest using k_function = 'clara'. This is a k-medoid
# clustering function for large applications You can also play with additional parameters (see
# documentation for HTODemux()) to adjust the threshold for classification Here we are using
# the default settings

Sample_28Jul21 <- HTODemux(Sample_28Jul21, assay = "ADT", positive.quantile = 0.99)

table(Sample_28Jul21$ADT_classification.global)

# Group cells based on the max HTO signal
Idents(Sample_28Jul21) <- "ADT_maxID"
RidgePlot(Sample_28Jul21, assay = "ADT", features = rownames(Sample_28Jul21[["ADT"]])[1:3], ncol = 1)

hto <- rownames(Sample_28Jul21)

p1 <- FeatureScatter(Sample_28Jul21, feature1 = "B0305", feature2 = "B0306")

p2 <- FeatureScatter(Sample_28Jul21, feature1 = "B0305", feature2 = "B0307")

p1/p2

Idents(Sample_28Jul21) <- "ADT_classification.global"
VlnPlot(Sample_28Jul21, features = "nCount_RNA", pt.size = 0.1, log = TRUE)


# First, we will remove negative cells from the object
Sample_28Jul21.subset <- subset(Sample_28Jul21, idents = "Negative", invert = TRUE)

# Calculate a tSNE embedding of the HTO data
DefaultAssay(Sample_28Jul21.subset) <- "ADT"
Sample_28Jul21.subset <- ScaleData(Sample_28Jul21.subset, features = rownames(Sample_28Jul21.subset),
                                 verbose = FALSE)
Sample_28Jul21.subset <- RunPCA(Sample_28Jul21.subset, features = rownames(Sample_28Jul21.subset), approx = FALSE)
Sample_28Jul21.subset <- RunTSNE(Sample_28Jul21.subset, dims = 1:8, perplexity = 100)
DimPlot(Sample_28Jul21.subset)

# To increase the efficiency of plotting, you can subsample cells using the num.cells argument
HTOHeatmap(Sample_28Jul21.subset, assay = "ADT")

hto <- rownames(Sample_28Jul21)

Idents(Sample_28Jul21.subset) <- "ADT_maxID"

p1 <- FeatureScatter(Sample_28Jul21.subset, feature1 = "B0305", feature2 = "B0306")

p2 <- FeatureScatter(Sample_28Jul21.subset, feature1 = "B0305", feature2 = "B0307")

p1/p2

RidgePlot(Sample_28Jul21.subset, assay = "ADT", features = rownames(Sample_28Jul21[["ADT"]])[1:3], ncol = 1)

# Extract the singlets
Idents(Sample_28Jul21.subset) <- "ADT_classification.global"

Retina <- subset(Sample_28Jul21.subset, idents = "Singlet")

DefaultAssay(Retina) <- "RNA"

# Select the top 1000 most variable features
Retina <- FindVariableFeatures(Retina, selection.method = "vst", nfeatures = 2000)

# Scaling RNA data, we only scale the variable features here for efficiency
Retina <- ScaleData(Retina, features = VariableFeatures(Retina))

# Run PCA
Retina <- RunPCA(Retina, features = VariableFeatures(Retina))

ElbowPlot(Retina)


# We select the top 10 PCs for clustering and tSNE based on PCElbowPlot
Retina <- FindNeighbors(Retina, reduction = "pca", dims = 1:14)
Retina <- FindClusters(Retina, resolution = 0.6, verbose = FALSE)
Retina <- RunUMAP(Retina, reduction = "pca", dims = 1:14)

# Projecting singlet identities on TSNE visualization
DimPlot(Retina, group.by = "ADT_classification")

#  Now to take out Blood cells

FeaturePlot(Retina,
            features = c("Ptprc", # Leukocyte lineage
                         "Cd14", # Macrophage
                         "Cd3e", # T cell                        
                         "Cd79a") # B cell
)

p1 <- DimPlot(Retina)

p1

Blood.cells <- CellSelector(p1, ident = "BloodCell")

# You can use a selector or clustering, or expression levels, to weed out blood cells 

Retina <- subset(Retina, cells = Blood.cells, invert = TRUE)

Retina <- subset(Retina, subset = Ptprc > 0 | Cd14 > 0 | Cd3e > 0 | Cd79a > 0, invert = TRUE)

FeaturePlot(Retina, features = c("Cdh5", "tdTomato-WPRE-sv40pA"))

Idents(Retina) <- "hash.ID"

DimPlot(Retina)

Condition <- c("Control", "Dll4LOF", "Dll4/MycLOF")

levels(Retina) -> names(Condition) 

Retina <- RenameIdents(Retina, Condition)

Retina@meta.data$Condition <- Retina@active.ident

# Lets look again at different PCAs and decide the best one

# For Loop, going for 9 dim settings: 2,4,7,10,13,17,21,25,30

dim <- c(4, 7, 10, 13, 17, 20, 23, 27, 30)
n <- as.data.frame(list(c(1:9)))


myplots <- vector('list', nrow(n))

for (i in 1:length(dim)){
  FindNeighbors(Retina, dims = 1:dim[i])
  Retina <- FindClusters(Retina, resolution = 0.3)
  Retina <- RunUMAP(Retina, dims = 1:dim[i])
  myplots[[i]] <- local({
    i <- i
    p1 <- DimPlot(Retina, group.by = "Condition") + ggtitle(paste(dim[i],"PCAs"))+ theme(plot.title = element_text(hjust = 0.5))
  })
  
}

patchwork::wrap_plots(myplots, ncol = 3)

# I will go for 19 dimensions since that is the same as the previous liver analysis, and the plots look all the same after 17 dimensions

ElbowPlot(Retina, ndims = 30)

# THe ElbowPlot would tell us that 15 dimensions is a good threshold

rui_palette <- c( "#D0C9D7", "#71718D", "#AAAA55", "#E2E21C", "#FFE200", "#FFAA00", "#FF7100", "#FF3800", "#FF0000")


FeaturePlot(Retina, features = c("Dll4", "Myc"), split.by = "Condition", cols = rui_palette, order = T)


DimPlot(Retina, reduction = "umap", label = T, label.box = T, label.size = 5)
DimPlot(Retina, reduction = "umap", split.by = "Condition", label.size = 5)
DimPlot(Retina, split.by = "hash.ID")


Retina <- FindVariableFeatures(Retina, 
                                       selection.method = "vst",
                                       nfeatures = 2000)
Retina <- RunPCA(Retina,
                         features = VariableFeatures(object = Retina))
Retina <- FindNeighbors(Retina, dims = 1:19)
Retina <- FindClusters(Retina, resolution = 0.1)
Retina <- FindClusters(Retina, resolution = 0.3)
Retina <- FindClusters(Retina, resolution = 0.5)

Retina <- RunUMAP(Retina, dims = 1:19)

# Now let's rename all HTO idents to each condition
# BC305: "Control", BC302: "Dll4KO+Vehicle", BC303: "Dll4KO+SL327", BC304: "Dll4/MycKO"


conditions_names <- c("Control", "Dll4KO", "Dll4/MycKO")

Idents(Retina) <- "hash.ID"

names(conditions_names) <- levels(Retina)
Retina <- RenameIdents(Retina, conditions_names)

Retina@meta.data$Condition <- Retina@active.ident

p1 <- DimPlot(Retina, split.by = "Condition")

p2 <- DimPlot(Retina, split.by = "hash.ID")

p1/p2

saveRDS(Retina, "./rds/FrozenRetina.ECs.Alvaro.v2.rds")

# Now onto plotting the transgenes present in the retina dataset

Bestholtz_palette <- c("#DEDAD6", "#FEE392", "#FEC44E", "#FE9929", "#ED6F11", "#CC4C17", "#993411", "#65260C")

rui_palette <- c( "#D0C9D7", "#71718D", "#AAAA55", "#E2E21C", "#FFE200", "#FFAA00", "#FF7100", "#FF3800", "#FF0000")

genes <-  list(c("WPRE", "Pdgfb-CreERT2-IRES-EGFP", "Cdh5-CreERT2-BGH-pA", "sv40pA", "PhiYFP", "Cdh5"))

genes <- list(c("Odc1", "Kcne3", "Mki67", "Wnt2", "Msr1", "Ltbp4", "Rspo3", "Gja5", "Cdkn1a", "Stmn1", "Vegfa", "Esm1", "Hes1", "Dll4", "Myc", "Mycn"))

genes <- as.data.frame(genes)


#  Script

myplots <- vector('list', nrow(genes))


for (i in 1:nrow(genes)) {
  
  p21 <- FeaturePlot(Retina, pt.size = 1.3, features = genes[i, 1], combine = F)
  
  p22 <- FeaturePlot(Retina, split.by = "Condition", features = genes[i, 1], combine = F)
  
  p21 <- lapply(X = p21, FUN = function(p) p + NoAxes()+  scale_colour_gradientn(colours = rui_palette))
  
  p22 <- lapply(X = p22, FUN = function(p) p + NoLegend() + NoAxes()+  scale_colour_gradientn(colours = rui_palette))
  
  p21 <- Reduce( `+`, p21 )+patchwork::plot_layout( ncol = 1 )
  
  p22 <- Reduce( `+`, p22 )+patchwork::plot_layout( ncol = 3 )
  
  p23 <- list(p21, p22)
  
  
  design <- c(patchwork::area(1, 1, 1, 1), patchwork::area(1, 2, 1, 4))
  
  p24 <- Reduce( `+`,  p23)+patchwork::plot_layout(design = design)
  
  
  myplots[[i]] <- local({
    i <- i
    p24
  })
  
  
}

p25 <- patchwork::wrap_plots(myplots, ncol = 1)

#myplots[1:nrow(genes)] + patchwork::plot_layout(byrow = T, widths = 25, heights = 5)
p25


#  Same with Violin plots


for (i in 1:nrow(genes)) {
  
  p21 <- VlnPlot(Retina, features = genes[i, 1], split.by = "Condition") + NoLegend()+  theme(axis.title.x = element_blank(), axis.text.x = element_blank())
  
  myplots[[i]] <- local({
    i <- i
    p21
  })
  
  plegend <- VlnPlot(Retina, features = genes[i, 1], split.by = "Condition")
}

legend <- get_legend(plegend + theme(legend.position = "bottom", legend.justification = "center"))

p25 <- patchwork::wrap_plots(myplots, ncol = 4)
p25
#myplots[1:nrow(genes)] + patchwork::plot_layout(byrow = T, widths = 25, heights = 5)
p26 <- plot_grid(p25, legend, ncol = 1, rel_heights = c(1, .1))

p26


# Trying out dotplots

genes <- list(c("Odc1", "Kcne3"))


genes <- list(c("Dll4", "Myc", "Mycn"))

genes <- as.data.frame(genes)

p1 <- DotPlot(Retina, features = genes, dot.scale = 20, scale.by = "radius", scale.max = 100, scale.min = 0, scale = F)+ scale_colour_gradientn(colours = rui_palette)+ 
  theme(axis.title.x.bottom = element_blank(), axis.title.y.left = element_blank(), plot.subtitle = element_blank())

p2 <- DotPlot(Retina, features = genes, group.by = "Condition", dot.scale = 20, scale.by = "radius", scale.max = 25, scale.min = 0, scale = F)+ scale_colour_gradientn(colours = rui_palette)

p1/p2
