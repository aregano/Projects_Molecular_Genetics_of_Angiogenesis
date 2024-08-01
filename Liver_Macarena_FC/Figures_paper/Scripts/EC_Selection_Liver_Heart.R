Liver_Heart <- subset(Liver_Heart, subset = Condition == "Dll4/MycLOF")

FeaturePlot(Liver_Heart, features = "WPRE-SV40PA")

ElbowPlot(Liver_Heart, ndims = 30)

table(Liver_Heart@meta.data$Condition)

Liver_Heart <- FindNeighbors(Liver_Heart, dims = 1:30)
Liver_Heart <- FindClusters(Liver_Heart, resolution = 0.1)
Liver_Heart <- FindClusters(Liver_Heart, resolution = 0.2)
Liver_Heart <- FindClusters(Liver_Heart, resolution = 0.3)

Liver_Heart <- RunUMAP(Liver_Heart, dims = 1:30)
DimPlot(Liver_Heart, reduction = "umap", group.by = "RNA_snn_res.0.3", label = T, label.size = 5)

DimPlot(Liver_Heart, group.by = "Condition")

FeaturePlot(Liver_Heart, features = "CDH5")

FeaturePlot(Liver_Heart,
            features = c("PTPRC", # Leukocyte lineage
                         "CD14", # Macrophage
                         "CD3E", # T cell                        
                         "CD79A") # B cell
)

VlnPlot(Liver_Heart, ncol = 2,
        features =  c("PTPRC",
                      "CD14", # Macrophage
                      "CD3E", # T cell                        
                      "CD79A"), # B cell
        group.by =  "RNA_snn_res.0.3", 
        pt.size = 0 ) + theme(legend.position="none")

x <- Liver_Heart@meta.data$RNA_snn_res.0.3  %in%  as.character(c(0:4,6))
p1 <- DimPlot(Liver_Heart, reduction = "umap", cells = x,
              group.by = "RNA_snn_res.0.3",
              label = T, pt.size = 0.5)
p2 <- DimPlot(Liver_Heart,   cells = x,group.by="DF.classifications_0.25_0.09_913",
              cols = c("red","grey"),
              reduction="umap", pt.size=0.3 )
p1 + p2


x <- Liver_Heart@meta.data$RNA_snn_res.0.3  %in%  as.character(c(0:4,6))
Dll4MycLOF <- subset(Liver_Heart, cells = colnames(Liver_Heart)[x])
Dll4MycLOF <- FindVariableFeatures(Dll4MycLOF, 
                                                     selection.method = "vst",
                                                     nfeatures = 2000)
Dll4MycLOF <- RunPCA(Dll4MycLOF,
                                       features = VariableFeatures(object = Dll4MycLOF))
Dll4MycLOF <- FindNeighbors(Dll4MycLOF, dims = 1:30)
Dll4MycLOF <- FindClusters(Dll4MycLOF, resolution = 0.1)
                                             
                                             
Dll4MycLOF <- FindClusters(Dll4MycLOF, resolution = 0.05)        

genes <- rownames(Liver_Heart)
genes <- as.data.frame(genes)
Dll4MycLOF <- RunUMAP(Dll4MycLOF, dims = 1:30)

DimPlot(Dll4MycLOF)

FeaturePlot(Dll4MycLOF, features = "MYC")

VlnPlot(Dll4MycLOF, features = "MYC")


saveRDS(Dll4MycLOF, "rds/Groups/Dll4MycLOF.Fresh.rds")
