# Link: https://bioconductor.org/packages/release/bioc/vignettes/ReactomeGSA/inst/doc/analysing-scRNAseq.html

library(ReactomeGSA)
library(ReactomeGSA.data)
library(Seurat)
library(dplyr)
library(openxlsx)

data(jerby_b_cells)

Liver <- readRDS("//Tierra/SC/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Liver_MacarenaFC/Figures_Paper/rds/Liver_20Jan21.Figures.All.Conditions.OldDll4Myc.rds")

Liver <- subset(Liver, subset = Condition == "Control" | Condition == "Dll4KO" | Condition == "NOTCH1KO" | Condition == "RBPJKO")

levels(Liver) <- c("Control", "Dll4KO", "NOTCH1KO", "RBPJKO")

Liver@active.ident -> Liver@meta.data$Condition

Idents(Liver) <- "FigClustering"

gsva_result <- analyse_sc_clusters(Liver, verbose = TRUE)

# The resulting object is a standard ReactomeAnalysisResult object.

gsva_result

# pathways returns the pathway-level expression values per cell cluster:

pathway_expression <- pathways(gsva_result)

class(pathway_expression)

# pathway_expression_ls <- as.list(pathway_expression)

# simplify the column names by removing the default dataset identifier
colnames(pathway_expression) <- gsub("\\.Seurat", "", colnames(pathway_expression))

pathway_expression[1:3,]

#  A simple approach to find the most relevant pathways is to assess the maximum difference in expression for every pathway:
class(pathway_expression)

# find the maximum differently expressed pathway
max_difference <- do.call(rbind, apply(pathway_expression, 1, function(row) {
  values <- as.numeric(row[2:length(row)])
  return(data.frame(name = row[1], min = min(values), max = max(values)))
}))

pathway_expression$C6__X

pathway_expression_C6 <- pathway_expression %>% arrange(desc(C6__X)) 

# find the maximum differently expressed pathway
max_difference <- do.call(rbind, apply(pathway_expression, 1, function(row) {
  values <- as.numeric(row[2:length(row)])
  C6 <- as.numeric(row[length(row)])
  C5p <- as.numeric(row[length(row)-1])
  return(data.frame(name = row[1], min = min(values), max = max(values), C5p = C5p, C6 = C6))
}))

max_difference$diff <- max_difference$max - max_difference$min

max_difference$diffC6 <- max_difference$C6 - max_difference$min

# sort based on the difference
max_difference <- max_difference[order(max_difference$diff, decreasing = T), ]

max_difference <- max_difference[order(max_difference$diff, decreasing = T), ]

# sort based on the difference
max_difference <- max_difference[order(max_difference$diffC6, decreasing = T), ]

head(max_difference)

# sort based on the difference of C6

max_difference_C6 <- filter(max_difference, max == C6)

max_difference_C6 <- max_difference_C6[order(max_difference_C6$diffC6, decreasing = T), ]

head(max_difference_C6)

# sort based on the difference for C6 cluster


pathways_C6 <- row.names(max_difference_C6)

pathways_C6[1:10]

# Plotting the Results

max_difference_C5p <- filter(max_difference, max == C5p)

max_difference_C5p$diff <- max_difference_C5p$max - max_difference_C5p$min

max_difference_C5p <- max_difference_C5p[order(max_difference_C5p$diff, decreasing = T), ]

head(max_difference_C5p)

plot_gsva_pathway(gsva_result, pathway_id = rownames(max_difference_C5p)[1])

# Additional parameters are directly passed to gplots heatmap.2 function
plot_gsva_heatmap(gsva_result, max_pathways = 15, margins = c(6,20))

relevant_pathways <- c("R-HSA-983170", "R-HSA-388841", "R-HSA-2132295", "R-HSA-983705", "R-HSA-5690714")
plot_gsva_heatmap(gsva_result, 
                  pathway_ids = pathways_C6[1:10], # limit to these pathways
                  margins = c(15,30), # adapt the figure margins in heatmap.2
                  dendrogram = "col", # only plot column dendrogram
                  scale = "row", # scale for each pathway
                  key = FALSE, # don't display the color key
                  lwid=c(0.1,4)) # remove the white space on the left


plot_gsva_pca(gsva_result)

write.xlsx(max_difference_C5p, "./Tables/Reactome.C5.pathways.xlsx")

write.xlsx(max_difference_C6, "./Tables/Reactome.C6.pathways.xlsx")
write.xlsx(max_difference, "./Tables/Reactome.pathways.xlsx")
