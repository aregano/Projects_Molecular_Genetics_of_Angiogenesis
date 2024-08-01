# Load merged samples and find old Myc Dataset

Liver <- readRDS("S:/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Liver_MacarenaFC/Figures_Paper/rds/Groups/Liver.Figures.merge.All.Conditions.rds")
Liver_2wdeletion.v5.ECs <- readRDS("S:/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Liver_MacarenaFC/Figures_Paper/rds/Liver_2wdeletion.v5.ECs.rds")

###################################################################################################################

# Color Palette

my_palette_Rui_colors <- c("#E95A74", "#50B6EF", "#45FF8E", "#F4A753", "#A80519", "#880088", "#E28CF4", "#C1B80C","#FC0808",  "#0E47D8")

my_palette_Rui_colors_B <- c("#F4A753", "#E95A74", "#50B6EF", "#FC0808",  "#0E47D8", "#45FF8E", "#A80519", "#E28CF4", "#880088", "#C1B80C")

custom_rdylblu_palette <- c("#005AD9", "#0070FF", "#FFE135", "#FFDA13", "#E60026", "#C60012")

macarena_palette <- c("#FDDC85", "#6E2914")

Bestholtz_palette <- c("#DEDAD6", "#FEE392", "#FEC44E", "#FE9929", "#ED6F11", "#CC4C17", "#993411", "#65260C")

palette_Maca_Vln <- c("#606060", "#FFA040", "#F94040", "#05BE78", "darkturquoise")


##################################################################################################################

table(Liver_2wdeletion.v5.ECs@meta.data$Condition)


Dll4MycOld <- subset(Liver_2wdeletion.v5.ECs, subset = Condition == "Dll4/MycKO")

condition <- c("Dll4/MycKO")

Idents(Dll4MycOld) <- "Condition"

table(Dll4MycOld@meta.data$Condition)

table(Liver@meta.data$Condition)

Liver1 <- subset(Liver, subset = Condition == "Control" | Condition == "Control(4d)" | Condition == "Dll4KO" | Condition == "Dll4/MycLOF")

Liver <- merge(Liver1, Dll4MycOld)

table(Liver@meta.data$Condition)

# Perform mapping

saveRDS(Liver.query, "rds/Myc.Testing.rds")

Liver <- readRDS("S:/LAB_RB/LAB/Alvaro/Bioinformatics/Analysis/Liver_MacarenaFC/Figures_Paper/rds/Myc.Testing.rds")

########################################################################################################################################################

# Figure UMAP

table(Liver@meta.data$Condition)

DimPlot(Liver, label = T, reduction = "ref.umap", group.by = "Condition")

levels(Liver) <- c("Control", "Control(4d)", "Dll4KO", "Dll4/MycKO", "Dll4/MycLOF")

Liver@active.ident -> Liver@meta.data$Condition

Idents(Liver) <- Liver@meta.data$FigClustering

p1 <- DimPlot(Liver, cols = my_palette_Rui_colors_B, pt.size = 1.2, reduction = "ref.umap", group.by = "FigClustering", split.by = "Condition", combine = T)

p1

new_labels <- c("Control" = "Control", "Control(4d)" = "Control4d", "Dll4KO" = "italic(Dll4)^iDEC", "Dll4/MycKO" = "italic(Dll4/Myc)^iDEC-Old", "Dll4/MycLOF" = "italic(Dll4/Myc)^iDEC-New")

p21 <- p1+facet_grid(~Condition, labeller = as_labeller(new_labels, default = label_parsed))+
  theme(strip.text.x = element_text(size = 24), plot.title = element_blank())+NoLegend()+xlab("UMAP_1")+ylab("UMAP_2")

p22 <- DimPlot(Liver, group.by = "FigClustering", reduction = "ref.umap", cols = my_palette_Rui_colors_B, pt.size = 1.2)+theme(legend.text = element_text(size = 14), axis.title.x = element_blank(), axis.line.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks = element_blank(), plot.title = element_blank())

p21
p22

p23 <- list(p21, p22)

design <- c(patchwork::area(1, 1, 1, 4), patchwork::area(1, 5, 1, 5.5))

p24 <- Reduce( `+`,  p23)+patchwork::plot_layout(design = design)


cairo_pdf("Plots/Dll4Myc_testing/Dll4-Myc_testing_UMAP.pdf",  width = 35, height = 6, family = "Arial")
p24
dev.off()


jpeg("Plots/Dll4Myc_testing/Dll4-Myc_testing_UMAP.jpeg", width = 35, height = 6, units = 'in', res = 800)
p24
dev.off()



###############################################################################################################################


# Bar Plot

table(Liver@meta.data$Condition)

p1 <- dittoBarPlot(Liver, var = "FigClustering", group.by = "Condition",  x.reorder = c(1,2,5,3,4), color.panel = my_palette_Rui_colors_B, xlab = NULL, main = NULL, y.breaks = c(0,0.25,0.5, 0.75, 1), min = 0, max = 1)+
  theme(axis.text=element_text(size=12), axis.title = element_text(size = 18), legend.text = element_text(size = 16), axis.text.x = element_text(angle = 25))+
  scale_x_discrete(labels=c("Control", "Control4d", bquote(italic(Dll4)^iDEC+Veh), bquote(italic(Dll4/Myc)^iDEC-Old), bquote(italic(Dll4/Myc)^iDEC-New)))+
  geom_col(width = 0.1)+
  # scale_y_discrete(labels=c("0", "25", "50", "75", "100"))
  scale_y_continuous(labels = c("0", "25", "50", "75", "100"))
p1

cairo_pdf("Plots/Dll4Myc_testing/Dll4-Myc_testing_BarPlot.pdf",  width = 8, height = 6, family = "Arial")
p1
dev.off()


jpeg("Plots/Dll4Myc_testing/Dll4-Myc_testing_BarPlot.jpeg", width = 8, height = 6, units = 'in', res = 800)
p1
dev.off()


###############################################################################################################################

genes <-rownames(Liver)

genes <- as.data.frame(genes)
# Figure 4N. Violin Plot

genes <- list(c("MYC", "DLL4", "ODC1", "WPRE-SV40PA"))
genes <- as.data.frame(genes)
gene.names <- list(c("Myc", "Dll4", "Odc1", "WPRE-Sv40pA"))
gene.names <- as.data.frame(gene.names)

#  Script

myplots <- vector('list', nrow(genes))


for (i in 1:nrow(genes)) {
  
  p21 <- VlnPlot(Liver, features = genes[i, 1], group.by = "Condition", cols = palette_Maca_Vln)+ NoLegend() + theme(text = element_text(family="Arial"), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.title.y=element_blank(), plot.title = element_text(face = "bold.italic"))+ggtitle(gene.names[i, 1])
  
  myplots[[i]] <- local({
    i <- i
    p21
  })
  
}

plegend <- VlnPlot(Liver, features = genes[i, 1], group.by = "Condition") +
  scale_fill_manual(values= palette_Maca_Vln, labels=c("Control", "Control4d", expression(italic("Dll4"^"iDEC")), expression(italic(Dll4/Myc)^"iDEC Old"), expression(italic(Dll4/Myc)^"iDEC New")))
plegend
legend <- get_legend(plegend + theme(legend.position = "bottom", legend.justification = "center", legend.direction = "horizontal", legend.text = element_text(size = 9)))

p25 <- patchwork::wrap_plots(myplots, ncol = 1)

p26 <- plot_grid(p25, legend, ncol = 1, rel_heights = c(1, .07))

x.axis <- ggdraw() + draw_label("Expression Level", fontface='bold', angle = 90, size = 24)

p27 <- plot_grid(x.axis, p26, ncol = 2, rel_widths = c(.1, 1))

cairo_pdf("Plots/Dll4Myc_testing/Dll4-Myc_testing_VlnPlot.pdf", width = 6, height = 12, family = "Arial")
p27
dev.off()

jpeg("Plots/Dll4Myc_testing/Dll4-Myc_testing_VlnPlot.jpeg", width = 6, height = 9, units = 'in', res = 800)
p27
dev.off()

##########################################################################

#  Count cells

mycells <- sum(GetAssayData(object = Liver, slot = "data")["ODC1",]>0)/nrow(Liver@meta.data$Condition)

sum(GetAssayData(object = Liver, slot = "data")["ODC1",]>0)/nrow(Liver@meta.data$Condition["Dll4KO"])

Liver@assays[["RNA"]]@data["ODC1",]

# Function

PrctCellExpringGene <- function(object, genes, group.by = "all"){
  if(group.by == "all"){
    prct = unlist(lapply(genes,calc_helper, object=object))
    result = data.frame(Markers = genes, Cell_proportion = prct)
    return(result)
  }
  
  else{        
    list = SplitObject(object, group.by)
    factors = names(list)
    
    results = lapply(list, PrctCellExpringGene, genes=genes)
    for(i in 1:length(factors)){
      results[[i]]$Feature = factors[i]
    }
    combined = do.call("rbind", results)
    return(combined)
  }
}

calc_helper <- function(object,genes){
  counts = object[['RNA']]@counts
  ncells = ncol(counts)
  if(genes %in% row.names(counts)){
    sum(counts[genes,]>0)/ncells
  }else{return(NA)}
}


GOI <- PrctCellExpringGene(Liver, c("ODC1", "DLL4", "MYC", "WPRE-SV40PA"), "Condition")

