library(openxlsx)
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)
library(forcats)

my_palette_Rui_colors_B <- c("#F4A753", "#E95A74", "#50B6EF", "#FC0808",  "#0E47D8", "#45FF8E", "#A80519", "#E28CF4", "#880088", "#C1B80C")

cluster.fgsea <- read.xlsx("./Tables/fgsea_Clusters_Ctrl_Dll4KO_RbpjKO_Notch1KO.xlsx", sheet = 13, sep.names = " ", rowNames = F)

cluster.fgsea[2:11] <- round(cluster.fgsea[2:11], digits = 2)

class(cluster.fgsea)

cluster.fgsea

cluster.fgsea$Color_palette <- my_palette_Rui_colors_B

Hallmarks <- quos(HALLMARK_E2F_TARGETS, HALLMARK_MYC_TARGETS_V1, HALLMARK_MYC_TARGETS_V2, HALLMARK_G2M_CHECKPOINT, HALLMARK_OXIDATIVE_PHOSPHORYLATION)

Hallmarks_quotes <- c("HALLMARK_E2F_TARGETS", "HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2", "HALLMARK_G2M_CHECKPOINT", "HALLMARK_OXIDATIVE_PHOSPHORYLATION")

Hallmarks_quotes_abs <- c("Abs_HALLMARK_E2F_TARGETS", "Abs_HALLMARK_MYC_TARGETS_V1", "Abs_HALLMARK_MYC_TARGETS_V2", "Abs_HALLMARK_G2M_CHECKPOINT", "Abs_HALLMARK_OXIDATIVE_PHOSPHORYLATION")

cluster.fgsea <- cluster.fgsea %>% arrange(desc(!!! Hallmarks[5])) 

cluster.fgsea <- cluster.fgsea %>% arrange(desc(HALLMARK_E2F_TARGETS)) 

cluster.fgsea

cluster.fgsea$To_Analyze = cluster.fgsea[Hallmarks_quotes[5]]

cluster.fgsea$To_Analyze = unlist(cluster.fgsea$To_Analyze)

cluster.fgsea$Abs_values = cluster.fgsea[Hallmarks_quotes_abs[5]]

cluster.fgsea$Abs_values = unlist(cluster.fgsea$Abs_values)

cluster.fgsea$Regulated <- vector("character", nrow(cluster.fgsea))

for(i in 1:nrow(cluster.fgsea)) {cluster.fgsea$Regulated[[i]] <- local({
  i <- i
  { if (cluster.fgsea$To_Analyze[i] > 0) { cluster.fgsea$Regulated[i] = "Upregulated"
  } else if (cluster.fgsea$To_Analyze[i] == 0) { cluster.fgsea$Regulated[i] = "0"
  } else {  cluster.fgsea$To_Analyze[i] = "Downregulated"
  }}
})}

class(cluster.fgsea$Regulated)

cluster.fgsea$X1 <- as.factor(cluster.fgsea$X1)

cluster.fgsea$Regulated <- as.factor(cluster.fgsea$Regulated)

cluster.fgsea$ymax <- vector("numeric", nrow(cluster.fgsea))

for(i in 1:nrow(cluster.fgsea)) {cluster.fgsea$ymax[[i]] <- local({
  i <- i
  cluster.fgsea$ymax[i] = cumsum(cluster.fgsea$Abs_values[1:i])[i]
})}

cluster.fgsea$ymin <- vector("numeric", nrow(cluster.fgsea))

cluster.fgsea$ymin[2:nrow(cluster.fgsea)] <- cluster.fgsea$ymax[1:rows]

rows <- nrow(cluster.fgsea)-1

# Go from here

donut_text = cluster.fgsea$Abs_values/2 + c(0,cumsum(cluster.fgsea$Abs_values)[-length(cluster.fgsea$Abs_values)])

donutpie <- ggplot(cluster.fgsea) + 
  geom_rect(aes(fill = X1, ymax=ymax, ymin=ymin, xmax=4, xmin=2), fill=cluster.fgsea$Color_palette) +
  geom_label(data = cluster.fgsea, aes(x = 3, y = donut_text, label = To_Analyze),alpha=0.5, label.r = unit(0.3, "lines"), label.padding = unit(0.15, "lines"), size = 4.8)+
  geom_rect(aes(fill=Regulated, ymax=ymax, ymin=ymin, xmax=2, xmin=0)) + scale_fill_manual(values = regulated)+
  # geom_text(aes(x = 3, y = donut_text, label = HALLMARK_E2F_TARGETS, group = X1, vjust = 0), vjust = 0, size = 5)+
  xlim(c(0, 4)) + 
  theme(aspect.ratio=1) +
  coord_polar(theta="y")+
  theme_minimal()


donutpie

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )

p1 <- donutpie+
  blank_theme +
  theme_void()+NoLegend()+ggtitle(gene.sets.names[5])+ theme(text = element_text(size = 20))

p1

cairo_pdf("Plots/Figure_6/Donut_Pie_Charts/Figure_6C_Donut_Pie_Chart_OXIDATIVE_PHOSPHORYLATION.pdf", width = 5, height = 5, family = "Arial")

p1

dev.off() 



