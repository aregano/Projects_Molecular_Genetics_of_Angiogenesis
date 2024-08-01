library(openxlsx)
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)
library(forcats)

my_palette_Rui_colors_B <- c("#F4A753", "#E95A74", "#50B6EF", "#FC0808",  "#0E47D8", "#45FF8E", "#A80519", "#E28CF4", "#880088", "#C1B80C")

cluster.fgsea <- read.xlsx("./Tables/fgsea_Clusters_Ctrl_Dll4KO_RbpjKO_Notch1KO.xlsx", sheet = 13, sep.names = " ", rowNames = F)

cluster.fgsea[2:11] <- round(cluster.fgsea[2:11], digits = 2)


#  Barplot of Gene Sets with positive and negative values stacked


cluster.fgsea.barchart <- expand_grid(
  HALMARK    = Hallmarks_quotes,  # Define all unique student names
  color = my_palette_Rui_colors_B,     # Define all unique HW assignments
  value      = NA)
cluster.fgsea.barchart

cluster.fgsea.barchart.long <- pivot_longer(cluster.fgsea, cols=c(2,3,5,6), names_to = "Hallmark", values_to = "Cluster")

cluster.fgsea.barchart <- cluster.fgsea.barchart.long[c(1, 8:10)] %>% arrange(Hallmark)

write.xlsx(cluster.fgsea.barchart, "./Tables/fgsea_Clusters_Ctrl_Dll4KO_RbpjKO_Notch1KO.Barchart.xlsx")


dat1 <- subset(cluster.fgsea.barchart,Cluster >= 0)
dat2 <- subset(cluster.fgsea.barchart,Cluster < 0)


p1 <- ggplot() + 
  geom_bar(data = dat1, aes(x=Hallmark, y=Cluster, fill=X1),stat = "identity") +
  geom_bar(data = dat2, aes(x=Hallmark, y=Cluster, fill=X1),stat = "identity") + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 12), axis.title = element_blank(), legend.title = element_blank())+
  scale_fill_manual(values= cluster.fgsea.barchart$Color_palette)+
  scale_y_continuous(breaks=seq(-90,70,10))

p1

cairo_pdf("Plots/Figure_6/Donut_Pie_Charts/Figure_6C_Barchart_HallMark_Gene_Sets.pdf", width = 5, height = 10, family = "Arial")
p1
dev.off()

p1 <- ggplot(cluster.fgsea.barchart, aes(x = Hallmark, y = Cluster, 
                                         fill = X1, label = Cluster)) +
  geom_bar(stat = "identity") + geom_text(
    size = 3, position = position_stack(vjust = 0.5),colour = "white")   +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 12), axis.title = element_blank(), legend.title = element_blank())+
  scale_fill_manual(values= cluster.fgsea.barchart$Color_palette)+
  scale_y_continuous(breaks=seq(-90,70,10))

cairo_pdf("Plots/Figure_6/Donut_Pie_Charts/Figure_6C_Barchart_HallMark_Gene_Sets_with_numbers.pdf", width = 5, height = 10, family = "Arial")
p1
dev.off()

#####################################

# Trying out different ways to get the Barchar in order of values (NOT SOLVED)

ggplot(cluster.fgsea.barchart, aes(x = Hallmark, y = reorder(Cluster,-Cluster), 
                                   fill = X1, label = Cluster)) +
  geom_bar(stat = "identity") + geom_text(
    size = 3, position = position_stack(vjust = 0.5),colour = "white")   +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 12), axis.title = element_blank(), legend.title = element_blank())+
  scale_fill_manual(values= cluster.fgsea.barchart$Color_palette)+
  scale_y_continuous(breaks=seq(-90,70,10))


ggplot(cluster.fgsea.barchart, aes(Hallmark, Cluster, fill = interaction(-Cluster, X1))) + 
  geom_col() +
  scale_fill_manual(values = cluster.fgsea.barchart$Color_palette,
                    labels = with(my.df, X1[aux]), 
                    breaks = with(my.df, interaction(-Cluster, Hallmark)[aux]))


my.df <- cluster.fgsea.barchart %>% 
  arrange(Hallmark, X1) %>% 
  mutate(type = factor(X1)) %>% 
  arrange(Hallmark, -Cluster) 

aux <- with(my.df, match(sort(unique(X1)), X1))

ggplot(my.df, aes(Hallmark, Cluster, fill = interaction(-Cluster, Hallmark))) + 
  geom_col() + 
  scale_fill_manual(values = cluster.fgsea.barchart$Color_palette,
                    labels = with(my.df, X1[aux]), 
                    breaks = with(my.df, interaction(-Cluster, Hallmark)[aux]))


ggplot(data = cluster.fgsea.barchart, 
       aes(x = Hallmark, 
           y = Cluster, 
           fill = fct_reorder(X1, Cluster),
           label = Cluster)) + 
  geom_bar(stat = "identity") +
  guides(fill = guide_legend(title = "substrate"))+
  geom_bar(stat = "identity") + geom_text(
    size = 3, position = position_stack(vjust = 0.5),colour = "white")   +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 12), axis.title = element_blank(), legend.title = element_blank())+
  scale_fill_manual(values= cluster.fgsea.barchart$Color_palette)+
  scale_y_continuous(breaks=seq(-90,70,10))