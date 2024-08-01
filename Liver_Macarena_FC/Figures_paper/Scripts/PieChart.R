# PieChart
library(ggfittext)
library(ggplot2)
library(ggrepel)
library(openxlsx)
library(lessR)
library(Seurat)
library(reshape2)
library(dplyr)
library(tibble)
library(tidyr)
library(webr)

my_palette_Rui_colors_B <- c("#F4A753", "#E95A74", "#50B6EF", "#FC0808",  "#0E47D8", "#45FF8E", "#A80519", "#E28CF4", "#880088", "#C1B80C")

regulated <- c( "royalblue", "chartreuse")


cluster.fgsea <- read.xlsx("./Tables/GSEA_Fig4/fgsea_Dll4KO_Clusters.xlsx ", sheet = 11, sep.names = " ", rowNames = T)

cluster.fgsea <- round(cluster.fgsea, digits = 2)

gene.sets <- c(Abs_HALLMARK_E2F_TARGETS, Abs_HALLMARK_MYC_TARGETS_V1, Abs_HALLMARK_MYC_TARGETS_V2, Abs_HALLMARK_G2M_CHECKPOINT, Abs_HALLMARK_OXIDATIVE_PHOSPHORYLATION)

gene.sets.names <- c("E2F TARGETS", "MYC TARGETS V1", "MYC TARGETS V2", "G2M CHECKPOINT","OXIDATIVE PHOSPHORYLATION")


class(cluster.fgsea)



slices <- cluster.fgsea$Abs_HALLMARK_E2F_TARGETS
lbls <- rownames(cluster.fgsea)
pie(slices, labels = lbls, main="HallMark E2F Targets", col = my_palette_Rui_colors_B)
par(xpd=TRUE)
text(1.1,0,"A,B,C,D=0")
par(xpd=FALSE)

Abs_HALLMARK_E2F_TARGETS

bp<- ggplot(cluster.fgsea, aes(x="", y=Abs_HALLMARK_OXIDATIVE_PHOSPHORYLATION, fill=lbls))+
  geom_bar(width = 1, stat = "identity")

bp

pie <- bp + coord_polar("y", start=0)
pie

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )

p1 <- pie + scale_fill_manual(values=my_palette_Rui_colors_B)+ 
  blank_theme +
    geom_text(aes(label = HALLMARK_OXIDATIVE_PHOSPHORYLATION),
            position = position_stack(vjust = 0.5), show.legend = F)+
  theme_void()+NoLegend()+ggtitle(gene.sets.names[5])

  
cairo_pdf("Plots/Figure_6/Pie_Charts/Figure_6C_Pie_Chart_OXIDATIVE_PHOSPHORYLATION.pdf", width = 5, height = 5, family = "Arial")

p1

dev.off()  


###################### MULITPLE PIE CHARTS ###########################

cluster.fgsea <- read.xlsx("./Tables/GSEA_Fig4/fgsea_Dll4KO_Clusters.xlsx ", sheet = 11, sep.names = " ", rowNames = F)
  
  
cluster.fgsea.m = melt(cluster.fgsea)

cluster.fgsea.d <- cluster.fgsea.m[51:100,]
  

bp<- ggplot(cluster.fgsea, aes(x="", y=Abs_HALLMARK_E2F_TARGETS, fill=lbls))+
  geom_bar(width = 1, stat = "identity")
  
ggplot(data=cluster.fgsea.d, aes(x=" ", y=value, group=Cluster, fill=Cluster)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start=0) + 
    facet_grid(.~ variable) +
    theme_void()+
    scale_fill_manual(values=my_palette_Rui_colors_B)

##############################PieChart()#######################################


# Creamos el data frame con la nueva variable
mis_datos <- data.frame(x = c(rep(count_2[1], count_2[1]),  # 15 veces 15
                              rep(count_2[2], count_2[2]),  # 27 veces 27
                              rep(count_2[3], count_2[3]),  # 25 veces 25
                              rep(count_2[4], count_2[4]))) # 10 veces 10

e2f

genero <- factor(c(rep("Hombre", 10), rep("Mujer", 20)))
genero

tabla_genero <- table(genero)
tabla_genero


genero <- data.frame(gen = genero)

PieChart(gen,
                hole = 0,
                values = "%", 
                data = genero,
         fill = c("lightblue", "pink"),
         main = "", theme = lessR::style())

e2f <- data.frame(slices = e2f)

e2f <- as.factor(cluster.fgsea$HALLMARK_E2F_TARGETS)

my_data <- data.frame(cluster.fgsea$HALLMARK_E2F_TARGETS)

PieChart(HALLMARK_E2F_TARGETS, y = HALLMARK_E2F_TARGETS, hole = 0,
                data = cluster.fgsea, values = "input", values_digits = 0,
                fill = my_palette_Rui_colors_B, main = "")

count_2 <- c(15, 27, 25, 10)

my_data <- data.frame(x = c(rep(count_2[1], count_2[1]),  # 15 times 15
                            rep(count_2[2], count_2[2]),  # 27 times 27
                            rep(count_2[3], count_2[3]),  # 25 times 25
                            rep(count_2[4], count_2[4]))) # 10 times 10

PieChart(x, hole = 0, values = "%", data = my_data, fill = 1:4, main = "")

?lessR::PieChart


###################################Doughnut Chart#############################################

# Create test data.
data <- data.frame(
  category=c("A", "B", "C"),
  count=c(10, 60, 30)
)

# Compute percentages
data$fraction = data$count / sum(data$count)

# Compute the cumulative percentages (top of each rectangle)
data$ymax = cumsum(data$fraction)

# Compute the bottom of each rectangle
data$ymin = c(0, head(data$ymax, n=-1))

# Make the plot
ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
  geom_rect() +
  coord_polar(theta="y")+  # Try to remove that to understand how the chart is built initially
  xlim(c(2, 4)) # Try to remove that to see how to make a pie chart



#########################################################################################

# IMPORTANT! 
# In order to go around the different Gene Sets, you just need to change the index order of: Hallmarks, Hallmarks_quotes, Hallmarks_quotes_abs and gene.sets.names
# Finally, remember to store each pdf file under a different name

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

############################################ TUTORIALS ####################################

metadata <- data.frame(tibble::tribble(
  ~Domain,  ~Total,                    ~Group,    ~Share,
  "Eukaryotes PR2", 2138457,               "Algae PR2", "1320398",
  "Bacteria Silva 16S", 2594766, "Cyanobacteria Silva 16S", "1524609",
  "Eukaryotes Silva 18S", 1237100,         "Algae Silva 18S",  "796201",
  "Archaea NCBI",   8e+05,                 "Nogroup", "Nogroup"
)
)

#Make Share numeric
metadata <- mutate(metadata, Share = ifelse(Share == "Nogroup", "0", Share)) %>% 
  mutate(Share = as.numeric(Share))

#Calculate the size and invent Group labels for the unnamed fraction of each Domain
metadata2 <- data.frame(Domain = metadata$Domain,
                        Total = metadata$Total,
                        Group = paste("Not", metadata$Group),
                        Share = metadata$Total - metadata$Share,
                        stringsAsFactors = FALSE)

#Make a data set where each Group is explicit
metadata2 <- rbind(metadata2, metadata) %>% 
  filter(Share != 0) %>% 
  mutate(Group = ifelse(Group == "Not Nogroup", Domain, Group)) %>%
  arrange(Domain, Group) 

metadata2 <- metadata2 %>%
  mutate(Tot = sum(Share)) %>% 
  group_by(Domain) %>% 
  mutate(CUM = cumsum(Share), DomSize = max(CUM))

#Calculate the bottom edge of the Domains when stacked
DomBot <- unique(select(metadata2, Domain, Tot, DomSize)) %>% ungroup() %>% 
  mutate(Bottom = Tot - cumsum(DomSize))

metadata2 <- inner_join(metadata2, select(DomBot, Domain, Bottom))
#> Joining, by = "Domain"
metadata2 <- mutate(metadata2, Pos = Bottom + CUM - Share/2)


plt <- ggplot() + geom_col(aes(x = 2, y = Total, fill = Domain), 
                           data = metadata, color = "black") + 
  geom_col(aes(x = 3, y = Share, fill = Domain), 
           data = metadata2, color = "black") +
  geom_text(aes(label = Group, x= 3, y = Pos), data = metadata2, size = 3)+
  xlim(0, 3.5) + labs(x = NULL, y = NULL) + 
  theme(axis.ticks=element_blank(),
        axis.text=element_blank(),
        axis.title=element_blank())

plt

plt + coord_polar(theta = "y") 



