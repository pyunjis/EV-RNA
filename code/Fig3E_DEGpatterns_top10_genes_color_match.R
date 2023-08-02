# Heatmap of top 10 genes from each cluster for figure 3E

library(RColorBrewer)
library(edgeR)
library(RUVSeq)
library(reshape2)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(tidyverse)
library(gplots)

# Set working directory
setwd("E:/4. EV-RNA/4. EV-RNA/4. Code/2023-07-23 FINAL/")

# Import Data
ALL_cluster <- read.table("./Files/degpatterns/cluster_groups_FDR0.05_degpatterns_ALL_FR_minc0.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t",
                         check.names = FALSE)

ALL_norm <- read.table("./Files/degpatterns/norm_FDR0.05_degpatterns_ALL_FR_minc0.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t",
                   check.names = FALSE)
cts <- read.table("./Files/set1df_unnorm_protein_coding_deseq_norm.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t",
                  check.names = FALSE)
md <- read.table("./Files/metadata_20200421.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t",
                 check.names = FALSE)

# plot heatmap of top 10 genes
RunHeatmap_top10 <- function(cluster, type, norm){
  
  # Select top 10 genes from each cluster
  top_cluster <- cluster %>%
    group_by(cluster) %>%
    top_n(n=-10, wt = FDR)

  # Look at for specific condition
  # remove low quality samples
  md <- filter(md, !SampleID %in% c("EV065", "EV066"))
  
  # Remove precancer types
  md <- filter(md, Condition == "HD" | Condition == "LG" | Condition == "LV" | Condition == "MM")
  rownames(md) <- md$SampleID

  # Filter
  keep <- rownames(cts)[rownames(cts) %in% top_cluster$genes]
  cts2 <- cts[keep,rownames(md)]
  cts3 <- log(cts2+1, 2)

  md2 <- md[,c(1,4)]

  # re-ordering column
  md2$Fraction <- factor(md2$Fraction, levels = c("FR14", "FR58", "FR912", "FR1619", "FR2326", "FR3033"))
  md3 <- md2 %>%
    arrange(Fraction) 
  #FR_levels <- md3$Fraction
  
  SampleID_levels <- md3$SampleID

  cts4 <- cts3[,SampleID_levels]
  rownames(md3) <- md3$SampleID
  
  # re-ordering rows
  gene_order <-  top_cluster %>%
    arrange(cluster)
  
  gene_order_level<- gene_order$genes
  gene_order$FDR <- NULL
  rownames(gene_order) <- gene_order$genes
  gene_order2 <- gene_order
  gene_order$genes<- NULL
  rownames(gene_order) <- gene_order2$genes
  
  #Cluster <- RColorBrewer::brewer.pal(5, "Set1")[c(1,2,3)]
  library(scales)
  #show_col(hue_pal()(12))
  cluster <- hue_pal()(12)
  #cluster <- brewer.pal(6, "Set2")
  names(cluster) <- unique(as.character(gene_order$cluster))
  
  ## Generate palette and map this to each Condition, and to each Fraction
  fractionCols <- brewer.pal(length(unique(md3$Fraction)), "Set1")
  names(fractionCols) <- unique(md3$Fraction)
  
  annotColors <- list(Fraction = fractionCols, cluster = cluster)
  breaksList = seq(-2.5, 2.5, by = 0.1)
  md3$SampleID <- NULL
  
  cts5 <- cts4
  cts5 <- cts5[rownames(gene_order),]
  
  # gene_order has to be a datframe
  df <- as.data.frame(gene_order)
  df$cluster<- as.character(df$cluster)
  #pheatmap 
  breaks <- seq(-2.5, 2.5, by =0.1)

  clusterHeat1 <- pheatmap(cts5, show_rownames=F, show_colnames=F,
                           breaks = breaks,
                           annotation_names_row = TRUE,
                           annotation_names_col = FALSE,
                           annotation_colors = annotColors,
                           cluster_rows = FALSE,
                           cluster_cols = FALSE,
                           color = colorRampPalette(c("white", "white", "black"))(50),
                           annotation_col = md3,
                           scale = "row",
                           fontsize_row=4, fontsize_col=6, fontsize=8,
                           border_color=NA,
                           annotation_row = df
  )
  pdf(paste0("./Figures/Figure_3E_DEGpattern_top10_heatmap_ALL_color_no_cluster_minc0.pdf"), 6, 4)
  print(clusterHeat1)
  dev.off()
  
  clusterHeat2 <- pheatmap(cts5, show_rownames=F, 
                           breaks = breaks, show_colnames=F, 
                           annotation_names_row = TRUE,
                           annotation_names_col = FALSE,
                           annotation_colors = annotColors,
                           #cluster_rows = FALSE,
                           cluster_cols = FALSE,
                           color = colorRampPalette(c("white", "white", "black"))(50),
                           annotation_col = md3,
                           scale = "row",
                           fontsize_row=4, fontsize_col=6, fontsize=8,
                           border_color=NA,
                           annotation_row = df
  )
  
  
  #pdf(paste0("./Figures/DEGpattern_top10_heatmap_ALL_color_cluster_minc0.pdf"), 6, 4)
  #print(clusterHeat2)
  #dev.off()
  

  
}

RunHeatmap_top10(cluster = ALL_cluster, norm = ALL_norm)
# RunHeatmap_top10(cluster = MM_cluster, type = "MM", norm = MM_norm)
# RunHeatmap_top10(cluster = LG_cluster, type = "LG", norm = LG_norm)
# RunHeatmap_top10(cluster = LCir_cluster, type = "LCir", norm = LCir_norm)
# RunHeatmap_top10(cluster = MG_cluster, type = "MG", norm = MG_norm)
