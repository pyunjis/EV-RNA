# Upset plot of how many DE genes are shared for Supp Fig 10
#library(VennDiagram) 
library(reshape2)
library(dplyr)
library(UpSetR)
library(ggplot2)


# set working directory
setwd("E:/4. EV-RNA/4. EV-RNA/4. Code/2023-07-23 FINAL/")

# set_LG <- c("FR3033", "FR2326", "FR1619","FR912", "FR58", "FR14")
# set_LV <- c("FR912", "FR58", "FR14")
# set_MM <- c("FR2326", "FR1619","FR912", "FR58", "FR14")




Upset_plot <- function(target, WIDTH){

  df1 <- read.table(paste0("./Files/DESeq2_fraction/unnorm_HD_FR14_vs_", target, "_FR14_deseq2_results.txt", sep=""), header = TRUE, stringsAsFactors = FALSE, sep = "\t",
                    check.names = FALSE)
  df2 <- read.table(paste0("./Files/DESeq2_fraction/unnorm_HD_FR58_vs_", target, "_FR58_deseq2_results.txt", sep=""), header = TRUE, stringsAsFactors = FALSE, sep = "\t",
                    check.names = FALSE)
  df3 <- read.table(paste0("./Files/DESeq2_fraction/unnorm_HD_FR912_vs_", target, "_FR912_deseq2_results.txt", sep=""), header = TRUE, stringsAsFactors = FALSE, sep = "\t",
                    check.names = FALSE)
  df4 <- read.table(paste0("./Files/DESeq2_fraction/unnorm_HD_FR1619_vs_", target, "_FR1619_deseq2_results.txt", sep=""), header = TRUE, stringsAsFactors = FALSE, sep = "\t",
                    check.names = FALSE)
  df5 <- read.table(paste0("./Files/DESeq2_fraction/unnorm_HD_FR2326_vs_", target, "_FR2326_deseq2_results.txt", sep=""), header = TRUE, stringsAsFactors = FALSE, sep = "\t",
                    check.names = FALSE)
  df6 <- read.table(paste0("./Files/DESeq2_fraction/unnorm_HD_FR3033_vs_", target, "_FR3033_deseq2_results.txt", sep=""), header = TRUE, stringsAsFactors = FALSE, sep = "\t",
                    check.names = FALSE)
  
  # Filter for sig.genes
  df1 <- df1[df1$padj <0.05 & df1$log2FoldChange>1 ,]
  df2 <- df2[df2$padj <0.05 & df2$log2FoldChange>1 ,]
  df3 <- df3[df3$padj <0.05 & df3$log2FoldChange>1 ,]
  df4 <- df4[df4$padj <0.05 & df4$log2FoldChange>1 ,]
  df5 <- df5[df5$padj <0.05 & df5$log2FoldChange>1 ,]
  df6 <- df6[df6$padj <0.05 & df6$log2FoldChange>1 ,]
  
  
  
  df1$ensID<- rownames(df1)
  df2$ensID<- rownames(df2)
  df3$ensID<- rownames(df3)
  df4$ensID<- rownames(df4)
  df5$ensID<- rownames(df5)
  df6$ensID<- rownames(df6)
  
  
  df1 <- df1$ensID
  df2 <- df2$ensID
  df3 <- df3$ensID
  df4 <- df4$ensID
  df5 <- df5$ensID
  df6 <- df6$ensID
  
  cluster2_list <- list(df1,df2,df3,df4,df5,df6)
  gene_list <- list()
  types <- c(1,2,3,4,5,6)
  names <- c("FR14","FR58","FR912","FR1619","FR2326", "FR3033")
  for (i in 1:(length(types))) {
    baseline <- types[i]
    #cname <- paste0(types[1])
    cname <- names[i]
    gene_list[[cname]] <- cluster2_list[[baseline]]
  }
  
  # confirm the number and elements
  #a <- intersect(df1, df2)
  
  # upset(fromList(gene_list),  sets = set_type, keep.order = TRUE,
  #       order.by = c("degree", "freq"), decreasing = c(TRUE, TRUE),  text.scale = 2)
  # To save pdf

  #pdf(paste0("./Figures/UpsetPlot_", target,".pdf",sep=""), 8,5)
  pdf(paste0("./Figures/Supp_Figure_10_UpsetPlot_", target,".pdf",sep=""), paper = "a4r", width = WIDTH, height = 5.5)
  #upset(fromList(gene_list), order.by = c("freq","degree"), decreasing = c(TRUE, TRUE), text.scale = 2)
  a <- upset(fromList(gene_list),  sets = c("FR3033", "FR2326", "FR1619","FR912", "FR58", "FR14"), keep.order = TRUE,
        order.by = c("degree", "freq"), decreasing = c(TRUE, TRUE),  text.scale = 2)
  print(a)
  dev.off()
  
  # # To see length of each vector for the groups
  # str(expressed)
}


Upset_plot(target = "MM", WIDTH = 8) # pdf 8, 5.5
Upset_plot(target = "LV", WIDTH = 4) # pdf 4, 5.5
Upset_plot(target = "LG", WIDTH = 8) # pdf 8, 5.5
       