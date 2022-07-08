# Figure 7A Selective packaging heatmap

# Load library
library(DESeq2)
library(dplyr)
library(pheatmap)
library(ggplot2)
library(gplots)
library(reshape2)
library(RColorBrewer)

# Set working directory
setwd("D:/Box drive update/4. EV-RNA/4. Code/2022-02-17 V7/")

# Import Data
# cts <- read.table("./df/2022-01-27/set1df_unnorm_protein_coding_deseq_norm.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t",
#                   check.names = FALSE) # Was normalized using 11609:11664

cts <- read.table("./Files/set1df_unnorm_protein_coding_deseq_norm.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t",
                  check.names = FALSE)

md <- read.table("./Files/metadata_20200421.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t",
                 check.names = FALSE)

# Grab row names of cts for ERCC
spikes <- rownames(cts)[grep("^ERCC", rownames(cts))]
genes <- rownames(cts)[grep("^ENS", rownames(cts))]

# how many genes and ERCC spikes are we including
length(genes) # 11609
length(spikes) #55

# Filter to include only genes
keep <- rownames(cts)[rownames(cts) %in% genes]
cts <- cts[keep,] 


cluster_heatmap <- function(Target, Baseline){
  # Load DESeq2 results for each fraction
  FR14 <- read.table(paste0("./Files/DESeq2_fraction/unnorm_",Baseline,"_FR14_vs_",Target,"_FR14_deseq2_results.txt", sep = ""), header = TRUE, stringsAsFactors = FALSE, sep = "\t",
                     check.names = FALSE)
  FR58 <- read.table(paste0("./Files/DESeq2_fraction/unnorm_",Baseline,"_FR58_vs_",Target,"_FR58_deseq2_results.txt", sep = ""), header = TRUE, stringsAsFactors = FALSE, sep = "\t",
                     check.names = FALSE)
  FR912 <- read.table(paste0("./Files/DESeq2_fraction/unnorm_",Baseline,"_FR912_vs_",Target,"_FR912_deseq2_results.txt", sep = ""), header = TRUE, stringsAsFactors = FALSE, sep = "\t",
                      check.names = FALSE)
  FR1619 <- read.table(paste0("./Files/DESeq2_fraction/unnorm_",Baseline,"_FR1619_vs_",Target,"_FR1619_deseq2_results.txt", sep = ""), header = TRUE, stringsAsFactors = FALSE, sep = "\t",
                       check.names = FALSE)
  FR2326 <- read.table(paste0("./Files/DESeq2_fraction/unnorm_",Baseline,"_FR2326_vs_",Target,"_FR2326_deseq2_results.txt", sep = ""), header = TRUE, stringsAsFactors = FALSE, sep = "\t",
                       check.names = FALSE)
  FR3033 <- read.table(paste0("./Files/DESeq2_fraction/unnorm_",Baseline,"_FR3033_vs_",Target,"_FR3033_deseq2_results.txt", sep = ""), header = TRUE, stringsAsFactors = FALSE, sep = "\t",
                       check.names = FALSE)
  dfList <- list(FR14, FR58, FR912, FR1619, FR2326, FR3033)
  names <- unique(md$Fraction)
  pos_list <- list()
  neg_list <- list()
  
  
  # Filter important gene based on log2FC & padj < 0.05 
  for (i in 1:6) {
    type <- dfList[[i]]  
    name <- names[[i]]
    # Grab log2FC>1 & log2FC<1
    pos <- type[type$log2FoldChange > 1 & type$padj < 0.05,]
    neg <- type[type$log2FoldChange < 1 & type$padj < 0.05,]
    pos$gene <- rownames(pos) 
    neg$gene <- rownames(neg)
    # Add cluster information by adding FR type into a column
    type$cluster <- paste(name)
    pos_list[[name]] <- pos
    neg_list[[name]] <- neg
    
  }  
  
  
  #Combine data frame into a list and extrapolate "union" by rowname
  idx <- Reduce(union, lapply(pos_list, rownames))
  cts_pos <- cts[idx, ]
  
  idx <- Reduce(union, lapply(neg_list, rownames))
  cts_neg <- cts[idx, ]
  
  # Filter md
  md <- dplyr::filter(md, Condition == "HD" | Condition == Target)
  rownames(md) <- md$SampleID
  md2 <- md
  
  # For LV065 & LV066 need to be removed
  md2 <- filter(md2, !SampleID %in% c("EV065", "EV066"))
  rownames(md2) <- md2$SampleID
  
  # keep only relevant SampleID in your count table
  filtered <- cts_pos
  keep <- colnames(filtered)[colnames(filtered) %in% rownames(md2)]
  filtered_hd <- filtered[,keep]  # no log2 counts
  
  # grab only gene names
  genes <- rownames(filtered_hd)[grep("^ENSG", rownames(filtered_hd))]
  #length(genes)
  cts_all <- filtered_hd[genes,] # counts with all fraction
  
  type <- "HD"
  HD_levels<- paste(type, c("FR14", "FR58","FR912", "FR1619", "FR2326", "FR3033"), sep = "_")
  type2 <- Target
  CX_levels<- paste(type2, c("FR14", "FR58","FR912", "FR1619", "FR2326", "FR3033"), sep = "_")
  
  # Keep the levels for md2
  md2$Condition_FR <- factor(md2$Condition_FR, levels=c(HD_levels, CX_levels))
  md2$Fraction <- factor(md2$Fraction, levels=c("FR14", "FR58", "FR912", "FR1619", "FR2326", "FR3033"))
  
  # Tidy the data
  cts_all$gene <- rownames(cts_all)
  ctsMelt <- melt(cts_all, by = "gene") %>%
    mutate(log2counts = log2(value + 1)) 
  
  # Rename columns
  colnames(ctsMelt)[2:4] <- c("SampleID","counts", "log2counts")
  
  # Add type info
  iv <- match(ctsMelt$SampleID, md2$SampleID)
  ctsMelt$Fraction <- md2[iv,]$Fraction
  ctsMelt$Condition <- md2[iv,]$Condition
  ctsMelt$Condition_FR <- md2[iv,]$Condition_FR
  
  fracs <- c("FR14","FR58","FR912","FR1619","FR2326","FR3033")
  comp_list <-  list()
  cts_MM <- ctsMelt
  
  # Calculate average healthy counts across each fraction and standardize type counts by this
  for (i in 1:length(fracs)) {
    f <- fracs[i]
    
    df_HD <- cts_MM %>%
      filter(Fraction == f, Condition == "HD") %>%
      group_by(gene, Condition_FR) %>%
      mutate(HD_mean = mean(log2counts)) %>%
      mutate(delta = log2counts-HD_mean)
    
    df_type <- cts_MM %>%
      filter(Fraction == f, Condition == Target)
    df_type$HD_mean <- df_HD$HD_mean
    df_type$delta <- df_type$log2counts - df_HD$HD_mean
    
    # Define variable names for each dataframe
    df_name_HD <- paste("cts","HD",f,sep = "_")
    df_name_type <- paste("cts","CX",f,sep = "_")
    
    comp_list[[df_name_HD]] <- df_HD
    comp_list[[df_name_type]] <- df_type
  }
  
  #combine all sub-dataframes with delta
  delta <- bind_rows(comp_list$cts_HD_FR14, comp_list$cts_HD_FR58, comp_list$cts_HD_FR912, comp_list$cts_HD_FR1619, comp_list$cts_HD_FR2326, comp_list$cts_HD_FR3033, 
                     comp_list$cts_CX_FR14, comp_list$cts_CX_FR58, comp_list$cts_CX_FR912, comp_list$cts_CX_FR1619, comp_list$cts_CX_FR2326, comp_list$cts_CX_FR3033)
  delta2 <- delta[,c(1,2,9)] # grab gene, sampleID, and delta
  
  # dcast the dataframe
  datac <- dcast(delta2, gene ~ SampleID, value.var="delta")
  rownames(datac) <- datac$gene
  datac$gene <- NULL
  
  # Re-order the SampleID for biological replicates to be together within fraction
  md_3 <- md2
  #md_3 <- md_3 %>% filter(Condition %in% c("HD", type))
  md_3$Condition_FR <- factor(md_3$Condition_FR, levels = c(HD_levels, CX_levels))
  md_3 <- md_3 %>% arrange(Condition_FR)
  rownames(md_3) <- md_3$SampleID
  md_3$SampleID <- NULL
  ordered_list <- rownames(md_3)
  
  plot1 <- datac[,ordered_list]
  
  md_4 <- md_3[,c(1,3)]
  
  ## Generate proper color scheme for heatmap
  fractionCols <- brewer.pal(length(unique(md_4$Fraction)), "Set1")
  names(fractionCols) <- unique(md_4$Fraction)
  
  cols <- read.delim("./Files/colors.txt", stringsAsFactors = F)
  
  conditionCols <- cols[cols$Group %in% c("HD",Target),]$Colour
  names(conditionCols) <- cols[cols$Group %in% c("HD",Target),]$Group
  
  ### annotation of gene clusters for pheatmap
  cluster <- brewer.pal(6, "Set1")
  names(cluster) <- sub("^", "Cluster_", seq(1,6))
  
  # Specify colour scheme
  #ann_colors = list(Cluster = cluster, Fraction = fractionCols, Condition = conditionCols)
  ann_colors = list(Fraction = fractionCols, Condition = conditionCols)
 
  # Heatmap without any cluster
  breaks <- seq(-3.5, 3.5, by =0.1)
  clusterHeat <- pheatmap(plot1, show_rownames=F, 
                          breaks = breaks,
                          #annotation_row = my_gene_col,
                          annotation_names_row = TRUE,
                          annotation_names_col = FALSE,
                          #annotation_colors = ann_colors,
                          #cluster_rows = FALSE, 
                          cluster_cols = FALSE,
                          redgreen(70),
                          annotation_col = md_4, 
                          scale = "row", 
                          fontsize_row=4, fontsize_col=6, fontsize=8,
                          border_color=NA
  )
  save_pheatmap_pdf <- function(x, filename) {
    stopifnot(!missing(x))
    stopifnot(!missing(filename))
    pdf(filename, height = 4)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
  }
  
  #save_pheatmap_pdf(clusterHeat, paste0("./Heatmap/2021-10-19/Figure_4C_ournorm_protein_coding_", Target, "_", Baseline, ".pdf"))
  
  # New - Assign cluster by order of DE identified
  md_4 <- md_3[,c(1,3)]
  
  ## Generate proper color scheme for heatmap
  fractionCols <- brewer.pal(length(unique(md_4$Fraction)), "Set1")
  names(fractionCols) <- unique(md_4$Fraction)
  
  cols <- read.delim("./Files/colors.txt", stringsAsFactors = F)
  
  conditionCols <- cols[cols$Group %in% c("HD",Target),]$Colour
  names(conditionCols) <- cols[cols$Group %in% c("HD",Target),]$Group
  
  ### annotation of gene clusters for pheatmap
  # cluster <- brewer.pal(6, "Set1")
  # names(cluster) <- sub("^", "Cluster_", seq(1,6))
  
  # New * add annotation for fraction DE genes identified
  FR14_DE <- brewer.pal(6, "Set1")[1]
  names(FR14_DE) <- list("FR14")
  FR58_DE <- brewer.pal(6, "Set1")[2]
  names(FR58_DE) <- list("FR58")
  FR912_DE <- brewer.pal(6, "Set1")[3]
  names(FR912_DE) <- list("FR912")
  FR1619_DE <- brewer.pal(6, "Set1")[4]
  names(FR1619_DE) <- list("FR1619")
  FR2326_DE <- brewer.pal(6, "Set1")[5]
  names(FR2326_DE) <- list("FR2326")
  # FR3033_DE <- brewer.pal(6, "Set1")[6]
  # names(FR3033_DE) <- list("FR3033")
  
  # Specify colour scheme
  #ann_colors = list(Cluster = cluster, Fraction = fractionCols, Condition = conditionCols)
  ann_colors = list(FR14_DE = FR14_DE, FR58_DE = FR58_DE, FR912_DE = FR912_DE, FR1619_DE = FR1619_DE, FR2326_DE = FR2326_DE, Fraction = fractionCols, Condition = conditionCols)
  
  # Let's add column of fraction DE genes identified
  # need to make a data frame with a gene in one column
  gene_df <- plot1
  gene_df$gene <- rownames(gene_df)
  geneDF <- as.data.frame(gene_df[,61])
  geneDF$gene <- geneDF$`gene_df[, 61]`
  rownames(geneDF) <- geneDF$gene
  geneDF$`gene_df[, 61]` <- NULL
  test <- geneDF 
  
  pos_list$FR14$Fraction <- "FR14"
  pos_list$FR58$Fraction <- "FR58"
  pos_list$FR912$Fraction <- "FR912"
  pos_list$FR1619$Fraction <- "FR1619"
  pos_list$FR2326$Fraction <- "FR2326"
  # pos_list$FR3033$Fraction <- "FR3033"
  
  # Add FR DE identified
  iv <- match(test$gene, pos_list$FR14$gene)
  test$FR14_DE <- pos_list$FR14[iv,]$Fraction
  
  iv <- match(test$gene, pos_list$FR58$gene)
  test$FR58_DE <- pos_list$FR58[iv,]$Fraction
  
  iv <- match(test$gene, pos_list$FR912$gene)
  test$FR912_DE <- pos_list$FR912[iv,]$Fraction
  
  iv <- match(test$gene, pos_list$FR1619$gene)
  test$FR1619_DE <- pos_list$FR1619[iv,]$Fraction
  
  iv <- match(test$gene, pos_list$FR2326$gene)
  test$FR2326_DE <- pos_list$FR2326[iv,]$Fraction
  
  # iv <- match(test$gene, pos_list$FR3033$gene)
  # test$FR3033_DE <- pos_list$FR3033[iv,]$Fraction
  
  
  # test order of least shared to most shared
  #test3 <- test2 
  #test4 <- order(test3, na.last=FALSE)
  #test5 <- test3[test4,]
  
  
  test2 <- test
  rownames(test2)<- test2$gene 
  test2$gene <- NULL
  
  # Arrange genes with the order of fraction identified (i.e. FR14_DE first)
  # LG has all fractions
  # MM doesn't have fraction 3033_DE
  # LV doesn't have fraction 1619_DE, 2326_DE, 3033_DE
  test3 <- test %>%
    arrange(FR14_DE, FR58_DE, FR912_DE, FR1619_DE, FR2326_DE)
  # test3_MM <- test %>%
  #   arrange(FR14_DE, FR58_DE, FR912_DE, FR1619_DE, FR2326_DE)
  # test3_LV <- test %>%
  #   arrange(FR14_DE, FR58_DE, FR912_DE)
  # test3 <- test2 %>%
  #   arrange(FR3033_DE, FR2326_DE, FR1619_DE, FR912_DE, FR58_DE, FR14_DE)

  gene_level <- test3$gene 
  #test2$gene <- factor(test2$gene, levels = unique(test2$gene))
  plot_order2 <- plot1[gene_level,]

  
  rownames(test3) <- test3$gene
  test3$gene <- NULL 
  
  clusterHeat2 <- pheatmap(plot_order2, show_rownames=F, 
                          breaks = breaks,
                          #annotation_row = my_gene_col,
                          annotation_row = test3,
                          annotation_names_row = TRUE,
                          annotation_names_col = FALSE,
                          annotation_colors = ann_colors,
                          cluster_rows = FALSE, 
                          cluster_cols = FALSE,
                          redgreen(70),
                          annotation_col = md_4, 
                          scale = "row", 
                          fontsize_row=4, fontsize_col=6, fontsize=8,
                          border_color=NA
  )
  save_pheatmap_pdf <- function(x, filename) {
    stopifnot(!missing(x))
    stopifnot(!missing(filename))
    pdf(filename, height = 4.5)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
  }
  write.table(test3, file=paste0("./Files/gene_list_protein_coding_unnorm_deseq_norm_HDvs", Target, "_cluster_up.txt",sep=""), sep="\t", quote=F) 
  save_pheatmap_pdf(clusterHeat2, paste0("./Figures/Selective packaging heatmap/Figure7A_protein_coding_", Target, "_", Baseline, ".pdf"))
 
}
#cluster_heatmap(Target= "LG", Baseline ="HD")
#cluster_heatmap(Target= "LV", Baseline ="HD")
cluster_heatmap(Target= "MM", Baseline ="HD")
