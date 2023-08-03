# Identify DE genes using pair-wise comparisons

# Load libraries
library(DESeq2)
library(dplyr)
library(pheatmap)
library(ggplot2)
library(gplots)


setwd("E:/4. EV-RNA/4. EV-RNA/4. Code/2023-07-23 FINAL/")

# Define tables to read-in 
cts <- read.table("Files/set1df_unnorm_protein_coding.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t",
                  check.names = FALSE)

md <- read.table("Files/metadata_20200421.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t",
                 check.names = FALSE)

# Remove low quality samples EV065 and EV066
cts = cts[, -which(names(cts) %in% c("EV065","EV066"))]
md <- md[ -which(md$SampleID %in% c("EV065","EV066")),]
rownames(md) <- md$SampleID

# Filter for types of interest
RunHeatmap <- function(Baseline, Target){
  md <- filter(md, Condition_FR == Baseline | Condition_FR == Target )

  # now that we have filtered, we can use our sampleIDs as the rownames
  rownames(md) <- md$SampleID
  md <- md[,c(1,2,10)]
  
  # Keep only the relevant sampleIDs in your counts table
  keep <- colnames(cts)[colnames(cts) %in% rownames(md)]
  cts <- cts[, keep]
  dim(cts)

  # Check md too
  md <- md[colnames(cts),]
  dim(md)
  md$SampleID <- NULL

  # Deseq 
  dds_design <- "Condition_FR" # Variable you want to compare

  # Create dds object from counts data and correct columns
  dds <- DESeqDataSetFromMatrix(countData=cts,
                                colData=md,
                                design= as.formula(paste('~',dds_design)))
  
  # Extract normalized count
  dds <- estimateSizeFactors(dds, controlGenes=11610:11664) 
 
  # Remove uninformative columns
  dds <- dds[ rowSums(counts(dds)) >= 1, ]
  
  # Normalization and pre-processing
  dds <- DESeq(dds)

  # results for pair-wise comparison
  res = results(dds, contrast=c("Condition_FR", Target, Baseline),  independentFiltering = FALSE, cooksCutoff = Inf)
 
  # shrink fold changes for lowly expressed genes
  res = lfcShrink(dds, contrast=c("Condition_FR", Target, Baseline), res = res)
  
  # Order by FDR
  results = as.data.frame(res[order(res$padj),])
  numSig <- sum(results$padj < 0.05 & results$log2FoldChange > 1, na.rm=TRUE)
  DEGenes <- results[results$padj <0.05 & results$log2FoldChange>1 ,]
  DEGene_list <- rownames(DEGenes)
  
  # write DESEq2 results to table
  #dir.create("./DESeq2_fraction/")
  filename = paste0("./Files/DESeq2_fraction/unnorm_",Baseline,"_vs_",Target,"_deseq2_results.txt")
  write.table(results, file=filename, sep = "\t", quote = F)

  # Obtain normalized counts (check if normalize by size factor)
  rld <- rlog(dds, blind=FALSE)
  
  # Reduce assay(rld) to include only the overlapping genes
  size <- length(DEGene_list)
  
  plot1 <- assay(rld)[DEGene_list,]

  # making a cluster column to match the row of datac
  datad <- as.data.frame(plot1)
  md$SampleID <- rownames(md) 
  md2 <- md %>%
    arrange(Condition_FR)
  ordered_list<- md2$SampleID
  plot2 <- datad[,ordered_list]
  rownames(md2) <- md2$SampleID
  md2$SampleID <- NULL

  title <- paste0(size," DE genes " ,Baseline," vs ",Target,"")
 
  heatmap1 <- pheatmap(plot2, show_rownames=F, 
                       cluster_cols = FALSE,
                       cluster_rows = FALSE, 
                       annotation_col = md2, scale = "row", border_color = "NA",
                       main = title,  col = redgreen(50),
                       fontsize_row=8, fontsize_col=8, fontsize=8)
 
  # pdf(paste0("./Figures/Heatmap/Heatmap_unnorm_",size, " DE ",Baseline,"_",Target,".pdf"))
  # print(heatmap1)
  # dev.off()

  
}

RunHeatmap(Baseline = "HD_FR14", Target= "LV_FR14")
RunHeatmap(Baseline = "HD_FR58", Target= "LV_FR58")
RunHeatmap(Baseline = "HD_FR912", Target= "LV_FR912")
RunHeatmap(Baseline = "HD_FR1619", Target= "LV_FR1619")
RunHeatmap(Baseline = "HD_FR2326", Target= "LV_FR2326")
RunHeatmap(Baseline = "HD_FR3033", Target= "LV_FR3033")

RunHeatmap(Baseline = "HD_FR14", Target= "LG_FR14")
RunHeatmap(Baseline = "HD_FR58", Target= "LG_FR58")
RunHeatmap(Baseline = "HD_FR912", Target= "LG_FR912")
RunHeatmap(Baseline = "HD_FR1619", Target= "LG_FR1619")
RunHeatmap(Baseline = "HD_FR2326", Target= "LG_FR2326")
RunHeatmap(Baseline = "HD_FR3033", Target= "LG_FR3033")

RunHeatmap(Baseline = "HD_FR14", Target= "MM_FR14")
RunHeatmap(Baseline = "HD_FR58", Target= "MM_FR58")
RunHeatmap(Baseline = "HD_FR912", Target= "MM_FR912")
RunHeatmap(Baseline = "HD_FR1619", Target= "MM_FR1619")
RunHeatmap(Baseline = "HD_FR2326", Target= "MM_FR2326")
RunHeatmap(Baseline = "HD_FR3033", Target= "MM_FR3033")
