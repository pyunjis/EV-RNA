# Identify DE genes using pair-wise comparisons

# Load libraries
library(DESeq2)
library(dplyr)
library(pheatmap)
library(ggplot2)
library(gplots)


setwd("E:/4. EV-RNA/4. EV-RNA/4. Code/2023-07-23 FINAL/")

# Define tables to read-in for 2ml
cts <- read.table("Files/set1df_unnorm_protein_coding.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t",
                  check.names = FALSE)

md <- read.table("Files/metadata_20200421.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t",
                 check.names = FALSE)
#ercc <- read.table("./others/ERCC_concentration.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t", check.names = FALSE)


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

  # Deseq needs integer, not decimals
  #int_cts <- cts
  #int_cts[] <- lapply(int_cts, round)

  dds_design <- "Condition_FR" # Variable you want to compare

  
  # Create dds object from counts data and correct columns
  dds <- DESeqDataSetFromMatrix(countData=cts,
                                colData=md,
                                design= as.formula(paste('~',dds_design)))
  
  # Extract normalized count
  dds <- estimateSizeFactors(dds, controlGenes=11610:11664) # was originally 11609
  #sizeFactors(dds)
  # Remove uninformative columns
  dds <- dds[ rowSums(counts(dds)) >= 1, ]
  
  # Normalization and pre-processing
  dds <- DESeq(dds)

  # results for pair-wise comparison
  res = results(dds, contrast=c("Condition_FR", Target, Baseline),  independentFiltering = FALSE, cooksCutoff = Inf)
  # shrink fold changes for lowly expressed genes
  #res = lfcShrink(dds, contrast=c("Condition_FR", Target, Baseline), res=res) # check again
  res = lfcShrink(dds, contrast=c("Condition_FR", Target, Baseline), res = res)
  # Likelihood Ratio test to look at differential expression across ALL types, and not just pairs of types (contrast)
  #dds.lrt <- DESeq(dds, test="LRT", reduced=~1)
  #res.lrt <- results(dds.lrt, cooksCutoff = Inf, independentFiltering=FALSE)

  #res <- results(dds.lrt)
  #topGenes <- head(order(res$padj), 50)
  
  
  # Order by FDR
  results = as.data.frame(res[order(res$padj),])
  numSig <- sum(results$padj < 0.05 & results$log2FoldChange > 1, na.rm=TRUE)
  DEGenes <- results[results$padj <0.05 & results$log2FoldChange>1 ,]
  DEGene_list <- rownames(DEGenes)
  # write DESEq2 results to table
  #dir.create("./DESeq2_fraction/")
  #filename = paste0("./DESeq2_fraction/2021-10-19/unnorm_",Baseline,"_vs_",Target,"_deseq2_results.txt")
  filename = paste0("./Files/DESeq2_fraction2/unnorm_",Baseline,"_vs_",Target,"_deseq2_results.txt")
  write.table(results, file=filename, sep = "\t", quote = F)
  #filename2 = paste0("./Files/Biological relevance/DESEQ_FDR_order/unnorm_",Baseline,"_vs_",Target,"_deseq2_padj0.05_results.txt")
  #write.table(DEGenes, file=filename2, sep = "\t", quote = F)


  # Obtain normalized counts (check if normalize by size factor)
  rld <- rlog(dds, blind=FALSE)
  
  

  # Reduce assay(rld) to include only the overlapping genes
  # plot1 <- assay(rld.lrt)[overlap,]
  size <- length(DEGene_list)
  #number_genes2 <- paste0(number_genes)
  
  plot1 <- assay(rld)[DEGene_list,]
  
  #plot1 <- assay(rld)[topGenes,]

  ## Replace ensemble id's with gene id's
  #gene_id = read.delim("/Users/khyu/Desktop/programming/geneCounts/biomart_ensembl_geneid.txt")
  #anno <- get(load("/Users/khyu/Desktop/Programming/20190812/geneCounts_20190812/hg38.Ens_94.biomaRt.geneAnno.Rdata"))

  ## Remove unique identifier .xx from heatmap data
  # rownames(plot1) <- sub("\\.[0-9]*", "", rownames(plot1))
  # iv <- match(rownames(plot1), gene_id$ensembl_gene_id)
  # head(gene_id[iv,])

  ## assign the rownames of plot_LG to their external gene name using a variable, where we have indexed the 
  ## row number for all matches between these two dataframes
  ## Use paste to get rid of factors of this column, and just paste the value of the gene name
  # rownames(plot1) <- paste(gene_id[iv, "external_gene_name"])
  
  # making a cluster column to match the row of datac
  datad <- as.data.frame(plot1)
  md$SampleID <- rownames(md) 
  md2 <- md %>%
    arrange(Condition_FR)
  ordered_list<- md2$SampleID
  plot2 <- datad[,ordered_list]
  rownames(md2) <- md2$SampleID
  md2$SampleID <- NULL
  #row.names(md) <- colnames(plot1)
  title <- paste0(size," DE genes " ,Baseline," vs ",Target,"")
  #title <- paste0("Top 50 FDR ",Baseline," vs ",Target,"")
  heatmap1 <- pheatmap(plot2, show_rownames=F, 
                       #clustering_distance_cols = "correlation", clustering_method = "ward.D2",
                       cluster_cols = FALSE,
                       cluster_rows = FALSE, #clustering_distance_rows = "correlation", 
                       #clustering_distance_cols = "correlation", clustering_method = "average", 
                       annotation_col = md2, scale = "row", border_color = "NA",
                       main = title,  col = redgreen(50),
                       #cellwidth =11, cellheight=7,
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

# RunHeatmap(Baseline = "HD_FR14", Target= "LCir_FR14")
# RunHeatmap(Baseline = "HD_FR58", Target= "LCir_FR58")
# RunHeatmap(Baseline = "HD_FR912", Target= "LCir_FR912")
# RunHeatmap(Baseline = "HD_FR1619", Target= "LCir_FR1619")
# RunHeatmap(Baseline = "HD_FR2326", Target= "LCir_FR2326")
# RunHeatmap(Baseline = "HD_FR3033", Target= "LCir_FR3033")
# 
# RunHeatmap(Baseline = "HD_FR14", Target= "MG_FR14")
# RunHeatmap(Baseline = "HD_FR58", Target= "MG_FR58")
# RunHeatmap(Baseline = "HD_FR912", Target= "MG_FR912")
# RunHeatmap(Baseline = "HD_FR1619", Target= "MG_FR1619")
# RunHeatmap(Baseline = "HD_FR2326", Target= "MG_FR2326")
# RunHeatmap(Baseline = "HD_FR3033", Target= "MG_FR3033")

###############################################################################



