library(DESeq2)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)

# Set working directory
setwd("D:/Box drive update/4. EV-RNA/4. Code/2022-02-17 V7/")

permutation_test <- function(Type1,Type2){ 
  counts <- read.table("./Files/set1df_unnorm_protein_coding.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t",
                  check.names = FALSE)

  md <- read.table("./Files/metadata_20200421.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t",
                   check.names = FALSE)
  
  contrast <- c(Type1, Type2)
  baseline <- contrast[[2]]
  target <- contrast[[1]]
  
  
  # ensID to rownames for set5df_ournorm_control2
  # rownames(counts) <- counts[,1]
  # counts[,1] <- NULL
  
  # Extract only the Types that we want in further analysis & only the PP_ID and Status informative columns
  md <- filter(md, Condition_FR == baseline | Condition_FR == target)
  
  # For LV065 & LV066 need to be removed
  md <- filter(md, !SampleID %in% c("EV065", "EV066"))
  rownames(md) <- md$SampleID
  md$SampleID <- NULL
  
  # keep only relevant SampleID in your count table
  keep <- colnames(counts)[colnames(counts) %in% rownames(md)]
  subdata <- counts[, keep]
  dim(subdata)
  
  # Deseq needs integer, not decimals
  #subdata <- subdata
  # subdata[] <- lapply(subdata, round)
  
  # Check
  stopifnot(all.equal(rownames(md), colnames(subdata)))
  
  # Get the number of Cancer samples and number of HD samples from md table
  num1 = sum(md$Condition_FR == baseline)
  num2 = sum(md$Condition_FR == target)
  
  # Create a vector for both HD and Can, with a 1 for every HD and a 2 for every Cancer sample
  One_vector = rep(c(1), times = num1)
  Two_vector = rep(c(2), times = num2)
  
  # Permutation
  # Concatenate the HD and Can vector to create your "start group" vector
  start_group = c(One_vector, Two_vector)
  cutoff=0.05 # was originally 0.01
  number_of_diff_genes=c()
  #number_of_diff_genes_up=c()
  group_list = list()
  number_of_try = 1000
  nulldist <- rep(0, number_of_try)
  
  for (i in 1:number_of_try) {
    print(i)
    group = data.frame(type=factor(sample(start_group)))
    
    dds = DESeqDataSetFromMatrix(countData = subdata,
                                 colData = group,
                                 design = ~ type)
    
    # Extract normalized counts
    #dds = estimateSizeFactors(dds)
    dds <- estimateSizeFactors(dds, controlGenes=11610:11664)
    
    # Remove genes with zero counts over all samples
    dds <- dds[ rowSums(counts(dds)) >= 1, ]
    
    # Make sure of reference, set it by rlevel
    dds$type = relevel(dds$type, ref = 1)
    
    # The standard differential expression analysis steps are wrapped into a single function, DESeq
    dds = DESeq(dds)
    
    # Extract results
    res = results(dds, contrast = c("type", "1", "2"), independentFiltering = FALSE,cooksCutoff = Inf)
    res = lfcShrink(dds, contrast=c("type", "1", "2"), res = res)
    #res2 <- as.data.frame(res)
    #num=res2[res2$padj < 0.05 & res2$log2FoldChange>1 ,]
    #tmp=sum(res$padj < cutoff , na.rm=TRUE)
    tmp=sum(res$padj < cutoff & res$log2FoldChange>1 , na.rm=TRUE)
    #sig_genes = res[res$padj<0.05 & res$log2FoldChange>1,]
    
    number_of_diff_genes = c(number_of_diff_genes,tmp)
    #nulldist[i] <- length(sig_genes)
    group_list[[i]] <- group
    
    # Save number of diff_genes
    
    
  }
  
  # Obtain the number of genes that meet padj<0.01 for reference line in histogram
  dds <- DESeqDataSetFromMatrix(countData=subdata,
                                colData=md,
                                design= ~ Condition_FR) # This is Type from md which compared correct HD vs Cx labels
  
  #dds <- estimateSizeFactors(dds)
  dds <- estimateSizeFactors(dds, controlGenes=11610:11664)
  
  # Remove uninformative columns
  dds <- dds[ rowSums(counts(dds)) >= 1, ]
  
  # Normalization and pre-processing
  dds <- DESeq(dds)
  
  # Extract results and the number of significant genes with padj<0.01
  res = results(dds, contrast = c("Condition_FR", target, baseline), independentFiltering = FALSE,cooksCutoff = Inf)
  res = lfcShrink(dds, contrast=c("Condition_FR", target, baseline), res = res)
  #numSig <- sum(results$padj < 0.05  & results$log2FoldChange > 1, na.rm=TRUE)
  numSig <- sum(res$padj < 0.05  & res$log2FoldChange > 1, na.rm=TRUE)
  
  #results2 <- as.data.frame(results)
  #sig_genes <- results2[results2$padj<0.05 & results2$log2FoldChange>1,]
  #numSig2 <- length(rownames(sig_genes))
  
  number_of_diff_genes <- as.data.frame(number_of_diff_genes)
  names(number_of_diff_genes) <- "NumDiffGenes"
  number_of_diff_genes$Actual <- numSig
  
  p <- ggplot(number_of_diff_genes, aes(x=NumDiffGenes)) +
    #xlim(c(0,500))+
    geom_histogram(bins=100) +
    geom_vline(data=number_of_diff_genes, mapping=aes(xintercept = numSig, color = "Correct Labels"), 
               linetype="longdash", size=0.6, show.legend = T) +
    scale_color_manual(values = "gray75", name = "Number of DE genes") +
    ggtitle(paste(number_of_try, "Random Permutations:", baseline, "vs", target)) +
    xlab("Number of significant genes") +
    theme(aspect.ratio=1,
          plot.title = element_text(hjust = 0.5),
          legend.title = element_text(size=10, hjust = 0.5))
  
  
  pdf(paste0("./Figures/Permutation_test/Permutation_test_histogram_set1df_norm_ctrl_gene_", Type1, "_", Type2, "_Actual DE_", numSig, ".pdf"))
  print({
    p
  })
  dev.off()
  
  # df <- data.frame(stringsAsFactors = FALSE)
  # 
  # for (i in 1:number_of_try) {
  #   if (i==1) {
  #     df = group_list[[i]]
  #   }
  #   else {
  #     df = cbind(df, group_list[[i]])
  #   }
  #   colnames(df)[i] = paste("perm",i, sep = "_")
  # }
  # 
  # write.csv(number_of_diff_genes, file = paste0("./Permutation_test/unnorm_norm_ctrl_gene_", Type1, "_", Type2, "_number.diff.genes.csv"))
  # write.csv(df, paste0("./Permutation_test/unnorm_norm_ctrl_gene_", Type1, "_", Type2, "_permutation.list.csv"))
}
permutation_test(Type1 = "LV_FR14", Type2 = "HD_FR14")
permutation_test(Type1 = "LV_FR58", Type2 = "HD_FR58")
permutation_test(Type1 = "LV_FR912", Type2 = "HD_FR912")
permutation_test(Type1 = "LV_FR1619", Type2 = "HD_FR1619")
permutation_test(Type1 = "LV_FR2326", Type2 = "HD_FR2326")
permutation_test(Type1 = "LV_FR3033", Type2 = "HD_FR3033")

permutation_test(Type1 = "MM_FR14", Type2 = "HD_FR14")
permutation_test(Type1 = "MM_FR58", Type2 = "HD_FR58")
permutation_test(Type1 = "MM_FR912", Type2 = "HD_FR912")
permutation_test(Type1 = "MM_FR1619", Type2 = "HD_FR1619")
permutation_test(Type1 = "MM_FR2326", Type2 = "HD_FR2326")
permutation_test(Type1 = "MM_FR3033", Type2 = "HD_FR3033")

permutation_test(Type1 = "LG_FR14", Type2 = "HD_FR14")
permutation_test(Type1 = "LG_FR58", Type2 = "HD_FR58")
permutation_test(Type1 = "LG_FR912", Type2 = "HD_FR912")
permutation_test(Type1 = "LG_FR1619", Type2 = "HD_FR1619")
permutation_test(Type1 = "LG_FR2326", Type2 = "HD_FR2326")
permutation_test(Type1 = "LG_FR3033", Type2 = "HD_FR3033")