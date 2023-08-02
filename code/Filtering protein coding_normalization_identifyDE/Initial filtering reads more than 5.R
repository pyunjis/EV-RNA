# Filtering total cell-free RNA/spike-ins 
# with more than five reads in at least one fraction per sample


# Load Library
library(edgeR)
library(RUVSeq)
library(reshape2)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(tidyverse)
library(RColorBrewer)

# Set working directory 

setwd("E:/4. EV-RNA/4. EV-RNA/4. Code/2023-07-23 FINAL/")
# Define tables to read-in for 2ml
cts <- read.table("./Files/NOVO20200421AK_counts.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t",
                  check.names = FALSE)

md <- read.table("./Files/metadata_20200421.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t",
                 check.names = FALSE)

rownames(md) <- md$SampleID
rownames(cts) <- cts[,1]
cts[,1] <- NULL

# Filter out non-expressed genes, by requiring more than 5 counts in 28 samples (one fraction out of 6 fractions per sample)
filter <- apply(cts, 1, function(x) length(x[x>5])>=28) # considering genes present in at least plasma and EV replicates
filtCounts <- cts[filter,]
filtered <- filtCounts + 1

# Grab the rownames of our filtered matrix for ENS gene IDs or gene names
genes <- rownames(filtered)[grep("^ENS", rownames(filtered))]
spikes <- rownames(filtered)[grep("^ERCC", rownames(filtered))]

length(genes) #12671
length(spikes) #55

# Set1 = unnormalized cts
set1 <- newSeqExpressionSet(as.matrix(filtered), phenoData = md)

# Take set1 seqexpression set into dataframe, and take only ERCC for plotting RLE and PCA
set1ct <- counts(set1)
set1df <- as.data.frame(set1ct)
write.table(set1df, file="./Files/set1df_unnorm.txt", sep="\t", quote=F)
