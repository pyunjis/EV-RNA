# Normalize using control genes from deseq

# Load library
library(DESeq2)
library(dplyr)
library(pheatmap)
library(ggplot2)
library(gplots)
library(reshape2)
library(RColorBrewer)

# Set working directory
setwd("E:/4. EV-RNA/4. EV-RNA/4. Code/2023-07-23 FINAL/")

# Import Data
cts <- read.table("./Files/set1df_unnorm_protein_coding.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t",
                  check.names = FALSE)
md <- read.table("./Files/metadata_20200421.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t",
                 check.names = FALSE)


# Remove low quality samples EV065 and EV066
cts = cts[, -which(names(cts) %in% c("EV065","EV066"))]
md <- md[ -which(md$SampleID %in% c("EV065","EV066")),]


# Remove precancer types
md <- filter(md, Condition == "HD" | Condition == "LG" | Condition == "LV" | Condition == "MM")
rownames(md) <- md$SampleID

# Keep only relevant SampleIDs in count table
keep <- colnames(cts)[colnames(cts) %in% rownames(md)]
cts <- cts[,keep]

# What attribute of metadat data will we be running DE on
dds_design <- "Condition_FR" # Variable you want to compare


# Create dds object from counts data and correct columns
dds <- DESeqDataSetFromMatrix(countData=cts,
                              colData=md,
                              design= as.formula(paste('~',dds_design)))

# Extract normalized count
dds <- estimateSizeFactors(dds, controlGenes=11610:11664)
# Check size factor
sizeFactors(dds)

# Remove uninformative columns
dds <- dds[ rowSums(counts(dds)) >= 1, ]

# Normalization and pre-processing
dds <- DESeq(dds)

# obtain normalized count from calculating counts divided by size factor
deseq_norm2 <- counts(dds, normalized = TRUE)
dim(deseq_norm2)

#write.table(deseq_norm2, file="./df/set1df_unnorm_protein_coding_deseq_norm.txt", sep="\t", quote=F) 
#write.table(deseq_norm2, file="./df/2022-01-27/set1df_unnorm_protein_coding_deseq_norm.txt", sep="\t", quote=F) 
write.table(deseq_norm2, file="./Files/set1df_unnorm_protein_coding_deseq_norm.txt", sep="\t", quote=F) 

