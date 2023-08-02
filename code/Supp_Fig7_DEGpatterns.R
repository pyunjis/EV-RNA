# Perform degpatterns for Supp Fig 7

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

# Define tables to read-in
cts <- read.table("./Files/set1df_unnorm_protein_coding_deseq_norm.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t",
                  check.names = FALSE)
md <- read.table("./Files/metadata_20200421.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t",
                 check.names = FALSE)

# remove low quality samples
md <- filter(md, !SampleID %in% c("EV065", "EV066"))

# Remove precancer types
md <- filter(md, Condition == "HD" | Condition == "LG" | Condition == "LV" | Condition == "MM")
rownames(md) <- md$SampleID

# Keep only relevant SampleIDs in count table
keep <- colnames(cts)[colnames(cts) %in% rownames(md)]
cts <- cts[,keep]

# type fraction levels
# type <- "HD"
# type_levels<- paste(type, c("FR14", "FR58","FR912", "FR1619", "FR2326", "FR3033"), sep = "_")


# Make into a set
#md$Condition_FR <- factor(md$Condition_FR, levels=c(type_levels))
md$Fraction <- factor(md$Fraction, levels=c("FR14", "FR58", "FR912", "FR1619", "FR2326", "FR3033"))




set <- newSeqExpressionSet(as.matrix(cts), phenoData = md)
x <- as.factor(md$Fraction)

# design to run DGE 
design <- model.matrix(~Fraction, data=pData(set)) # was originally Condition_FR to run per type, now run across Fraction
y <- DGEList(counts=counts(set), group =x, lib.size = rep(c(1), times = 118))#from normalized set3 ***
y <- calcNormFactors(y, method ="none") # change calcNormFactor method = upperquartile to none or RLE
y <- estimateGLMCommonDisp(y, design, verbose = TRUE) # dispersion vs biological dispersion sq root dispersion gives CV
y <- estimateGLMTagwiseDisp(y, design) # therefore y contains norm factor = 1

#colors <- brewer.pal(8, "Set1")
#plotMDS(y, pch=21, col="black", oma=c(3,3,3,15), bg=colors[x], cex=2, main ="MDS", cex.axis =1.3, cex.lab =1.5, cex.main =2)
fit <- glmFit(y, design)
#fit_F <- glmQLFit(y, design)
lrt <- glmLRT(fit, coef=2:6) # This is finding difference between any of 6 groups
res <- as.data.frame(topTags(lrt, n=nrow(y)))
sig_LRT <- res[res$FDR <0.05 ,] # Filter based on significance

logcounts <- log(cts+1, 2)
logcounts2 <- logcounts[rownames(sig_LRT),]

library(DEGreport)
clusters <- degPatterns(logcounts2, metadata = md, time = "Fraction", col=NULL, minc = 0)


#colors <- brewer.pal(8, "Set1")


# What type of data structure is "cluster" output?
class(clusters) # list

# Let's see what is stored in the "df" component
head(clusters$df) # shows gene with cluster components

# Extract the Group 1 genes
cluster_groups <- clusters$df

# Match ensID to obtain FDR from sig_LRT
iv <- match(rownames(cluster_groups), rownames(sig_LRT))
cluster_groups$FDR <- paste(sig_LRT[iv, "FDR"])

# Make into an excel
write.table(cluster_groups, file="./Files/degPatterns/cluster_groups_FDR0.05_degpatterns_ALL_FR_minc0.txt", sep="\t", quote=F) 

# Normalize
norm <- clusters$normalized
iv <- match(norm$genes, rownames(sig_LRT))
norm$FDR <- paste(sig_LRT[iv, "FDR"])


# Make into an excel
write.table(norm, file="./Files/degPatterns/norm_FDR0.05_degpatterns_ALL_FR_minc0.txt", sep="\t", quote=F) 

counts <- table(distinct(norm, genes, cluster)[["cluster"]])
#names(counts) <- factor(names(counts), levels=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26","27"))

table <- inner_join(norm, data.frame(cluster = as.integer(names(counts)), 
                                     tittle = paste(names(counts), "- genes:", counts), 
                                     stringsAsFactors = FALSE), by = "cluster")
# Make into an excel
write.table(table, file="./Files/degPatterns/table_FDR0.05_degpatterns_ALL_FR_minc0.txt", sep="\t", quote=F) 

table2 <- table %>%
  arrange(cluster)
tittle_order <- unique(table2$tittle) # sort_unique defines unique levels
table2$tittle <- factor(table2$tittle, levels = tittle_order)   
color_cluster <- as.factor(table2$cluster)
# Use the data yourself for custom figures
#clusters[["normalized"]] instead of table
a <- ggplot(table2,
       aes(Fraction, value, color = color_cluster))+#, fill = color_cluster)) +
  geom_boxplot(alpha = 0, outlier.size = 0, outlier.shape = NA) +
  facet_wrap(~tittle) +
  geom_point(alpha= 0.4, size =1, position = position_jitterdodge(dodge.width = 0.9)) +
  # change the method to make it smoother
  stat_smooth(aes_string(x = "Fraction", y = "value", 
                         group = color_cluster, color = color_cluster), se = FALSE, method = "lm", 
              formula = y ~ poly(x, splan)) +
  geom_line(aes_string(group = "genes"), 
            alpha = 0.1) +
  stat_summary(fun=mean, geom="line", aes(group=1), size = 0.5) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, 
                                   vjust = 0.5)) + ylab("Z-score of gene abundance") + xlab("") +
  theme(panel.grid.major = element_line(colour = "grey95", size = 1)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  theme(strip.background = element_rect(colour = "grey50", fill = "grey80")) +
  theme(text = element_text(size=12), axis.text.x = element_text(angle=90, hjust=1))


pdf("./Figures/Supp_Figure_7_degpatterns_ALL_minc0.pdf", width = 8, height =6)
print(a)
dev.off()
