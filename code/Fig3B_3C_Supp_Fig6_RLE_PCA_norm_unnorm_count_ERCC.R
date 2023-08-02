# Figure 3B RLE plot
# Figure 3C PCA plot
# Supp Figure 6A-D, unnorm vs norm counts and ERCC

# Load Library
library(RUVSeq)
library(reshape2)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)

# Set working directory 
setwd("E:/4. EV-RNA/4. EV-RNA/4. Code/2023-07-23 FINAL/")

# Define tables to read-in
cts <- read.table("./Files/set1df_unnorm_protein_coding_deseq_norm.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t",
                  check.names = FALSE)
unnorm <- read.table("./Files/set1df_unnorm_protein_coding.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t",
                     check.names = FALSE)
md <- read.table("./Files/metadata_20200421.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t",
                 check.names = FALSE)

# Remove low quality samples
md <- md[ -which(md$SampleID %in% c("EV065","EV066")),]

# Remove precancer types
md <- filter(md, Condition == "HD" | Condition == "LG" | Condition == "LV" | Condition == "MM")
rownames(md) <- md$SampleID

# Remove precancer on unnorm count table
keep <- colnames(unnorm)[colnames(unnorm) %in% rownames(md)]
unnorm <- unnorm[,keep]  # no log2 counts


# Keep only ERCC for normalized cts
keep <- rownames(cts)[grep("^ERCC", rownames(cts))]
cts_ercc <- cts[keep,]

# Keep only genes by removing ERCC for cts
keep <- rownames(cts)[grep("^ENS", rownames(cts))]
cts <- cts[keep,]

# Keep only ERCC for unnorm cts
keep <- rownames(unnorm)[grep("^ERCC", rownames(unnorm))]
unnorm_ercc <- unnorm[keep,]

# Keep only genes by removing ERCC for unnorm cts
keep <- rownames(unnorm)[grep("^ENS", rownames(unnorm))]
unnorm <- unnorm[keep,]


# re-order by Fraction
md$Fraction <- factor(md$Fraction, levels = c("FR14", "FR58", "FR912", "FR1619", "FR2326", "FR3033"))
md <- md %>%
  arrange(Condition) %>%
  arrange(Fraction)

# Making new sample list
md$SampleID <- factor(md$SampleID, levels = md$SampleID)

# Order counts table the same as md
cts <- cts[,levels(md$SampleID)]
unnorm <- unnorm[,levels(md$SampleID)]
cts_ercc <- cts_ercc[,levels(md$SampleID)]
unnorm_ercc <- unnorm_ercc[,levels(md$SampleID)]
rownames(md) <- md$SampleID

# Create set 0 for unnorm cts
set0 <- newSeqExpressionSet(as.matrix(unnorm), phenoData = md)
set0$Fraction <- factor(set0$Fraction, levels = c("FR14", "FR58", "FR912", "FR1619", "FR2326", "FR3033"))
d <- as.factor(set0$Fraction)
colors <- brewer.pal(8, "Set1")

# Create set 1 for norm cts
set1 <- newSeqExpressionSet(as.matrix(cts), phenoData = md)
set1$Fraction <- factor(set1$Fraction, levels = c("FR14", "FR58", "FR912", "FR1619", "FR2326", "FR3033"))
d <- as.factor(set1$Fraction)
colors <- brewer.pal(8, "Set1")

# Plot log expression instead of relative
# make tidy
cts2 <- cts
cts2$gene <- rownames(cts2)
ctsMelt <- melt(cts2, by = "gene") %>%
  mutate(log2counts = log2(value + 1)) 

colnames(ctsMelt)[2:4] <- c("SampleID","counts", "log2counts")

# Add type info
iv <- match(ctsMelt$SampleID, md$SampleID)
ctsMelt$Fraction <- md[iv,]$Fraction


# plotting median expression for all detected genes
a <- ggplot(ctsMelt, aes(x=SampleID, y=log2counts, fill=Fraction)) +
  geom_bar(stat="summary",fun = "median", colour="black", size =0.1, show.legend = FALSE) +
  #geom_boxplot(alpha=0.7, coef=0, fill = colors[d]) +
  ylab("log2(counts+1)") +
  #xlab("SampleID") +
  scale_fill_brewer( palette = "Set1" ) +
  theme(#text = element_text(family = "Arial"),
    plot.margin = unit(c(0.3,0.3,0.3,0.3),"cm"),
    axis.text.x=element_text(angle = 90, size = 5),
    axis.text.y=element_text(size=16),
    axis.title.y=element_text(size=16),
    axis.title.x=element_text(size=10),
    legend.title=element_blank(),
    legend.position= "none",
    legend.text=element_text(size=10),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"))
pdf("./Figures/Figure_3B_norm_med_gene_exp.pdf", width = 5, height =4.5)
print(a)
dev.off()

# Find min, max, and median of each fraction
FR14 <- ctsMelt %>%
  filter(Fraction == "FR14")
FR58 <- ctsMelt %>%
  filter(Fraction == "FR58")
FR912 <- ctsMelt %>%
  filter(Fraction == "FR912")
FR1619 <- ctsMelt %>%
  filter(Fraction == "FR1619")
FR2326 <- ctsMelt %>%
  filter(Fraction == "FR2326")
FR3033 <- ctsMelt %>%
  filter(Fraction == "FR3033")

median(FR14$log2counts)
min(FR14$log2counts)
max(FR14$log2counts)
median(FR58$log2counts)
min(FR58$log2counts)
max(FR58$log2counts)
median(FR912$log2counts)
min(FR912$log2counts)
max(FR912$log2counts)
median(FR1619$log2counts)
min(FR1619$log2counts)
max(FR1619$log2counts)
median(FR2326$log2counts)
min(FR2326$log2counts)
max(FR2326$log2counts)
median(FR3033$log2counts)
min(FR3033$log2counts)
max(FR3033$log2counts)

# count how many genes detected with log2expression higher than 10
length(unique(FR14$SampleID)) #20
length(unique(FR58$SampleID)) #20
length(unique(FR912$SampleID)) #20
length(unique(FR1619$SampleID)) #20
length(unique(FR2326$SampleID)) #19
length(unique(FR3033$SampleID)) #19
sum(FR14$log2counts > 10) /20
sum(FR58$log2counts > 10) /20
sum(FR912$log2counts > 10) /20
sum(FR1619$log2counts > 10) /19
sum(FR2326$log2counts > 10) /19
sum(FR3033$log2counts > 10) /19

# plotting median expression of all detected genes
# ggplot(ctsMelt, aes(x=SampleID, y=log2counts)) +
#   geom_boxplot(alpha=0.7, outlier.shape = NA, coef=0, fill = colors[d]) +
#   #stat_summary(fun=median, geom="point", shape=21, size=2, color="black", fill=colors[d]) +
#   theme(legend.position="none") +
#   scale_fill_brewer(palette="Set1") +
#   ylim(-3,15) +
#   #geom_hline(yintercept=0, color= "grey50")+
#   #ylab("RLE") +
#   #ggtitle("RLE of unnormalized ERCC count") +
#   theme(panel.grid.major = element_line(colour = "white", size = 1)) +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
#   theme(strip.background = element_rect(colour = "grey50", fill = "grey50")) +
#   theme(text = element_text(size=15), axis.text.x = element_text(angle=90, hjust=1)) 

# Create set 2 for unnorm_ercc
set2 <- newSeqExpressionSet(as.matrix(unnorm_ercc), phenoData = md)
set2$Fraction <- factor(set2$Fraction, levels = c("FR14", "FR58", "FR912", "FR1619", "FR2326", "FR3033"))
d <- as.factor(set2$Fraction)
colors <- brewer.pal(8, "Set1")

# Create set 3 for norm_ercc
set3 <- newSeqExpressionSet(as.matrix(cts_ercc), phenoData = md)
set3$Fraction <- factor(set3$Fraction, levels = c("FR14", "FR58", "FR912", "FR1619", "FR2326", "FR3033"))
d <- as.factor(set3$Fraction)
colors <- brewer.pal(8, "Set1")


# plotRLE with new order
pdf("./Figures/Supp_Figure_6A_Unnorm_RLE.pdf", width = 5, height =5)
plotRLE(set0, col=colors[d], legend=TRUE, outline=FALSE, ylim=c(-2, 6), ylab ="RLE",
        cex.axis =1.3, cex.lab =1.5, cex.main =2, whisklty = 0, staplelty = 0)
dev.off()

# plotRLE with new order
pdf("./Figures/Supp_Figure_6B_Norm_RLE.pdf", width = 5, height =5)
plotRLE(set1, col=colors[d], legend=TRUE, outline=FALSE, ylim=c(-2, 6), ylab ="RLE",
        cex.axis =1.3, cex.lab =1.5, cex.main =2, whisklty = 0, staplelty = 0)
dev.off()

#las = 2 makes text x-axis rotate 90.

# plotRLE with new order
pdf("./Figures/Supp_Figure_6C_Unnorm_ERCC.pdf", width = 5, height =5)
plotRLE(set2, col=colors[d], legend=TRUE, outline=FALSE, ylim=c(-3, 3), ylab ="RLE",
        cex.axis =1.3, cex.lab =1.5, cex.main =2, whisklty = 0, staplelty = 0)
dev.off()

# plotRLE with new order
pdf("./Figures/Supp_Figure_6D_Norm_ERCC.pdf", width = 5, height =5)
plotRLE(set3, col=colors[d], legend=TRUE, outline=FALSE, ylim=c(-3, 3), ylab ="RLE",
        cex.axis =1.3, cex.lab =1.5, cex.main =2, whisklty = 0, staplelty = 0)
dev.off()

# Plot PCA for norm cts
pdf("./Figures/Figure_3C_Norm_PCA.pdf", width = 5, height = 5)
plotPCA(set1, bg=colors[d], labels= FALSE, pch=21, cex=1.5, oma=c(3,3,3,15), cex.axis =1.2, cex.lab =1.2, cex.main =1.2)
dev.off()
