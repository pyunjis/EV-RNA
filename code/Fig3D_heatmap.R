# Figure 3D_heatmap

# Load Library
library(RUVSeq)
library(pheatmap)
library(RColorBrewer)
library(gplots)
library(dplyr)

# Set working directory 
setwd("E:/Box drive update/4. EV-RNA/4. Code/2022-02-17 V7/")
# Define tables to read-in
set4df <- read.table("Files/set1df_unnorm_protein_coding_deseq_norm.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t",
                     check.names = FALSE)
md <- read.table("Files/metadata_20200421.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t",
                 check.names = FALSE)
cols <- read.delim("Files/colors.txt", stringsAsFactors = F)
rownames(md) <- md$SampleID

# Remove low quality samples EV065 and EV066
#set4df = set4df[, -which(names(set4df) %in% c("EV065","EV066"))]
md <- md[ -which(md$SampleID %in% c("EV065","EV066")),]

# Remove precancer types
md <- filter(md, Condition == "HD" | Condition == "LG" | Condition == "LV" | Condition == "MM")
rownames(md) <- md$SampleID


# Keep only the gene names and exclude the ERCC transcripts
keep <- rownames(set4df)[grep("^ENS", rownames(set4df))]
filtered <- set4df[keep,]

# Make into a set
set4 <- newSeqExpressionSet(as.matrix(filtered), phenoData = md)

md$Fraction <- factor(md$Fraction, levels=c("FR14", "FR58", "FR912", "FR1619", "FR2326", "FR3033"))
x <- as.factor(md$Fraction)

# Empirical
design <- model.matrix(~x, data=pData(set4)) 
y <- DGEList(counts=counts(set4), group =x, lib.size = rep(c(1),each = nrow(md))) #from normalized set3 ***
y <- calcNormFactors(y, method ="none") # change calcNormFactor method = upperquartile to none or RLE
y <- estimateGLMCommonDisp(y, design, verbose = TRUE) # dispersion vs biological dispersion sq root dispersion gives CV
y <- estimateGLMTagwiseDisp(y, design) # therefore y contains norm factor = 1

# Get colour palette
colors <- brewer.pal(8, "Set1")

# Transform normalized counts to log form
logcounts<- cpm(y, normalized.lib.sizes = FALSE, log = TRUE)

# Find variable genes
var_genes <- apply(logcounts,1, var)
select_var <- names(sort(var_genes, decreasing=TRUE)) # include all genes
highly_variable_lcpm <- logcounts[select_var,]

# Format metadata
md <- md[,c(2,4)]

# Generate palette and map this to each Condition, and to each Fraction
fractionCols <- brewer.pal(length(unique(md$Fraction)), "Set1")
names(fractionCols) <- unique(md$Fraction)

## Remove precancer types or not
cols <- filter(cols, Group == "HD" | Group == "LG" | Group == "LV" | Group == "MM")
conditionCols <- cols$Colour

names(conditionCols) <- cols$Group
annotColors <- list(Fraction = fractionCols, Condition = conditionCols)

# no clustering by column (organize across fraction)
md$Fraction <- factor(md$Fraction, levels=c("FR14", "FR58", "FR912", "FR1619", "FR2326", "FR3033"))
md$Condition <- factor(md$Condition, levels = cols$Group)
md$sampleID <- rownames(md)
md2_arrange <- md %>%
  arrange(Fraction)
rownames(md2_arrange) <- md2_arrange$sampleID
sample_order <- rownames(md2_arrange)
highly_variable <- highly_variable_lcpm[,sample_order]
md2_arrange$sampleID <- NULL
breaksList = seq(-4.5, 4.5, by = 0.1)

hm_2 <- pheatmap(highly_variable, trace = "none", main = "All genes (n=11,609) by variance across fractions", show_rownames=F, show_colnames = F,
                 annotation_col = md2_arrange, cluster_cols = FALSE, scale = "row", col=redgreen(90), annotation_colors = annotColors,
                 breaks =breaksList, fontsize_row=8, fontsize_col=5, fontsize=8)

save_pheatmap_pdf <- function(x, filename) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, height = 4)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_pdf(hm_2, "Figures/Figure_3D_heatmap.pdf")
