# Figure 2E biotype
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(dplyr)
library(data.table)
library(grid)

# Set working directory
setwd("E:/Box drive update/4. EV-RNA/4. Code/2022-02-17 V7/")

# Import Data
anno <- get(load("Files/hg38.Ens_94.biomaRt.geneAnno.Rdata"))
counts <- read.table("Files/set1df_unnorm.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t",
                     check.names = FALSE, row.names = 1)
coverage <- read.delim("Files/read_distribution/NOVO20200421AK_read_coverage.txt", stringsAsFactors = F)
md <- read.table("Files/metadata_20200421.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t",
                 check.names = FALSE)
# Grab the rownames of our filtered matrix for ENS vs ERCC
genes <- rownames(counts)[grep("^ENS", rownames(counts))]
spikes <- rownames(counts)[grep("^ERCC", rownames(counts))]

length(genes) #12671
length(spikes) #55

# Remove ERCC and consider only genes to test enrichment
keep <- rownames(counts)[rownames(counts) %in% genes]
filtered <- counts[keep,]

# Filter out precancer
md <- filter(md, Condition == "HD" | Condition == "LG" | Condition == "LV" | Condition == "MM")
rownames(md) <- md$SampleID

# keep only relevant SampleID in your count table
keep <- colnames(filtered)[colnames(filtered) %in% rownames(md)]
filtered <- filtered[,keep]  # no log2 counts


################################################
# Annotation of biotypes for gene counts
################################################

# Remove unique gene identifiers from counts table
fD <- as.data.frame(rownames(filtered))
names(fD) <- "ensID"
rownames(fD) <- fD$ensID
iv <- match(rownames(fD), anno$ensembl_gene_id)
fD$biotype <- anno[iv,]$gene_biotype

# generate new dataframe expr which includes gene information
# here presence or absence of gene expression is denoted by a 1/0
# Stacked barplot with all other biotypes marked as "other"
expr <- filtered
expr <- expr[rowSums(expr)>0,]
expr <- expr!=1
expr[expr=="FALSE"] <- 0
expr <- as.data.frame(expr)
expr$ensID <- rownames(expr)

# Add geneID and biotype information to this matrix
m <- match(rownames(expr), rownames(fD))
expr <- left_join(expr, fD[m,])
rownames(expr) <- expr$ensID

# Melt this matrix so it can be plotted
geneTable <- melt(expr, id=c("ensID","biotype"))
geneTable <- geneTable[geneTable$value>0,] # values with count

# Let's filter out some biotypes so that we don't have such a messy graph
ordered <- as.data.table(geneTable)[, .N, by = biotype][order(-N)]
# Grab the first 6 options
filtBio <- as.vector(ordered[6:1,]$biotype) # 6 types: Processed_transcript, transcribed_unprocessed_pseudogene, lincRNA, processed_pseudogene, lincRNA, antisesnse, protein_coding
other <- as.vector(ordered[7:nrow(ordered)]$biotype) # 21 others

# Filter for top 6 biotypes of interest
# Assign all biotypes that are not in the top 6 to be defined as "other"
geneTable$biotype[geneTable$biotype %in% other] <- "other" # Others were assigned for other 21 types
geneTable$biotype <- factor(geneTable$biotype, levels=c("other", filtBio)) # factor so that our graph is in a comprehensive order

# Sum counts across all biotype/sampleID combinations
sampleTable <- as.data.table(geneTable)[, .N, by= .(variable)] # Sample ID with N of genes
geneTable <- as.data.table(geneTable)[, .N, by = .(biotype, variable)] # 7 biotypes with sample ID

# Add fraction imformation
g <- match(sampleTable$variable, md$SampleID)
sampleTable$Type <- md[g,]$Fraction

# Calculate a fraction for each biotype per sample
s <- match(geneTable$variable, sampleTable$variable) # 168 x 7 = 1176 match
geneTable$Frac <- geneTable$N/sampleTable[s,]$N # count for each biotype / 12671 (number of genes)

# Format metadata ad add Condition_FR info into geneTable
rownames(md)<- md$SampleID
c <- match(geneTable$variable, rownames(md))
geneTable$Type <- paste(md[c,]$Condition_FR) # Condition_FR as type

biotypeBar2 <- ggplot(geneTable, aes(x=variable, y=Frac, fill=biotype)) +
  geom_bar(stat="identity") +
  ylab("Fraction of each biotype") +
  xlab("SampleID") +
  scale_fill_brewer( palette = "YlGnBu" ) +
  theme(#text = element_text(family = "Arial"),
    plot.margin = unit(c(0.3,0.3,0.3,0.3),"cm"),
    axis.text.x=element_text(angle = 90, size = 25),
    axis.text.y=element_text(size=30),
    axis.title.y=element_text(size=40),
    axis.title.x=element_text(size=40),
    legend.title=element_blank(),
    legend.text=element_text(size=30),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"))
biotypeBar2

# Mean, min, and max of each biotype
protein_coding <- geneTable %>%
  filter(biotype == "protein_coding")
trans_pseudo <- geneTable %>%
  filter(biotype == "transcribed_unprocessed_pseudogene")
process_trans <- geneTable %>%
  filter(biotype == "processed_transcript")
process_pseudo <- geneTable %>%
  filter(biotype == "processed_pseudogene")
lincRNA <- geneTable %>%
  filter(biotype == "lincRNA")
antisense <- geneTable %>%
  filter(biotype == "antisense")
other <- geneTable %>%
  filter(biotype == "other")

mean(protein_coding$Frac)
min(protein_coding$Frac)
max(protein_coding$Frac)
sd(protein_coding$Frac)

mean(trans_pseudo$Frac)
min(trans_pseudo$Frac)
max(trans_pseudo$Frac)

mean(process_trans$Frac)
min(process_trans$Frac)
max(process_trans$Frac)

mean(trans_pseudo$Frac)
min(trans_pseudo$Frac)
max(trans_pseudo$Frac)

mean(process_pseudo$Frac)
min(process_pseudo$Frac)
max(process_pseudo$Frac)

mean(lincRNA$Frac)
min(lincRNA$Frac)
max(lincRNA$Frac)

mean(antisense$Frac)
min(antisense$Frac)
max(antisense$Frac)

mean(other$Frac)
min(other$Frac)
max(other$Frac)


pdf("Figures/Figure_2E_biotypes.pdf",width = 40, height = 10)


#jpeg("Figures/Figure_2E_biotypes.jpg", width = 4000, height = 1000)
print(biotypeBar2)
dev.off()
