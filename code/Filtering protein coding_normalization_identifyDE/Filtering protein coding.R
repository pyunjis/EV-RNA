# Filtering unnorm count for protein coding

# Set working directory
setwd("E:/4. EV-RNA/4. EV-RNA/4. Code/2023-07-23 FINAL/Files/")

annoFile <- get(load("hg38.Ens_94.biomaRt.geneAnno.Rdata"))

# Select biotype
unique(anno$transcript_biotype)
biotypes <- "protein_coding"

# Set mito <-1 to remove mitochondrial genes
# (set mito <- 0 if you do not want to remove mitochondrial genes)
mito <- 0 

##----------load counts------------#
print("Loading counts table")
#print(countsFile)

counts <- read.table("set1df_unnorm.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t",
                     check.names = FALSE)
##----------load anno------------#
print("Loading annotation table")
print(annoFile)

if(strsplit(biotypes, split='\\,')[[1]]!=""){
  anno.sub <- anno[paste(anno$gene_biotype) %in% strsplit(biotypes, split='\\,')[[1]] ,]
  #counts.sub <- counts[paste(counts$Genes) %in% unique(paste(anno.sub$external_gene_name)) , ]
  counts.sub <- counts[paste(rownames(counts)) %in% unique(paste(anno.sub$ensembl_gene_id)) , ]
}else{
  print("no biotypes provided")
  counts.sub <- counts
}

if(mito==1){
  print("tossing MT- genes")
  counts.sub <- counts.sub[grep("^MT-", paste(counts.sub$Genes), invert=TRUE), ]
}

# Add ERCC back to protein_coding table
# Grab ERCC
spikes <- rownames(counts)[grep("^ERCC", rownames(counts))]
length(spikes)

keep <- rownames(counts)[rownames(counts) %in% spikes]
ercconly <- counts[keep,]

unnorm_ercc <- rbind(counts.sub, ercconly)
write.table(unnorm_ercc, file="set1df_unnorm_protein_coding.txt", sep="\t", quote=F) 
