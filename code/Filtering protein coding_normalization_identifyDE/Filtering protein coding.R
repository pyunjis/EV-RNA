# Filtering unnorm count for protein coding

# Set working directory
setwd("D:/Box drive update/4. EV-RNA/4. Code/2022-02-17 V7/")

countsFile <- "set1df_unnorm"

annoFile <- get(load("Files/hg38.Ens_94.biomaRt.geneAnno.Rdata"))

# Select biotype
unique(anno$transcript_biotype)
biotypes <- "protein_coding"

# Set mito <-1 to remove mitochondrial genes
# (set mito <- 0 if you do not want to remove mitochondrial genes)
mito <- 0 
##----------load counts------------#
print("Loading counts table")
#print(countsFile)

## check if an rda file or tab sep
if(grepl('rda|RData|Rdata',countsFile)){
  counts <- get(load(file=countsFile))
}
if(grepl('txt|tsv',countsFile)){
  counts <- read.delim(file=countsFile)
}

##----------load counts------------#
print("Loading counts table")
print(countsFile)

## must be a tsv or txt tab sep file

#counts <- read.delim(file=countsFile)
counts <- read.table("./Files/set1df_unnorm.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t",
                     check.names = FALSE)
rownames(counts) <- counts[,1]
counts[,1] <- NULL
##----------load anno------------#
print("Loading annotation table")
print(annoFile)

## check if an rda file or tab sep
if(grepl('rda|RData|Rdata',annoFile)){
  anno <- get(load(file=annoFile))
}
if(grepl('txt|tsv',annoFile)){
  anno <- read.delim(file=annoFile)
}

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

#write.table(counts.sub, file=sub(".txt", ".filt.txt", countsFile), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
#write.table(counts.sub, file="Files/set1df_unnorm_protein_coding.txt", sep="\t", quote=F) 
 
# Add ERCC back to protein_coding table
# Grab ERCC
spikes <- rownames(counts)[grep("^ERCC", rownames(counts))]
length(spikes)

keep <- rownames(counts)[rownames(counts) %in% spikes]
ercconly <- counts[keep,]

unnorm_ercc <- rbind(counts.sub, ercconly)
write.table(unnorm_ercc, file="Files/set1df_unnorm_protein_coding.txt", sep="\t", quote=F) 


