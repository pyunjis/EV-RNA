# Figure 2B Mean Input Table 
# Figure 2D number of input reads
# Supp Figure 5A, 5B, and 5C %unique reads, exon fraction, and intron fraction 

# load libraries
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(dplyr)
library(data.table)
library(tidyverse)
library(scales)

setwd("E:/4. EV-RNA/4. EV-RNA/4. Code/2023-07-23 FINAL/")

coverage <- read.table("Files/read_distribution/NOVO20200421AK_read_coverage.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t",
                       check.names = FALSE)
star <- read.table("Files/read_distribution/NOVO20200421AK_star_filter.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t",
                   check.names = FALSE)
md <- read.table("Files/metadata_20200421.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t",
                 check.names = FALSE)
cols <- read.delim("Files/colors.txt", stringsAsFactors = F)
unique <- read.table("Files/read_distribution/unique_reads_unmapped_reads.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t",
                     check.names = FALSE)

# Filter out precancer
md <- filter(md, Condition == "HD" | Condition == "LG" | Condition == "LV" | Condition == "MM")
rownames(md) <- md$SampleID

# keep only relevant SampleID in your count table
filtered <- coverage
rownames(filtered) <- filtered$Sample

# keep only relevant SampleID in your count table
keep <- rownames(filtered)[rownames(filtered) %in% rownames(md)]
filtered_hd <- filtered[keep,]  # no log2 counts
coverage <- filtered_hd

# Add type and fraction information
coverage$Condition <- md$Condition
coverage$Fraction <- md$Fraction

cover_order <- coverage %>%
  gather(types, value, 2:4)

# Factor with a specific order
cover_order$types <- as.character(cover_order$types)
cover_order$types <- factor(cover_order$types, levels=c("Intron", "Intergenic", "Exon"))
cover_order$Fraction <- factor(cover_order$Fraction, levels=c("FR14", "FR58", "FR912", "FR1619", "FR2326", "FR3033"))

dodge <- position_dodge(width = 0.6)
blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA),
    text = element_text(color='black', size=18),
    panel.grid.major.y = element_line(colour = "black"),
    panel.grid.minor.y = element_line(colour = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x=element_text(angle=90, hjust=1, colour="black"),
    axis.text.y=element_text(colour="black"),
    plot.title=element_text(size=18, face="bold", hjust = 0.5)
  )

# Violin plot for intergenic, exon and intron separately
##only intergenic 
cover_order_inter <- cover_order %>%
  filter(types == "Intergenic")

p1 <- ggplot(cover_order_inter, aes(x=Condition, y=value, fill=Condition)) +
  geom_violin() +
  geom_boxplot(width=.3, outlier.colour=NA, position= dodge, color="black") +
  blank_theme +
  scale_fill_manual(breaks=cols$Group,
                    values=cols$Colour) +
  ylab("Intergenic Fraction") +
  facet_grid(.~Fraction) 

##only intron
cover_order_intron <- cover_order %>%
  filter(types == "Intron")

p2 <- ggplot(cover_order_intron, aes(x=Condition, y=value, fill=Condition)) +
  geom_violin()+
  geom_boxplot(width=.1, outlier.colour=NA, position= dodge, color="black")+
  blank_theme +
  scale_fill_manual(breaks=cols$Group,
                    values=cols$Colour) +
  ylab("Intron Fraction") +
  facet_grid(.~Fraction)

# Get mean, min, and maximum of intron fraction
mean(cover_order_intron$value)
min(cover_order_intron$value)
max(cover_order_intron$value)
# Get mean, min, and maximum of intergenic fraction
mean(cover_order_inter$value)
min(cover_order_inter$value)
max(cover_order_inter$value)

##only Exon
cover_order_exon <- cover_order %>%
  filter(types == "Exon")

p3 <- ggplot(cover_order_exon, aes(x=Condition, y=value, fill=Condition)) +
  geom_violin()+
  geom_boxplot(width=.1, outlier.colour=NA, position= dodge, color="black")+
  blank_theme +
  scale_fill_manual(breaks=cols$Group,
                    values=cols$Colour) +
  ylab("Exon Fraction") +
  facet_grid(.~Fraction)

test <- cover_order_exon %>%
  group_by(Fraction) %>%
  summarise(mean_exon = mean(value))

# Get mean, min, and maximum of intron fraction
mean(cover_order_exon$value)
min(cover_order_exon$value)
max(cover_order_exon$value)

# STAR number of input reads
star.df <- star
rownames(star.df) <- star.df$SampleID

# keep only relevant SampleID in your star table
keep <- rownames(star.df)[rownames(star.df) %in% rownames(md)]
star.df <- star.df[keep,]  # no log2 counts

star.df$Condition <- md$Condition
star.df$Fraction <- md$Fraction

# Arrange in the order of types (HD --> MM)
star.df$Condition <- factor(star.df$Condition, levels=c("HD", "LG", "LV", "MM"))
star_tidy <- star.df %>%
  arrange(Condition)

## Tidying data
star_uni_percent <- star_tidy %>%
  gather(type, value, 2:4) %>%
  filter(type=="%UniqueReads")

star_tidy$SampleID <- factor(star_tidy$SampleID, levels = star_tidy$SampleID)


p4 <- ggplot(star_tidy, aes(x=SampleID, y=InputRead, fill=Condition)) +
  geom_bar(stat = "identity") +
  ylab("Number of Input Reads") +
  scale_fill_manual(breaks=cols$Group,
                    values=cols$Colour) +
  theme(plot.margin = unit(c(0.2,0.2,0.2,0.2),"cm"),
        axis.text.x=element_text(angle = 90, size = 25),
        panel.border = element_rect(colour = "black", fill=NA),
        axis.text.y=element_text(size=35),
        axis.title.y=element_text(size=45),
        axis.title.x=element_text(size=45),
        legend.title=element_blank(),
        legend.text=element_text(size=40),
        panel.grid.major.y = element_line(colour = "black"),
        panel.grid.minor.y = element_line(colour = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

test2 <- star_tidy %>%
  group_by(Fraction) %>%
  summarise(mean_percentunique = mean(`%UniqueReads`))

# Get mean, min, and maxium of input / unique reads
mean(star_tidy$InputRead)
min(star_tidy$InputRead)
max(star_tidy$InputRead)
mean(star_tidy$UniqueRead)
min(star_tidy$UniqueRead)
max(star_tidy$UniqueRead)




p5 <- ggplot(star_uni_percent, aes(x=Condition, y=value, fill=Condition)) +
  geom_violin() +
  geom_boxplot(width=.3, outlier.colour=NA, position= dodge, color="black") +
  ylab("%Unique Read") +
  blank_theme +
  scale_fill_manual(breaks=cols$Group,
                    values=cols$Colour) +
  facet_grid(.~Fraction)


pdf("Figures/Figure_2D_Input_Reads.pdf", width = 40, height = 10)
print(p4)
dev.off()


pdf("Figures/Supp_Figure_5A_percent_unique_reads.pdf", width =9, height = 3)
print(p5)
dev.off()


pdf("Figures/Supp_Figure_5B_Exon_Fraction.pdf", width =9, height = 3)
print(p3)
dev.off()


pdf("Figures/Supp_Figure_5C_Intron_Fraction.pdf", width =9, height = 3)
print(p2)
dev.off()


