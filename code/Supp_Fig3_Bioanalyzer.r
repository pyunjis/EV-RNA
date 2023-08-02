# Supplementary Figure 3a

library(vctrs)
library(cli)
library(tidyverse)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(gplots)
library(RColorBrewer)
library(stringr)
library(readxl)
library(scales)

#set WD
setwd("E:/4. EV-RNA/4. EV-RNA/4. Code/2023-07-23 FINAL/")
getwd()
#import and tidy bioanalyzer trace
BAtrace <- read_excel("./Files/BA traces.xlsx", sheet = 1)
BAtrace$nrow <- rownames(BAtrace)
BAtrace_tidy <- BAtrace %>% gather(type, value, 4:5)


#splitting type into 2 dataframes and finding smooth line between
{
BAtrace_cfrna <- BAtrace_tidy %>% filter(type == "cfRNA")
spline_int_cfrna <- as.data.frame(spline(BAtrace_cfrna$nrow, BAtrace_cfrna$value))
spline_int_cfrna$type <- BAtrace_cfrna$type
BAtrace_tissue <- BAtrace_tidy %>% filter(type == "Tissue RNA")
spline_int_tissue <- as.data.frame(spline(BAtrace_tissue$nrow, BAtrace_tissue$value))
spline_int_tissue$type <- BAtrace_tissue$type
#recombine both traces
spline_int <- rbind(spline_int_cfrna, spline_int_tissue)
}

#processing to get reference nucleotide values and positions
{
bp_marks <- unique(BAtrace$bp)
bp_marks <- bp_marks[-1]
bp_ref <- BAtrace %>% filter(bp %in% c(bp_marks))
bp_ref$nrow <- as.numeric(bp_ref$nrow)
bp_label <- bp_ref$bp
bp_break <- bp_ref$nrow
maxlim <- as.numeric(length(BAtrace$nrow))
}

#plot with axis labeled by known nucleotide elution positions from BAtrace$bp
Figure_S3 <- ggplot(BAtrace_tidy, aes(x=nrow, y= value, group = type)) +
  geom_line(data = spline_int, aes(x = x, y = y, linetype=type)) +
  ylab("FU") +
  xlab("Fragment Size (nt)") +
  scale_y_continuous(breaks = c(-0.5,0,0.5,1,1.5,2,2.5,3), limits = c(-0.5,3), expand = c(0, 0)) +
  scale_x_continuous(labels = c(bp_label), 
                     breaks = c(bp_break), 
                     limits = c(1,maxlim), expand = c(0, 0)) +
  theme_classic()

Figure_S3

#Output .pdf
pdf("./Figures/Supp_Figure_3A_Bioanalyzer.pdf", width =10, height = 5)
print(Figure_S3)
dev.off()
