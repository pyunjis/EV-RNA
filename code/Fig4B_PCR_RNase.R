# Figure 4B_PCR_RNase

#loading library
library(tidyverse)
library(readxl)
library(janitor)
library(skimr)
library(GGally)
library(broom)
library(moderndive)
library(here)
library(ggplot2)
library(pheatmap)
library(dplyr)

setwd("D:/Box drive update/4. EV-RNA/4. Code/2022-02-17 V7/")

# cts <- read_xlsx(("2019-07-12 124918 Triton and RNase treated plasma.xlsx"), sheet = 6) %>%
#   clean_names()


cts <- read_xlsx(("Files/2020-09-10 161238 RNase_3.xlsx"), sheet = 2) %>%
  clean_names()

# Make into average
cts2 <- cts %>%
  group_by(sample_name, target_name) %>%
  dplyr::mutate(avg_ct = mean(ctcorr)) %>%
  dplyr::mutate(sd_ct = sd(ctcorr))

# ggplot ERCC corrected cts
# Filter for desired genes
pcr <- cts2 %>% filter(target_name %in% c("ALB","APOE", "B2M", "CORO1C", "CPOX", "RPS6")) %>%  
  ggplot(., aes(x=sample_name, y=-avg_ct, color=target_name, group=target_name)) +
  geom_point(size=2, color = "black") +
  geom_line(size=1, color = "black") +
  geom_errorbar(aes(ymin = -avg_ct - sd_ct, ymax = -avg_ct + sd_ct), color = "black", width=1, position=position_dodge(0.05)) +
  facet_wrap(target_name ~ ., scales= "free_y") +
  ylab("raw -ct") +
  #ggtitle("RNase and Triton Experiment") +
  theme(panel.grid.major = element_line(colour = "white", size = 1)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  #theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  #theme(strip.background = element_rect(colour = "grey50", fill = "white")) + 
  #theme(text = element_text(size=12), axis.text.x = element_text(angle=90, hjust=1)) +
  theme_bw() +
  theme(axis.title.x=element_blank(),
  axis.text.x=element_blank()) # removing  x-axis labels

pdf("Figures/Figure_4B_PCR_RNase.pdf", width =5, height = 3)
print(pcr)
dev.off()

