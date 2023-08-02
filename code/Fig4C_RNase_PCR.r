# Figure_4B_Triton_RNase_Tukeys Revised

#loading library
{
library(tidyverse)
library(readxl)
library(janitor)
library(skimr)
library(GGally)
library(broom)
library(moderndive)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(ggpubr)
library(rstatix)
library(lemon)
}

#set WD for wherever is needed
setwd("E:/4. EV-RNA/4. EV-RNA/4. Code/2023-07-23 FINAL/")
getwd()

#import data
cts <- read_xlsx(("Files/2020-09-10 161238 RNase.xlsx"), sheet = 1) %>%
  clean_names()

{
# Make into average
cts2 <- cts %>%
  group_by(sample_name, target_name) %>%
  dplyr::mutate(avg_ct = mean(ct)) %>%
  dplyr::mutate(sd_ct = sd(ct))
}

#All 4 genes plotted at once using facet
{
data <- cts2 %>%
  #filter(sample == "S1" | sample == "S2FR" | sample == "S1FRS2") %>%
  filter(target_name %in% c("ALB", "B2M", "CORO1C", "RPS6"))
plot_colors <- c("black", "black")
color <- "white"
Figure_4C <- ggplot(data, aes(x = sample_name, y = -ct, group = sample, color = sample)) +
  geom_boxplot(aes(group = sample_name), linewidth =0.5, fill = color) +
  geom_line(aes(y = -avg_ct, group = sample, color = sample), alpha =0.5) +
  # geom_point(data = data, aes(x=sample_name, y=-ctcorr, color = sample), size = 1) +
  scale_color_manual(values = plot_colors) +
  facet_rep_grid(target_name ~ ., scales= "free_y") +
  ylab("raw -ct") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "white", size = 1)) +
  theme(panel.grid.minor = element_line(colour = "white", size = 1)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  labs(x = NULL, y = NULL) +
theme(axis.text.y = element_text(size = 12),
      axis.ticks.x = element_line(linewidth = 0.75),
      axis.ticks.y = element_line(linewidth = 0.75),
      panel.border = element_rect(linewidth = 1.25, fill = NA))

}



pdf("Figures/Figure_4C_Triton_RNase.pdf", width =6, height = 12)
print(Figure_4C)
dev.off()



