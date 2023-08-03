#Figure 1B

#set WD
setwd("E:/4. EV-RNA/4. EV-RNA/4. Code/2023-07-23 FINAL/")
getwd()

#load libraries:
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

#manual organization to ensure proper tidy data for R visualization
Qnano_tidy <- read_excel("./Files/Fig1B_qNano_A280.xlsx", sheet = 1)


#generating background for fraction ranges
data_breaks <- data.frame(start = c(0, 4, 8, 15, 22, 29),  # Create data with breaks
                          end = c(4, 8, 12, 19, 26, 33),
                          FR_ranges = c("FR14", "FR58", "FR912", "FR1619", "FR2326", "FR3033"))
data_breaks$FR_ranges <- factor(data_breaks$FR_ranges, levels = c("FR14", "FR58", "FR912", "FR1619", "FR2326", "FR3033"))

#visual scaling factor to plot both graphs on 1 plot
scaleFactor <- 89000000000 / max(Qnano_tidy$A280)

#calculating mean and sd for error_bar visualization
Qnano_tidy_stats <- Qnano_tidy %>% group_by(Fraction) %>% mutate(avg_a280 = mean(A280), sd_a280 = sd(A280), avg_concentration = mean(Concentration), sd_concentration = sd(Concentration)) %>% ungroup()

#plotting all together with bar, points, and error_bars
Figure_1B <- ggplot() +
  #adding background range colors
  geom_rect(data = data_breaks,
            aes(xmin = start,
                xmax = end,
                ymin = - Inf,
                ymax = Inf,
                fill = FR_ranges),
            alpha = 0.8) +
  scale_fill_manual(values = c("#ffabab", "#a7c6e0", "#b6deb4", "#d5b6d9", "#ffca96", "#ffffad")) +
  #both values of data
  geom_errorbar(data = Qnano_tidy_stats, aes(x=Fraction-0.25, ymin=avg_concentration-sd_concentration, ymax=avg_concentration+sd_concentration), width=.2, color="black") +
  geom_bar(data = Qnano_tidy, aes(x=Fraction-0.25, y=Concentration/3), stat="identity", width=0.5, color="#616161", fill="#616161") +
  geom_point(data = Qnano_tidy, aes(x=Fraction-0.25, y=Concentration), color="black") +
 
  geom_errorbar(data = Qnano_tidy_stats, aes(x=Fraction+0.25, ymin=(avg_a280-sd_a280)*scaleFactor, ymax=(avg_a280+sd_a280)*scaleFactor), width=.2, color = "red") +
  geom_bar(data = Qnano_tidy, aes(x=Fraction+0.25, y=A280/3 * scaleFactor), stat="identity", width=0.5, color="red", fill="red") +
  geom_point(data = Qnano_tidy, aes(x=Fraction+0.25, y=A280 * scaleFactor), color="darkred") +
  #adding second axis and themeing
  scale_x_continuous(breaks=seq(0,40,by=5), expand = c(0, 0)) +
  xlab("Fraction") +
  scale_y_continuous(name="Concentration (particles/ml)", sec.axis=sec_axis(~./scaleFactor, name="Plasma Protein (A280)"), breaks=seq(0,120000000000,by=20000000000), limits = c(0, 120000000000), expand = c(0, 0)) +
  theme_classic() +
  theme(
    legend.position = "none",
    panel.border = element_rect(linewidth = 1.5, fill = NA),
    axis.title.x=element_text(color="black", face = "bold"),
    axis.text.x=element_text(color="black", face = "bold"),
    axis.title.y.left=element_text(color="black", face = "bold"),
    axis.text.y.left=element_text(color="black", face = "bold"),
    axis.title.y.right=element_text(color="black", face = "bold"),
    axis.text.y.right=element_text(color="black", face = "bold"))

Figure_1B

#Output .pdf
pdf("./Figures/Figure_1B_qNano_Abs.pdf", width =5, height = 2.5)
print(Figure_1B)
dev.off()

