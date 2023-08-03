#Supplementary Figure S2

#set WD
setwd("E:/4. EV-RNA/4. EV-RNA/4. Code/2023-07-23 FINAL/")

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
library(scales)


#import data:
raw_tem <- read_excel("./Files/Fig1D_and_FigS2_TEM_data.xlsx", sheet = 1)

#turning into tidy data:
hist_dat <- raw_tem %>% gather(Fraction, Size, 1:5) %>% na.omit()
hist_dat$Fraction <- factor(hist_dat$Fraction, levels = c("FR4", "FR6", "FR8", "FR10", "FR12"))

#calculating the total in each size range
pie_dat <- hist_dat
pie_dat <- pie_dat %>% 
  mutate(Group = ifelse(Size < 50, "< 50 nm",
                      ifelse(Size > 100, "> 100 nm", "> 50 nm & < 100 nm")))
pie_dat$Group <- factor(pie_dat$Group, levels = c("< 50 nm", "> 50 nm & < 100 nm", "> 100 nm"))
pie_dat <- pie_dat %>% 
  group_by(Fraction, Group) %>%
  mutate(value = n()) %>% dplyr::slice(1) 

french_percent <- label_percent(
  accuracy = 1,
  decimal.mark = ".",
  suffix = " %"
)

#adding labels for the size range percentages
pie_FR_per <- pie_dat %>% group_by(Fraction) %>%
  mutate(per_value = round((value/sum(value))*100, digits = 1)) %>%
  mutate(per_val_text = value/sum(value))
pie_FR_per$per_val_text <- french_percent(pie_FR_per$per_val_text)
df_FR_per <- pie_FR_per

#generating formatting theme
theme_pie <- theme(axis.title.y=element_blank(),
                   axis.text.y=element_blank(),
                   axis.ticks.y=element_blank(),
                   axis.title.x=element_blank(),
                   axis.text.x=element_blank(),
                   axis.ticks.x=element_blank(),
                   panel.border = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   axis.line = element_line(colour = "white"))

pie_font_size <- 7
{
##-------------------------- FR4 --------------------------
pie_colors <- c("#02ba3b", "#0d7fbf", "#ea0000")

bp_FR4 <- df_FR_per %>% filter(Fraction == "FR4") %>% 
  ggplot(., aes(x="", y=value, fill=Group))+
  geom_bar(width = 1, stat = "identity", colour="white", size=1.25#adjust size for outline thickness
  )
bp_FR4
pie_FR4 <- bp_FR4 + coord_polar("y", start=0)
pie_FR4_plot <- pie_FR4 +  scale_fill_manual(values = pie_colors) + theme_minimal() +
  theme(axis.text.x=element_blank()) +
    geom_text(aes(y = value + c(1,+6,-5), #modify c(y1,y2,y3) for the desired label position (note this is scaled to raw data)
                  label = per_val_text), size=pie_font_size, colour="black")+
  theme_pie +
  theme(legend.position="none")
pie_FR4_plot

##-------------------------- FR6 --------------------------
pie_colors <- c("#02ba3b", "#0d7fbf")
bp_FR6 <- df_FR_per %>% filter(Fraction == "FR6") %>% 
  ggplot(., aes(x="", y=value, fill=Group))+
  geom_bar(width = 1, stat = "identity", colour="white", size=1.25#adjust size for outline thickness
  )
bp_FR6
pie_FR6 <- bp_FR6 + coord_polar("y", start=0)
pie_FR6_plot <- pie_FR6 +  scale_fill_manual(values = pie_colors) + theme_minimal() +
  theme(axis.text.x=element_blank()) +
  geom_text(aes(y = value + c(23,-15), #modify c(y1,y2,y3) for the desired label position (note this is scaled to raw data)
                label = per_val_text), size=pie_font_size, colour="black")+
  theme_pie +
theme(legend.position="none")
pie_FR6_plot

##-------------------------- FR8 --------------------------
pie_colors <- c("#02ba3b", "#0d7fbf")
bp_FR8 <- df_FR_per %>% filter(Fraction == "FR8") %>% 
  ggplot(., aes(x="", y=value, fill=Group))+
  geom_bar(width = 1, stat = "identity", colour="white", size=1.25#adjust size for outline thickness
  )
bp_FR8
pie_FR8 <- bp_FR8 + coord_polar("y", start=0)
pie_FR8_plot <- pie_FR8 +  scale_fill_manual(values = pie_colors) + theme_minimal() +
  theme(axis.text.x=element_blank()) +
  geom_text(aes(y = value + c(-65,-18), #modify c(y1,y2,y3) for the desired label position (note this is scaled to raw data)
                label = per_val_text), size=pie_font_size, colour="black")+
  theme_pie +
  theme(legend.position="none")
pie_FR8_plot

##-------------------------- FR10 --------------------------
pie_colors <- c("#02ba3b", "#0d7fbf")
bp_FR10 <- df_FR_per %>% filter(Fraction == "FR10") %>% 
  ggplot(., aes(x="", y=value, fill=Group))+
  geom_bar(width = 1, stat = "identity", colour="white", size=1.25#adjust size for outline thickness
  )
bp_FR10
pie_FR10 <- bp_FR10 + coord_polar("y", start=0)
pie_FR10_plot <- pie_FR10 +  scale_fill_manual(values = pie_colors) + theme_minimal() +
  theme(axis.text.x=element_blank()) +
  geom_text(aes(y = value + c(-550,100), #modify c(y1,y2,y3) for the desired label position (note this is scaled to raw data)
                label = per_val_text), size=pie_font_size, colour="black")+
  theme_pie +
  theme(legend.position="none")
pie_FR10_plot

##-------------------------- FR12 --------------------------
pie_colors <- c("#02ba3b", "#0d7fbf")
bp_FR12 <- df_FR_per %>% filter(Fraction == "FR12") %>% 
  ggplot(., aes(x="", y=value, fill=Group))+
  geom_bar(width = 1, stat = "identity", colour="white", size=1.25#adjust size for outline thickness
  )
bp_FR12
pie_FR12 <- bp_FR12 + coord_polar("y", start=0)
pie_FR12_plot <- pie_FR12 +  scale_fill_manual(values = pie_colors) + theme_minimal() +
  theme(axis.text.x=element_blank()) +
  geom_text(aes(y = value + c(-900,140), #modify c(y1,y2,y3) for the desired label position (note this is scaled to raw data)
                label = per_val_text), size=pie_font_size, colour="black")+
  theme_pie +
  theme(legend.position="none")
pie_FR12_plot

#making the legend
df_FR_per_legend <- df_FR_per
df_FR_per_legend$Group <- factor(df_FR_per_legend$Group, levels = c("> 100 nm", "> 50 nm & < 100 nm", "< 50 nm"))
legend <- ggplot(df_FR_per_legend, aes(x=Fraction, y= per_value, color = Group, fill = Group)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values =  c("#ea0000", "#0d7fbf", "#02ba3b")) +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank()) 
legend <- cowplot::get_legend(legend)
as_ggplot(legend)

#putting all plots together
ggarrange(pie_FR4_plot, pie_FR6_plot, pie_FR8_plot, pie_FR10_plot, pie_FR12_plot, NA, NA, legend, NA, NA, ncol = 5, nrow = 2)

Figure_1D <- ggarrange(pie_FR4_plot, pie_FR6_plot, pie_FR8_plot, pie_FR10_plot, pie_FR12_plot, ncol = 5, nrow = 1)
Figure_1D
}

#Output .pdf
pdf("./Figures/Figure_1D_TEM.pdf", width =15, height = 4.5)
print(Figure_1D)
dev.off()



#=============================================================
print("Plotting particle size distribution as a histogram for Supplemental Figure S2")

#plotting unnormalized particle size counts together
segm <- hist_dat %>% 
  ggplot(., aes(x=Size, fill=Fraction, color = Fraction)) +
  geom_histogram(binwidth=5, position="identity", alpha = 0.4) +
  facet_grid(Fraction~.) +
  xlim(c(20,400)) +
  theme(legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  theme_classic()
segm

##-----------<<<<<<<<< Normalized Histogram >>>>>>>>>>>---------------
#Using unnormalized plot values to generate normalized
pg <- ggplot_build(segm)
head(pg)
unnorm_y <- pg$data[[1]]$y # these are the histogram heights for each bar
unnorm_x <- pg$data[[1]]$x # these are the x-values for each histogram bar
unnorm_x.df <- as.data.frame(unnorm_x)
unnorm_y.df <- as.data.frame(unnorm_y)
Fraction <- c(rep("FR4", 77), rep("FR6", 77), rep("FR8", 77), rep("FR10", 77), rep("FR12", 77))
unnorm <- cbind(unnorm_x.df, unnorm_y.df, Fraction)


#normalized histogram:
norm <- unnorm
norm <- norm %>% group_by(Fraction) %>%
  mutate(Norm = (unnorm_y-min(unnorm_y))/(max(unnorm_y)-min(unnorm_y))*100)
norm$Size <- norm$unnorm_x
norm$Count <- norm$unnorm_y
norm$Fraction <- factor(norm$Fraction, levels = c("FR4", "FR6", "FR8", "FR10", "FR12"))


#plotting normalized histogram
seg <- norm %>% 
  ggplot(., aes(x=Size, y=Norm, fill=Fraction, color=Fraction)) +
  geom_bar(stat="identity", alpha =0.4, color = "darkgrey", fill = "darkgrey") +
  facet_grid(Fraction~.) +
  xlim(c(15,400)) +
  ylab("Normalized Particle Counts") +
  xlab("Size Distribution (nm)") +
  geom_vline(xintercept = 47.5) +
  geom_vline(xintercept = 97.5) +
  theme(legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  theme_classic()
seg

pdf(paste0("./Figures/Supp_Figure_2_TEM.pdf"))
print(seg)
dev.off()
