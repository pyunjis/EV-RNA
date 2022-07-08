# Figure 1D pie chart for TEM
# Supp Fig 1 histogram for TEM

#Load libraries
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(gplots)
library(RColorBrewer)
library(stringr)
library(reshape2)
library(dplyr)
library(scales)
library(readxl)
library(scales)
library(gg.gap) #for cutting ggplot axes into segment ranges

# Set working directory
setwd("E:/Box drive update/4. EV-RNA/4. Code/2022-02-17 V7")

# Import data
raw_tem <- read_excel("./Files/EV-RNA-TEM_Distributions3.xlsx", sheet = 2)
raw_tem_step5 <- read_excel("./Files/EV-RNA-TEM_Distributions3.xlsx", sheet = 3)

# Tidy data
hist_dat <- raw_tem %>% gather(Fraction, Size, 1:5) %>% na.omit()
hist_dat$Fraction <- factor(hist_dat$Fraction, levels = c("FR4", "FR6", "FR8", "FR10", "FR12"))

hist_dat5 <- raw_tem_step5 %>% gather(Fraction, Size, 1:5) %>% na.omit()
hist_dat5$Fraction <- factor(hist_dat5$Fraction, levels = c("FR4", "FR6", "FR8", "FR10", "FR12"))
hist_dat5 <- hist_dat5 %>% filter(Size != "200")

#combining top and bottom into 1 plot by segmenting y axis
segm <- hist_dat5 %>% 
  #filter(Fraction == "FR12") %>%
  ggplot(., aes(x=Size, fill=Fraction, color = Fraction)) +
  geom_histogram(binwidth=5, position="identity", alpha = 0.4) +
  facet_grid(Fraction~.) +
  xlim(c(20,400)) +
  theme(legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  theme_classic()
segm

# Normalized histogram
# Extract values of observations for normalization
pg <- ggplot_build(segm)
head(pg)
unnorm_y <- pg$data[[1]]$y # these are the histogram heights for each bar
unnorm_x <- pg$data[[1]]$x # these are the x-values for each histogram bar
unnorm_x.df <- as.data.frame(unnorm_x)
unnorm_y.df <- as.data.frame(unnorm_y)
Fraction <- c(rep("FR4", 77), rep("FR6", 77), rep("FR8", 77), rep("FR10", 77), rep("FR12", 77))
unnorm <- cbind(unnorm_x.df, unnorm_y.df, Fraction)


# Normalized histogram
norm5 <- unnorm
norm5 <- norm5 %>% group_by(Fraction) %>%
  mutate(Norm = (unnorm_y-min(unnorm_y))/(max(unnorm_y)-min(unnorm_y))*100)
norm5$Size <- norm5$unnorm_x
norm5$Count <- norm5$unnorm_y
norm5$Fraction <- factor(norm5$Fraction, levels = c("FR4", "FR6", "FR8", "FR10", "FR12"))

# Plotting normalized histogram
hist <- norm5 %>% 
  #filter(Fraction == "FR12") %>%
  ggplot(., aes(x=Size, y=Norm, fill=Fraction, color=Fraction)) +
  geom_bar(size=0.5, stat="identity", alpha =0.4, fill="#a9aaac", color = "#959698") +
  facet_grid(Fraction~.) +
  xlim(c(20,400)) +
  ylab("Normalized particle counts") +
  theme(legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))+ 
  #geom_vline(xintercept = 47.5) +
  #geom_vline(xintercept = 97.5) +
  theme_classic()
hist

pdf(paste0("./Figures/Supp_Figure_1_TEM_histogram2.pdf"), width = 5, height = 3.5)
print(hist)
dev.off()

##-------------------------- converting histogram by 5 into pie charts --------------------------
norm5 <- norm5 %>% select(-unnorm_x, -unnorm_y)
norm5 <- norm5 %>% mutate(Group = ifelse(Size < 50, "<50", 
                                         ifelse(Size > 45 & Size < 100, "50-100", "\u2265 100")))
norm5_pie <- norm5 %>%
  group_by(Group, Fraction) %>%
  mutate(value = sum(Count))
#FR4
df_FR4_g1 <- norm5_pie %>%
  filter(Fraction == "FR4", Group == "<50") %>% slice(1)
df_FR4_g2 <- norm5_pie %>%
  filter(Fraction == "FR4", Group == "50-100") %>% slice(1)
df_FR4_g3 <- norm5_pie %>%
  filter(Fraction == "FR4", Group == "\u2265 100") %>% slice(1)
df_FR4 <- rbind(df_FR4_g1, df_FR4_g2, df_FR4_g3)
#FR6
df_FR6_g1 <- norm5_pie %>%
  filter(Fraction == "FR6", Group == "<50") %>% slice(1)
df_FR6_g2 <- norm5_pie %>%
  filter(Fraction == "FR6", Group == "50-100") %>% slice(1)
df_FR6_g3 <- norm5_pie %>%
  filter(Fraction == "FR6", Group == "\u2265 100") %>% slice(1)
df_FR6 <- rbind(df_FR6_g1, df_FR6_g2, df_FR6_g3)
#FR8
df_FR8_g1 <- norm5_pie %>%
  filter(Fraction == "FR8", Group == "<50") %>% slice(1)
df_FR8_g2 <- norm5_pie %>%
  filter(Fraction == "FR8", Group == "50-100") %>% slice(1)
df_FR8_g3 <- norm5_pie %>%
  filter(Fraction == "FR8", Group == "\u2265 100") %>% slice(1)
df_FR8 <- rbind(df_FR8_g1, df_FR8_g2, df_FR8_g3)
#FR10
df_FR10_g1 <- norm5_pie %>%
  filter(Fraction == "FR10", Group == "<50") %>% slice(1)
df_FR10_g2 <- norm5_pie %>%
  filter(Fraction == "FR10", Group == "50-100") %>% slice(1)
df_FR10_g3 <- norm5_pie %>%
  filter(Fraction == "FR10", Group == "\u2265 100") %>% slice(1)
df_FR10 <- rbind(df_FR10_g1, df_FR10_g2, df_FR10_g3)
#FR12
df_FR12_g1 <- norm5_pie %>%
  filter(Fraction == "FR12", Group == "<50") %>% slice(1)
df_FR12_g2 <- norm5_pie %>%
  filter(Fraction == "FR12", Group == "50-100") %>% slice(1)
df_FR12_g3 <- norm5_pie %>%
  filter(Fraction == "FR12", Group == "\u2265 100") %>% slice(1)
df_FR12 <- rbind(df_FR12_g1, df_FR12_g2, df_FR12_g3)

df_FR <- rbind(df_FR4, df_FR6, df_FR8, df_FR10, df_FR12)

#calculating percentages
df_FR_per <- df_FR %>% group_by(Fraction) %>%
  mutate(per_value = round((value/sum(value))*100, digits = 1)) %>%
  mutate(per_val_text = value/sum(value))

#Function for modifying percent values displayed in pie charts:
french_percent <- label_percent(
  accuracy = 0.1,
  decimal.mark = ".",
  suffix = " %"
)

df_FR_per$per_val_text <- french_percent(df_FR_per$per_val_text)

##-------------------------- FR4 --------------------------
pie_colors <- c("#4daf4a", "#e41a1c", "#377eb8" )
bp_FR4 <- df_FR_per %>% filter(Fraction == "FR4") %>% 
  ggplot(., aes(x="", y=value, fill=Group))+
  geom_bar(width = 1, stat = "identity", colour="white", size=1.25#adjust size for outline thickness
  )
pie_FR4 <- bp_FR4 + coord_polar("y", start=3.14)
pie_FR4_plot <- pie_FR4 +  scale_fill_manual(values = pie_colors) + theme_minimal() +
  theme(axis.text.x=element_blank()) +
  geom_text(aes(y = value*0.5 + c(30#modify for the desired position (note this is scaled to raw data)
                                  , cumsum(value*0.8)[-length(value)]), 
                label = per_val_text), size=10, colour="black")+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "white"),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.position = "bottom"
      
  )
pie_FR4_plot

#saving plots FR4
# pdf(paste0("./FR4_TEM_pie_chart5.pdf"))
# print(pie_FR4_plot)
# dev.off()

##-------------------------- FR6 --------------------------
pie_colors <- c("#4daf4a", "#e41a1c", "#377eb8" )
bp_FR6 <- df_FR_per %>% filter(Fraction == "FR6") %>% 
  ggplot(., aes(x="", y=value, fill=Group))+
  geom_bar(width = 1, stat = "identity", colour="white", size=1.25#adjust size for outline thickness
  )
pie_FR6 <- bp_FR6 + coord_polar("y", start=3.14)
pie_FR6_plot <- pie_FR6 +  scale_fill_manual(values = pie_colors) + theme_minimal() +
  theme(axis.text.x=element_blank()) +
  geom_text(aes(y = value*0.4 + c(22#modify for the desired position (note this is scaled to raw data)
                                  , cumsum(value*0.8)[-length(value)]), 
                label = per_val_text), size=10, colour="black")+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "white"),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.position = "bottom"
  )
pie_FR6_plot
#saving plots FR6
# pdf(paste0("./FR6_TEM_pie_chart5.pdf"))
# print(pie_FR6_plot)
# dev.off()


##-------------------------- FR8 --------------------------
pie_colors <- c("#4daf4a", "#377eb8")
bp_FR8 <- df_FR_per %>% filter(Fraction == "FR8") %>% 
  filter(Group != "\u2265 100") %>%
  ggplot(., aes(x="", y=value, fill=Group))+
  geom_bar(width = 1, stat = "identity", colour="white", size=1.25#adjust size for outline thickness
  )
pie_FR8 <- bp_FR8 + coord_polar("y", start=3.14)
pie_FR8_plot <- pie_FR8 +  scale_fill_manual(values = pie_colors) + theme_minimal() +
  theme(axis.text.x=element_blank()) +
  geom_text(aes(y = value*0.45 + c(15#modify for the desired position (note this is scaled to raw data)
                                   , cumsum(value*1.1)[-length(value)]), 
                label = per_val_text), size=10, colour="black")+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "white"),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.position = "bottom"
  )
pie_FR8_plot
#saving plots FR8
# pdf(paste0("./FR8_TEM_pie_chart5.pdf"))
# print(pie_FR8_plot)
# dev.off()


##-------------------------- FR10 --------------------------
pie_colors <- c("#4daf4a", "#377eb8")
bp_FR10 <- df_FR_per %>% filter(Fraction == "FR10") %>% 
  filter(Group != "\u2265 100") %>%
  ggplot(., aes(x="", y=value, fill=Group))+
  geom_bar(width = 1, stat = "identity", colour="white", size=1.25#adjust size for outline thickness
  )
pie_FR10 <- bp_FR10 + coord_polar("y", start=3.14)
pie_FR10_plot <- pie_FR10 +  scale_fill_manual(values = pie_colors) + theme_minimal() +
  theme(axis.text.x=element_blank()) +
  geom_text(aes(y = value*0.15 + c(300#modify for the desired position (note this is scaled to raw data)
                                   , cumsum(value*1.02)[-length(value)]), 
                label = per_val_text), size=10, colour="black")+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "white"),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.position="bottom"
  )
pie_FR10_plot
#saving plots FR10
# pdf(paste0("./FR10_TEM_pie_chart5.pdf"))
# print(pie_FR10_plot)
# dev.off()


##-------------------------- FR12 --------------------------
pie_colors <- c("#4daf4a", "#377eb8")
bp_FR12 <- df_FR_per %>% filter(Fraction == "FR12") %>% 
  filter(Group != "\u2265 100") %>%
  ggplot(., aes(x="", y=value, fill=Group))+
  geom_bar(width = 1, stat = "identity", colour="white", size=1.25#adjust size for outline thickness
  )
pie_FR12 <- bp_FR12 + coord_polar("y", start=3.14)
pie_FR12_plot <- pie_FR12 +  scale_fill_manual(values = pie_colors) + theme_minimal() +
  theme(axis.text.x=element_blank()) +
  geom_text(aes(y = value*0.1 + c(160#modify for the desired position (note this is scaled to raw data)
                                  , cumsum(value*1.0)[-length(value)]), 
                label = per_val_text), size=10, colour="black")+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "white"),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.position="bottom"
  )
pie_FR12_plot
#saving plots FR12
# pdf(paste0("./FR12_TEM_pie_chart5.pdf"))
# print(pie_FR12_plot)
# dev.off()

# Save all pie charts into a figure
pdf(("./Figures/Figure_1D_TEM_pie_chart.pdf"), 30, 5, useDingbats = FALSE)
cowplot::plot_grid(pie_FR4_plot, pie_FR6_plot, pie_FR8_plot, pie_FR10_plot, pie_FR12_plot, nrow = 1,  labels = c('FR4', 'FR6', 'FR8', 'FR10', 'FR12'), label_size = 25)
dev.off()
