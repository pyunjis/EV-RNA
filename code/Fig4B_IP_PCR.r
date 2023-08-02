library(ggplot2)
library(readxl)
library(ggpubr)
library(dplyr)
library(rstatix)

#set working directory
setwd("E:/4. EV-RNA/4. EV-RNA/4. Code/2023-07-23 FINAL/")

#importing data
data <- read_excel("./Files/032723_qPCR_for_stats_comparisons_Figure_4B.xlsx", sheet = 1)
data$Replicate <- as.character(data$Replicate)
data$Pulldown <- factor(data$Pulldown, levels = c("IgG", "ApoB", "CD9"))

#loading themes
{
  blank_theme <- theme_minimal()+
    theme(
      axis.title.x = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA),
      text = element_text(color='black', size=14),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.y=element_text(size=15, colour="black"),
      axis.text.x=element_text(size=14, colour="black"),
      plot.title=element_text(size=14, hjust = 0.5))
  color = "black"
}

#==========================================================================================================
print("------------- Below is the Tukey's from raw CT ------------------")
#==========================================================================================================

  #============================================================
  paste0("manual scaling  plot_ALB")
  #============================================================
  {
    target = "ALB"
    data2 <- data %>%
    filter(Target_name == target)
    # ymax <- max(-(data$CT))
    # ymax <- as.numeric(ymax)
    ymax <- -18.79
    stat.test <- aov(-CT~Pulldown, data=data2) %>%
      tukey_hsd()
    plot_ALB_2 <- ggplot(data2, aes(x = Pulldown, y = -CT)) +
      geom_boxplot(aes(group = Pulldown), fill = "white", color = color, outlier.shape = NA) +
      # geom_point() +
      stat_pvalue_manual(
        stat.test, label= "p.adj.signif",
        y.position = c((ymax + 0.4*ymax), (ymax + 0.3*ymax),  NA)) +
      scale_y_continuous(expand=expansion(mult=c(0.1,0.15))) +
      labs(x = 'sample', y = NULL) +
      ggtitle(target) +
      blank_theme
    plot_ALB_2
  }
  
  paste0("manual plot_B2M")
  #============================================================
  {
    target = "B2M"
    data2 <- data %>%
      filter(Target_name == target)
    # ymax <- max(-(data$CT))
    # ymax <- as.numeric(ymax)
    ymax <- -18.79
    stat.test <- aov(-CT~Pulldown, data=data2) %>%
      tukey_hsd()
    plot_B2M_2 <- ggplot(data2, aes(x = Pulldown, y = -CT)) +
      geom_boxplot(aes(group = Pulldown), fill = "white", color = color) +
      # geom_point() +
      stat_pvalue_manual(
        stat.test, label= "p.adj.signif",
        y.position = c((ymax + 0.2*ymax), (ymax + 0.1*ymax),  NA)) +
      scale_y_continuous(expand=expansion(mult=c(0.1,0.15))) +
      labs(x = 'sample', y = NULL) +
      ggtitle(target) +
      blank_theme
    plot_B2M_2
  }
  
  paste0("manual plot_CORO1C")
  #============================================================
  {
    target = "CORO1C"
    data2 <- data %>%
      filter(Target_name == target)
    # ymax <- max(-(data$CT))
    # ymax <- as.numeric(ymax)
    ymax <- -18.79
    stat.test <- aov(-CT~Pulldown, data=data2) %>%
      tukey_hsd()
    plot_CORO1C_2 <- ggplot(data2, aes(x = Pulldown, y = -CT)) +
      geom_boxplot(aes(group = Pulldown), fill = "white", color = color) +
      # geom_point() +
      stat_pvalue_manual(
        stat.test, label= "p.adj.signif",
        y.position = c((ymax + 0.26*ymax), (ymax + 0.18*ymax),  NA)) +
      scale_y_continuous(expand=expansion(mult=c(0.1,0.15))) +
      labs(x = 'sample', y = NULL) +
      ggtitle(target) +
      blank_theme
    plot_CORO1C_2
  }

  paste0("manual plot_RPS6")
  #============================================================
  {
    target = "RPS6"
    data2 <- data %>%
      filter(Target_name == target)
    # ymax <- max(-(data$CT))
    # ymax <- as.numeric(ymax)
    ymax <- -18.79
    stat.test <- aov(-CT~Pulldown, data=data2) %>%
      tukey_hsd()
    plot_RPS6_2 <- ggplot(data2, aes(x = Pulldown, y = -CT)) +
      geom_boxplot(aes(group = Pulldown), fill = "white", color = color, outlier.shape = NA) +
      # geom_point() +
      stat_pvalue_manual(
        stat.test, label= "p.adj.signif",
        y.position = c((ymax + 0.32*ymax), (ymax + 0.18*ymax),  NA)) +
      scale_y_continuous(expand=expansion(mult=c(0.1,0.15))) +
      labs(x = 'sample', y = NULL) +
      ggtitle(target) +
      blank_theme
    plot_RPS6_2
  }
  
#putting all 4 into 1 plot
Figure_4B <- ggarrange(plot_ALB_2, plot_B2M_2, plot_CORO1C_2, plot_RPS6_2, 
            ncol = 2, nrow = 2,  
            common.legend = TRUE)
Figure_4B

#Output .pdf
pdf("./Figures/Figure_4B_IP_PCR.pdf", width =4.465, height = 4.7)
print(Figure_4B)
dev.off()  
  
  
  
  
  