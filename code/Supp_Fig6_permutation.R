library(dplyr)
library(data.table)
library(purrr)
library(readr)
library(stringr)
library(grid)
library(gridExtra)
library(ggplot2)

#setwd("C:/Users/doe/projects/new_norm/permutation/")
setwd("E:/Box drive update/4. EV-RNA/4. Code/2022-02-17 V7/Files/updated_output/")
read_plus <- function(flnm) {
  read_csv(flnm) %>% 
    mutate(filename = flnm)
}

tbl <- 
  list.files(path = "updated_output/", "permutation.list", full.names = T) %>%
  map_df(~read_plus(.))

diff.genes_files <- list.files(path = "updated_output/", "diff.genes.csv", full.names = T)

permutation_files <- list.files(path = "updated_output/", "permutation.list", full.names = T)

#tab1 <- read_plus(diff.genes_files[1])

#baseline <- paste(str_split(diff.genes_files[1], "_")[[1]][7], str_split(diff.genes_files[1], "_")[[1]][8], sep = "_")
#target <- paste(str_split(diff.genes_files[1], "_")[[1]][7], str_split(diff.genes_files[1], "_")[[1]][6], sep = "_")

#p <- ggplot(subset(tab1, NumDiffGenes < 1000), aes(x=NumDiffGenes)) +
  #xlim(c(0,500))+
#  geom_histogram(bins=100, fill = bar_color) +
#  geom_vline(data=tab1, mapping=aes(xintercept = Actual, color = paste0("Correct Label: ", Actual, " DE genes")),
#             linetype="longdash", size=0.6, show.legend = T) +
#  scale_color_manual(values = "blue", name = "Number of DE genes") +
  #geom_text(aes(x=Actual, label=paste0(Actual, "DE genes"), y=0), colour="blue", angle=90, size=5) +
  #geom_text(aes(x=Actual+80, label = paste0(Actual, "DE genes"), y = 500), color = "blue", size = 5) +
#  ggtitle(paste(1000, "Random Permutations:", baseline, "vs", target)) +#, subtitle = paste0("Number of permutations with more than 1000 DE genes:", nrow(subset(tab1, NumDiffGenes > 1000)))) +
  #labs(caption = paste0("Percentage of permutations with more DE genes than actual:", nrow(subset(tab1, NumDiffGenes < Actual))/nrow(tab1)*100,"%")) +
#  xlab("Number of significant genes") +
#  theme_bw() +
#  theme(plot.margin = unit(c(1,1,2,1), "cm"),
#        aspect.ratio=1,
#        plot.title = element_text(hjust = 0.5),
#        legend.title = element_text(size=10, hjust = 0.5),
#        legend.position = "none") +
#  annotation_custom(grob = textGrob(label = tab1$Actual[1], hjust = .5, gp = gpar(cex = 1.2, col = "blue")),
#                    ymin = -6, ymax = -6, xmin = tab1$Actual[1], xmax = tab1$Actual[1]) +
#  coord_cartesian(clip = 'off')

#grob <- grobTree(textGrob(paste0(nrow(subset(tab1, NumDiffGenes < Actual))/nrow(tab1)*100,"%"), x=0.1,  y=0.95, hjust=0,
#                          gp=gpar(col="red", fontsize=15, fontface="italic")))

#p + annotation_custom(grob)

colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33")
diff.genes_files <- diff.genes_files[c(1,7,13,5,11,17,6,12,18,2,8,14,3,9,15,4,10,16)]
p_list <- list()
for (i in 1:length(diff.genes_files)) {
  tab1 <- read_plus(diff.genes_files[i])
  
  baseline <- paste(str_split(diff.genes_files[i], "_")[[1]][7], str_split(diff.genes_files[i], "_")[[1]][8], sep = "_")
  target <- paste(str_split(diff.genes_files[i], "_")[[1]][7], str_split(diff.genes_files[i], "_")[[1]][6], sep = "_")
  
  if (grepl("FR14", baseline)) {
    bar_color <- colors[1]
  }
  if (grepl("FR58", baseline)) {
    bar_color <- colors[2]
  }
  if (grepl("FR912", baseline)) {
    bar_color <- colors[3]
  }
  if (grepl("FR1619", baseline)) {
    bar_color <- colors[4]
  }
  if (grepl("FR2326", baseline)) {
    bar_color <- colors[5]
  }
  if (grepl("FR3033", baseline)) {
    bar_color <- colors[6]
  }
  
  p <- ggplot(subset(tab1, NumDiffGenes < 200), aes(x=NumDiffGenes)) +
    #xlim(c(0,200))+
    geom_histogram(bins=100, fill = bar_color) +
    geom_vline(data=tab1, mapping=aes(xintercept = Actual, color = paste0("Correct Label: ", Actual, " DE genes")),
               linetype="longdash", size=0.6, show.legend = T) +
    scale_color_manual(values = "blue", name = "Number of DE genes") +
    #geom_text(aes(x=Actual+15, label=paste0(Actual, "DE genes"), y=70), colour="blue", angle=90, size=5) +
    #geom_text(aes(x=Actual+20, label = Actual, y = 60), color = "blue", size = 5) +
    #ggtitle(paste(1000, "Random Permutations:", baseline, "vs", target)) +#, subtitle = paste0("Number of permutations with more than 1000 DE genes:", nrow(subset(tab1, NumDiffGenes > 200)))) +
    #labs(caption = paste0("Percentage of permutations with more DE genes than actual:", nrow(subset(tab1, NumDiffGenes > 200))/nrow(tab1)*100,"%")) +
    xlab("Number of significant genes") +
    theme_bw() +
    theme(aspect.ratio=1,
          plot.title = element_text(hjust = 0.5),
          legend.title = element_text(size=10, hjust = 0.5),
          legend.position = "none",
          axis.title.x = element_blank())
  
  grob1 <- grobTree(textGrob(paste0(nrow(subset(tab1, NumDiffGenes < Actual))/nrow(tab1)*100,"%"), x=0.1,  y=0.95, hjust=0,
                           gp=gpar(col="red", fontsize=15, fontface="italic")))
  grob2 <- grobTree(textGrob(tab1$Actual[1], x=0.8,  y=0.95, hjust=0,
                           gp=gpar(col="blue", fontsize=15)))

  #pdf(paste0("updated_output/Permutation_test_histogram_", baseline, "_vs_", target, ".pdf"))
  #print(p + annotation_custom(grob))
  #dev.off()
  p_list <- append(p_list,list(p+annotation_custom(grob1)+annotation_custom(grob2)))
  #p_list <- append(p_list,list(p) + annotation_custom(grob1) + annotation_custom(grob2))
}

#pdf("C:/Users/doe/projects/new_norm/permutation/concatenated_permutation_plots_noxlim_transpose_noNumbers.pdf", height = 8, width = 15)
pdf("E:/Box drive update/4. EV-RNA/4. Code/2022-02-17 V7/Figures/concatenated_permutation_plots_noxlim_transpose_noNumbers.pdf", height = 8, width = 16)
grid.arrange(grobs = p_list, ncol = 6, as.table = F, right = "1")
dev.off()

