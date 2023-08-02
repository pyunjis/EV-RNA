
#For Coverage Figures
library(trackViewer)
library(ggcoverage)
library(data.table)
library(ggplot2)
library(rtracklayer)

gtf.gr <- rtracklayer::import.gff(con = "E:/4. EV-RNA/4. EV-RNA/4. Code/2023-07-23 FINAL/Files/gencode.v42.chr_patch_hapl_scaff.annotation.gtf")
track.df <- LoadTrackFile(track.folder = "E:/4. EV-RNA/4. EV-RNA/4. Code/2023-07-23 FINAL/Files/", format = "bw")
track.df[which(track.df$Group %like% "P1a.bw"),6] <- "Plasma 1"
track.df[which(track.df$Group %like% "p2b.bw"),6] <- "Plasma 2"
track.df[which(track.df$Group %like% "p3c.bw"),6] <- "Plasma 3"
track.df$Type <- track.df$Group

ALB <- ggcoverage(data = track.df, gtf.gr = gtf.gr, gene.name = "ALB", region = "chr4:73404193-73421482", color = c("#8dd3c7", "#d8c965", "#ffed6f"), range.position = "out", extend = 250) + scale_y_continuous(breaks=c(0,25,50,75,100)) + geom_gene(gtf.gr, label.size = 4, arrow.num = 0)

ACTB <- ggcoverage(data = track.df, gtf.gr = gtf.gr, gene.name = "ACTB", region = "chr7:5527123-5530653", color = c("#8dd3c7", "#d8c965", "#ffed6f"), range.position = "out", extend=50) + scale_y_continuous(breaks = c(0,250, 500, 750, 1000, 1500, 2000)) + geom_gene(gtf.gr, label.size = 4, arrow.num = 0)


pdf("E:/4. EV-RNA/4. EV-RNA/4. Code/2023-07-23 FINAL/Figures/Supp_Figure_4B_ALB_coverage.pdf", width = 12, height = 12)
print(ALB)
dev.off()

pdf("E:/4. EV-RNA/4. EV-RNA/4. Code/2023-07-23 FINAL/Figures/Supp_Figure_4A_ACTB_coverage.pdf", width = 12, height = 12)
print(ACTB)
dev.off()
