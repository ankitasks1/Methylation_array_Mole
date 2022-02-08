library(ggplot2)
library(gplots)
setwd("/media/ankitv/Archivio2/ankit/Mole/trial_backup_files/Boxplots/hg38/arrayWithseq/All")
data <- read.table("Boxplot_maternalgDMRs_hannahway.txt",header=TRUE)
rownames(data)
rownames(data)=data[,1]
rownames(data)
colnames(data)
data = data[,-1]
data = as.matrix(data)
head(data)
dim(data)
colfunc <- colorRampPalette(c("#291BEB","white", "red"))
heatmap.2(data, Colv = "NA",trace = "none", col = colfunc , lmat=rbind( c(5, 4, 2), c(6,1,3 ) ), lwid=c(1, 3, 2 ),lhei = c(1,7), keysize=1.2, key.par = list(cex=0.5), density.info=c("none"), dendrogram="none", scale = "none", sepwidth=c(0.001, 0.001), cexRow=0.5, font=3, cexCol = 1, margins =c(5,8), breaks = seq(0,1, length.out = 100))




#Add separation color
heatmap.2(data, Colv="NA", Rowv = "NA", trace = "none", colsep = 1:length(data), rowsep = 1:length(data), col = colfunc ,  dendrogram="none", scale = "none", sepwidth=c(0.001, 0.001), sepcolor ="none" , cexRow=1, font=3, cexCol = 1, margins =c(5,8), breaks = seq(0,1, length.out = 100))

#Let it cluster
heatmap.2(data, trace = "none", col = colfunc , lmat=rbind( c(5, 4, 2), c(6,1,3 ) ), lwid=c(1, 3, 2 ),lhei = c(1,7), keysize=1.2, key.par = list(cex=0.5), density.info=c("none"), dendrogram="none", scale = "none", sepwidth=c(0.001, 0.001), cexRow=0.8, font=3, cexCol = 1, margins =c(5,8), breaks = seq(0,1, length.out = 100))
