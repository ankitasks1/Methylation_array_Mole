library(ggplot2)
library(plyr)
setwd("/media/ankitv/Archivio2/ankit/Mole/trial_backup_files/Scatters/hg38/Scatter_Plot_20Kb_RepCpGremoved")
#awk '{print $1"%"$2"%"$3"\t"$4"\t"$5"\t"$6}' /media/ankitv/Archivio2/ankit/Mole/bisulfite-seq/bismark_cov/bedgraph/hg38/array_data_hg38/All/CRS_merge_hg38.txt | grep %% -v > CRS_merge_hg38_coordsmerge%.txt
#>745232 CRS_merge_hg38_coordsmerge%.txt
#fgrep -f ./../20kb+1/CpGmayrepeated_cordmerge.txt CRS_merge_hg38_coordsmerge%.txt -v | awk -F'%' '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' | grep chr > CRS_merge_hg38_repCpGminus.txt
#>745197 CRS_merge_hg38_repCpGminus.txt
#bedtools intersect -wa -wb -a ./../20kb+1/hg38.20kb+1_usethis.txt -b CRS_merge_hg38_repCpGminus.txt > CRS_merge_hg38_repCpGminus.20kb+1.intersected.txt
#awk '{print $1"%"$2"%"$3"\t"$7"\t"$8"\t"$9}' CRS_merge_hg38_repCpGminus.20kb+1.intersected.txt > CRS_merge_hg38_repCpGminus.20kb+1.intersected_rearranged.txt
#20kb+1
CRS_merge_hg38_repCpGminus.20kb_1 =read.table("CRS_merge_hg38_repCpGminus.20kb+1.intersected_rearranged.txt", header = FALSE, stringsAsFactors = FALSE)
head(CRS_merge_hg38_repCpGminus.20kb_1)
aggregate_CRS_merge_hg38_repCpGminus.20kb_1 = aggregate(CRS_merge_hg38_repCpGminus.20kb_1[,2:4],by=list(CRS_merge_hg38_repCpGminus.20kb_1$V1), mean)
head(aggregate_CRS_merge_hg38_repCpGminus.20kb_1)
aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted = aggregate_CRS_merge_hg38_repCpGminus.20kb_1[order(aggregate_CRS_merge_hg38_repCpGminus.20kb_1$Group.1),]
head(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted)
rownames(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted)
rownames(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted)=aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted[,1]
rownames(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted)
colnames(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted)
aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted = aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted[,-1]
aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted = as.matrix(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted)
head(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted)
dim(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted)
Count_CRS_merge_hg38_repCpGminus.20kb_1 <- count(CRS_merge_hg38_repCpGminus.20kb_1, "V1")
aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted <- cbind(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted, Count_CRS_merge_hg38_repCpGminus.20kb_1)
head(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted)
dim(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted)
write.table(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted, "aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE)

#CpG cutoff 3
aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted_cutoff3 <- aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted[which(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted$freq >= 3),]
head(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted_cutoff3)
dim(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted_cutoff3)
write.table(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted_cutoff3, "aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted_cutoff3.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE)
aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted1_cutoff3 <- aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted_cutoff3[,1:3]
colnames(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted1_cutoff3) <- c("Control_Placenta","KHDC3L_BiCHM","AnCHM")
head(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted1_cutoff3)
dim(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted1_cutoff3)
write.table(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted1_cutoff3, "aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted1_cutoff3.txt" , quote = FALSE, append = FALSE, sep = "\t")

summary(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted1_cutoff3)
basicStats(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted1_cutoff3)
wilcox.test(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted1_cutoff3$Control_Placenta, aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted1_cutoff3$KHDC3L_BiCHM)
wilcox.test(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted1_cutoff3$Control_Placenta, aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted1_cutoff3$AnCHM)
aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted2_cutoff3 <- aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted1_cutoff3[,1:2]
write.table(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted2_cutoff3, "aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted2_cutoff3.txt" , quote = FALSE, append = FALSE, sep = "\t")

colfunc <- colorRampPalette(c("#291BEB", "red"))
ggplot(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted2_cutoff3) + geom_hex(aes(x=Control_Placenta, y=KHDC3L_BiCHM), bins = 300) +
  geom_abline(intercept = 0, color="grey")+
  scale_fill_gradientn("", colours = c("blue", "green", "red"),
                       values = scales::rescale(c(0.0, 0.25, 0.5, 0.75, 1.0)))+ theme_bw()
write.table(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted2_cutoff3, "aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted2_cutoff3.txt", quote = FALSE, append = FALSE, sep = "\t")
ggsave("CRS_merge_hg38_repCpGminus.20kb+1.intersected.scatter_cutoff3.svg", width=15*1.25, height=15*1.25, units="cm", dpi=96)

aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted3_cutoff3 <- cbind(rownames(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted1_cutoff3), aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted1_cutoff3)
head(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted3_cutoff3)
dim(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted3_cutoff3)
colnames(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted3_cutoff3) <- c("regions", "Control_Placenta","KHDC3L_BiCHM","AnCHM")
head(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted3_cutoff3)
aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted3_cutoff3[,1]
rownames(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted3_cutoff3)=aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted3_cutoff3[,1]
rownames(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted3_cutoff3)
colnames(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted3_cutoff3)
aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted3_cutoff3_re = aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted3_cutoff3[,-1]
aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted3_cutoff3_re <- data.frame(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted3_cutoff3_re)
head(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted3_cutoff3_re)
mean.placenta <- aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted3_cutoff3_re$Control_Placenta
mean.bichm <- aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted3_cutoff3_re$KHDC3L_BiCHM
plot(mean.placenta, mean.bichm, col="grey", ylim=c(0,1), xlim=c(0,1), pch=16, cex = 0.2, ylab="", xlab="", bty = 'n')
lines(x = c(0,1), y = c(0,1), col = "grey")

#Find Low, Intermediate and High methylated regions in KHDC3L mole
aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted3_cutoff3_ratio 
#CpG Nocutoff
aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted_Nocutoff <- aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted
head(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted_Nocutoff)
dim(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted_Nocutoff)
write.table(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted_Nocutoff, "aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted_Nocutoff.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE)
aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted1_Nocutoff <- aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted_Nocutoff[,1:3]
colnames(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted1_Nocutoff) <- c("Control_Placenta","KHDC3L_BiCHM","AnCHM")
head(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted1_Nocutoff)
dim(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted1_Nocutoff)
summary(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted1_Nocutoff)
wilcox.test(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted1_Nocutoff$Control_Placenta, aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted1_Nocutoff$KHDC3L_BiCHM)
wilcox.test(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted1_Nocutoff$Control_Placenta, aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted1_Nocutoff$AnCHM)
aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted2_Nocutoff <- aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted1_Nocutoff[,1:2]
colfunc <- colorRampPalette(c("#291BEB", "red"))
ggplot(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted2_Nocutoff) + geom_hex(aes(x=Control_Placenta, y=KHDC3L_BiCHM), bins = 300) +
  geom_abline(intercept = 0, color="grey")+
  scale_fill_gradientn("", colours = c("blue", "green", "red"),
                       values = scales::rescale(c(0.0, 0.25, 0.5, 0.75, 1.0)))+ theme_bw()
write.table(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted2_Nocutoff, "aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted2_Nocutoff.txt", quote = FALSE, append = FALSE, sep = "\t")
ggsave("CRS_merge_hg38_repCpGminus.20kb+1.intersected.scatter_Nocutoff.svg", width=15*1.25, height=15*1.25, units="cm", dpi=96)

aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted3_Nocutoff <- cbind(rownames(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted1_Nocutoff), aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted1_Nocutoff)
head(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted3_Nocutoff)
dim(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted3_Nocutoff)
colnames(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted3_Nocutoff) <- c("regions", "Control_Placenta","KHDC3L_BiCHM","AnCHM")
head(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted3_Nocutoff)
aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted3_Nocutoff[,1]
rownames(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted3_Nocutoff)=aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted3_Nocutoff[,1]
rownames(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted3_Nocutoff)
colnames(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted3_Nocutoff)
aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted3_Nocutoff_re = aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted3_Nocutoff[,-1]
aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted3_Nocutoff_re <- data.frame(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted3_Nocutoff_re)
head(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted3_Nocutoff_re)
mean.placenta <- aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted3_Nocutoff_re$Control_Placenta
mean.bichm <- aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted3_Nocutoff_re$KHDC3L_BiCHM
plot(mean.placenta, mean.bichm, col="grey", ylim=c(0,1), xlim=c(0,1), pch=16, cex = 0.2, ylab="", xlab="", bty = 'n')
lines(x = c(0,1), y = c(0,1), col = "grey")


#CpG cutoff 10
aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted_cutoff10 <- aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted[which(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted$freq >= 10),]
head(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted_cutoff10)
dim(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted_cutoff10)
write.table(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted_cutoff10, "aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted_cutoff10.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE)
aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted1_cutoff10 <- aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted_cutoff10[,1:3]
colnames(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted1_cutoff10) <- c("Control_Placenta","KHDC3L_BiCHM","AnCHM")
head(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted1_cutoff10)
dim(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted1_cutoff10)
summary(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted1_cutoff10)
wilcox.test(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted1_cutoff10$Control_Placenta, aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted1_cutoff10$KHDC3L_BiCHM)
wilcox.test(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted1_cutoff10$Control_Placenta, aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted1_cutoff10$AnCHM)
aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted2_cutoff10 <- aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted1_cutoff10[,1:2]
colfunc <- colorRampPalette(c("#291BEB", "red"))
ggplot(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted2_cutoff10) + geom_hex(aes(x=Control_Placenta, y=KHDC3L_BiCHM), bins = 300) +
  geom_abline(intercept = 0, color="grey")+
  scale_fill_gradientn("", colours = c("blue", "green", "red"),
                       values = scales::rescale(c(0.0, 0.25, 0.5, 0.75, 1.0)))+ theme_bw()
# write.table(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted2_cutoff10, "aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted2_cutoff10.txt", quote = FALSE, append = FALSE, sep = "\t")
ggsave("CRS_merge_hg38_repCpGminus.20kb+1.intersected.scatter_cutoff10.svg", width=15*1.25, height=15*1.25, units="cm", dpi=96)

aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted3_cutoff10 <- cbind(rownames(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted1_cutoff10), aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted1_cutoff10)
head(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted3_cutoff10)
dim(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted3_cutoff10)
colnames(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted3_cutoff10) <- c("regions", "Control_Placenta","KHDC3L_BiCHM","AnCHM")
head(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted3_cutoff10)
aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted3_cutoff10[,1]
rownames(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted3_cutoff10)=aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted3_cutoff10[,1]
rownames(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted3_cutoff10)
colnames(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted3_cutoff10)
aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted3_cutoff10_re = aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted3_cutoff10[,-1]
aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted3_cutoff10_re <- data.frame(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted3_cutoff10_re)
head(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted3_cutoff10_re)
mean.placenta <- aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted3_cutoff10_re$Control_Placenta
mean.bichm <- aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted3_cutoff10_re$KHDC3L_BiCHM
plot(mean.placenta, mean.bichm, col="grey", ylim=c(0,1), xlim=c(0,1), pch=16, cex = 0.2, ylab="", xlab="", bty = 'n')
lines(x = c(0,1), y = c(0,1), col = "grey")
write.table(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted3_cutoff10_re, "aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted3_cutoff10_re.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE)

#--------------------------------------- Overlap with DMR -----------------------------------------------#

#Two files DMRs_avg_meth.txt and classic_gDMRs_avg_meth.txt were created using aggregated_updated_filteredCpGs_donebyR.R and selecting classical ones usinng excel sheet
#Now first import this file DMRs_avg_meth.txt : All DMRs: 43 selected >=3CpGs
DMRs_avg_meth =read.table("DMRs_avg_meth.txt", header = TRUE, stringsAsFactors = FALSE)
head(DMRs_avg_meth)
dim(DMRs_avg_meth)
DMRs_avg_meth[,1]
rownames(DMRs_avg_meth)=DMRs_avg_meth[,1]
rownames(DMRs_avg_meth)
colnames(DMRs_avg_meth)
DMRs_avg_meth = DMRs_avg_meth[,-1]
DMRs_avg_meth <- data.frame(DMRs_avg_meth)
head(DMRs_avg_meth)
dim(DMRs_avg_meth)
dmrs_placenta <- DMRs_avg_meth$Control_Placenta
dmrs_bichm <- DMRs_avg_meth$KHDC3L_BiCHM
par(new=T) # overlap with previous plot (need to define the same xlim and ylim for both plots)
plot(dmrs_placenta, dmrs_bichm, col="red", ylim=c(0,1), xlim=c(0,1), pch=16,cex = 0.5, ylab="", xlab="", bty = 'n')

#Next import this file classic_gDMRs_avg_meth.txt : All classical DMRs: 32 selected >=3CpGs
classical_DMRs_avg_meth =read.table("classic_gDMRs_avg_meth.txt", header = TRUE, stringsAsFactors = FALSE)
head(classical_DMRs_avg_meth)
dim(classical_DMRs_avg_meth)
classical_DMRs_avg_meth[,1]
rownames(classical_DMRs_avg_meth)= classical_DMRs_avg_meth[,1]
rownames(classical_DMRs_avg_meth)
colnames(classical_DMRs_avg_meth)
classical_DMRs_avg_meth = classical_DMRs_avg_meth[,-1]
classical_DMRs_avg_meth <- data.frame(classical_DMRs_avg_meth)
head(classical_DMRs_avg_meth)
cdmrs_placenta <- classical_DMRs_avg_meth$Control_Placenta
cdmrs_bichm <- classical_DMRs_avg_meth$KHDC3L_BiCHM
par(new=T) # overlap with previous plot (need to define the same xlim and ylim for both plots)
plot(cdmrs_placenta, cdmrs_bichm, col="darkblue", ylim=c(0,1), xlim=c(0,1), pch=16,cex = 0.5, ylab="KHDC3L BiCHM", xlab="Control Placenta", main ="Methylation of Control Placenta v/s KHDC3L BiCHM", bty = 'n')

#awk '{print $1"%"$2"%"$3"%"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' CR_more1.5_Hyper_GeneNOTgene.txt > CR_more1.5_Hyper_GeneNOTgene%.txt
#Next import this file CR_more1.5_Hyper_GeneNOTgene.txt : 89 selected >=10CpGs
CR_more1.5_Hyper_GeneNOTgene =read.table("CR_more1.5_Hyper_GeneNOTgene%.txt", header = TRUE, stringsAsFactors = FALSE)
head(CR_more1.5_Hyper_GeneNOTgene)
dim(CR_more1.5_Hyper_GeneNOTgene)
CR_more1.5_Hyper_GeneNOTgene[,1]
rownames(CR_more1.5_Hyper_GeneNOTgene)= CR_more1.5_Hyper_GeneNOTgene[,1]
rownames(CR_more1.5_Hyper_GeneNOTgene)
colnames(CR_more1.5_Hyper_GeneNOTgene)
CR_more1.5_Hyper_GeneNOTgene = CR_more1.5_Hyper_GeneNOTgene[,-1]
CR_more1.5_Hyper_GeneNOTgene <- data.frame(CR_more1.5_Hyper_GeneNOTgene)
head(CR_more1.5_Hyper_GeneNOTgene)
hyperOnepointFive_placenta <- CR_more1.5_Hyper_GeneNOTgene$Control_Placenta
hyperOnepointFive_bichm <- CR_more1.5_Hyper_GeneNOTgene$KHDC3L_BiCHM
par(new=T) # overlap with previous plot (need to define the same xlim and ylim for both plots)
plot(hyperOnepointFive_placenta, hyperOnepointFive_bichm, col="black", ylim=c(0,1), xlim=c(0,1), pch=16,cex = 0.5, ylab="KHDC3L BiCHM", xlab="Control Placenta", main ="Methylation of Control Placenta v/s KHDC3L BiCHM", bty = 'n')

#awk '{print $1"%"$2"%"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' CR_more2.0_Hyper_GeneNOTgene.txt > CR_more2.0_Hyper_GeneNOTgene%.txt
#Next import this file CR_more2.0_Hyper_GeneNOTgene.txt : 6 selected >=10CpGs
CR_more2.0_Hyper_GeneNOTgene =read.table("CR_more2.0_Hyper_GeneNOTgene%.txt", header = TRUE, stringsAsFactors = FALSE)
head(CR_more2.0_Hyper_GeneNOTgene)
dim(CR_more2.0_Hyper_GeneNOTgene)
CR_more2.0_Hyper_GeneNOTgene[,1]
rownames(CR_more2.0_Hyper_GeneNOTgene)= CR_more2.0_Hyper_GeneNOTgene[,1]
rownames(CR_more2.0_Hyper_GeneNOTgene)
colnames(CR_more2.0_Hyper_GeneNOTgene)
CR_more2.0_Hyper_GeneNOTgene = CR_more2.0_Hyper_GeneNOTgene[,-1]
CR_more2.0_Hyper_GeneNOTgene <- data.frame(CR_more2.0_Hyper_GeneNOTgene)
head(CR_more2.0_Hyper_GeneNOTgene)
hyperTwo_placenta <- CR_more2.0_Hyper_GeneNOTgene$Control_Placenta
hyperTwo_bichm <- CR_more2.0_Hyper_GeneNOTgene$KHDC3L_BiCHM
par(new=T) # overlap with previous plot (need to define the same xlim and ylim for both plots)
plot(hyperTwo_placenta, hyperTwo_bichm, col="orange", ylim=c(0,1), xlim=c(0,1), pch=16,cex = 0.5, ylab="KHDC3L BiCHM", xlab="Control Placenta", main ="Methylation of Control Placenta v/s KHDC3L BiCHM", bty = 'n')

####################---------------------Extra commands -----------------------------------------##

#DMRs 63 = Mat, Pat, Secondary, placenta
#bedtools intersect -wa -wb -a aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted2.cords -b /media/ankitv/Archivio2/ankit/Mole/bisulfite-seq/bismark_cov/dmrs/hg38/Monk_human_ICRs-oocytegDMRs_+_Hanna_Placenta_gDMRs+PatnSecondaryDMRs_hg38.bed | awk '{print $1"%"$2"%"$3"\t"$7}' | sort -k1,1 -u  > 20kb_regions_containing_DMRs.txt 
#DMRs_in_20kb_regions =read.table("20kb_regions_containing_DMRs.txt", header = FALSE, stringsAsFactors = FALSE)
#head(DMRs_in_20kb_regions)
#dim(DMRs_in_20kb_regions)
#colnames(DMRs_in_20kb_regions) <- c("regions","DMRs")
#head(DMRs_in_20kb_regions)
#DMRs_in_20kb_regions_1 <- DMRs_in_20kb_regions[,1]
x = aggregate_CRS_merge_hg38_repCpGminus.20kb_1_DMRs_others$Control_Placenta
y = aggregate_CRS_merge_hg38_repCpGminus.20kb_1_DMRs_others$KHDC3L_BiCHM
z = aggregate_CRS_merge_hg38_repCpGminus.20kb_1_DMRs_others$DMRs
ggplot(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_DMRs_others, aes(x=Control_Placenta, y=KHDC3L_BiCHM)) +
  geom_point(aes(color=DMRs, shape=DMRs, size=DMRs)) + 
  geom_abline(intercept = 0, color="grey")+
  scale_shape_manual(values=c(16, 16))+ 
  scale_color_manual(values=c('#CC0000','#C0C0C0'))+
  scale_size_manual(values=c(2,0.2))+
  theme(legend.position="top") + theme_bw()

ggplot(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_DMRs_others, aes(x=Control_Placenta, y=KHDC3L_BiCHM, color = "grey")) +
  geom_point(shape=16, color="grey", size=0.1)+
  geom_abline(intercept = 0, color="grey") + theme_classic()

ggsave("CR_merge.scatter.svg", width=20*1.25, height=20*1.25, units="cm", dpi=96)

aggregate_CRS_merge_hg38_repCpGminus.20kb_1_DMRs = merge(aggregate_CRS_merge_hg38_repCpGminus.20kb_1, DMRs_in_20kb_regions, by="regions", all.x=TRUE)
head(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_DMRs)
dim(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_DMRs)
write.table(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_DMRs, "aggregate_CRS_merge_hg38_repCpGminus.20kb_1_DMRs.txt" , quote = FALSE, append = FALSE, sep = "\t")
setwd("/media/ankitv/Archivio2/ankit/Mole/trial_backup_files/Scatters/hg38/Scatter_Plot_20Kb_RepCpGremoved")
aggregate_CRS_merge_hg38_repCpGminus.20kb_1_DMRs_others = read.table("aggregate_CRS_merge_hg38_repCpGminus.20kb_1_DMRs_others.txt", header = TRUE, stringsAsFactors = FALSE)
ggplot(aggregate_CRS_rearranged_mat_gDMRs_sorted_counted1, aes(x=Control_Placenta, y=KHDC3L_BiCHM, color = "red")) +
  geom_point(shape=16, size=1)+  theme_classic() + xlim(0,1) + ylim(0,1)

#20kb. Just for testing purpose
CRS_merge_hg38_repCpGminus.20kb =read.table("CRS_merge_hg38_repCpGminus.20kb.intersected_rearranged.txt", header = FALSE, stringsAsFactors = FALSE)
head(CRS_merge_hg38_repCpGminus.20kb)
aggregate_CRS_merge_hg38_repCpGminus.20kb = aggregate(CRS_merge_hg38_repCpGminus.20kb[,2:4],by=list(CRS_merge_hg38_repCpGminus.20kb$V1), mean)
head(aggregate_CRS_merge_hg38_repCpGminus.20kb)
aggregate_CRS_merge_hg38_repCpGminus.20kb_sorted = aggregate_CRS_merge_hg38_repCpGminus.20kb[order(aggregate_CRS_merge_hg38_repCpGminus.20kb$Group.1),]
head(aggregate_CRS_merge_hg38_repCpGminus.20kb_sorted)
rownames(aggregate_CRS_merge_hg38_repCpGminus.20kb_sorted)
rownames(aggregate_CRS_merge_hg38_repCpGminus.20kb_sorted)=aggregate_CRS_merge_hg38_repCpGminus.20kb_sorted[,1]
rownames(aggregate_CRS_merge_hg38_repCpGminus.20kb_sorted)
colnames(aggregate_CRS_merge_hg38_repCpGminus.20kb_sorted)
aggregate_CRS_merge_hg38_repCpGminus.20kb_sorted = aggregate_CRS_merge_hg38_repCpGminus.20kb_sorted[,-1]
aggregate_CRS_merge_hg38_repCpGminus.20kb_sorted = as.matrix(aggregate_CRS_merge_hg38_repCpGminus.20kb_sorted)
head(aggregate_CRS_merge_hg38_repCpGminus.20kb_sorted)
dim(aggregate_CRS_merge_hg38_repCpGminus.20kb_sorted)
Count_CRS_merge_hg38_repCpGminus.20kb <- count(CRS_merge_hg38_repCpGminus.20kb, "V1")
aggregate_CRS_merge_hg38_repCpGminus.20kb_sorted_counted <- cbind(aggregate_CRS_merge_hg38_repCpGminus.20kb_sorted, Count_CRS_merge_hg38_repCpGminus.20kb)
head(aggregate_CRS_merge_hg38_repCpGminus.20kb_sorted_counted)
dim(aggregate_CRS_merge_hg38_repCpGminus.20kb_sorted_counted)
aggregate_CRS_merge_hg38_repCpGminus.20kb_sorted_counted <- aggregate_CRS_merge_hg38_repCpGminus.20kb_sorted_counted[which(aggregate_CRS_merge_hg38_repCpGminus.20kb_sorted_counted$freq >= 3),]
aggregate_CRS_merge_hg38_repCpGminus.20kb_sorted_counted1 <- aggregate_CRS_merge_hg38_repCpGminus.20kb_sorted_counted[,1:3]
colnames(aggregate_CRS_merge_hg38_repCpGminus.20kb_sorted_counted1) <- c("Control_Placenta","KHDC3L_BiCHM","AnCHM")
head(aggregate_CRS_merge_hg38_repCpGminus.20kb_sorted_counted1)
dim(aggregate_CRS_merge_hg38_repCpGminus.20kb_sorted_counted1)
aggregate_CRS_merge_hg38_repCpGminus.20kb_sorted_counted2 <- aggregate_CRS_merge_hg38_repCpGminus.20kb_sorted_counted1[,1:2]
colfunc <- colorRampPalette(c("#291BEB", "red"))
ggplot(aggregate_CRS_merge_hg38_repCpGminus.20kb_sorted_counted2) + geom_hex(aes(x=Control_Placenta, y=KHDC3L_BiCHM), bins = 300) +
  geom_abline(intercept = 0, color="grey")+
  scale_fill_gradientn("", colours = c("blue", "green", "red"),
                       values = scales::rescale(c(0.0, 0.25, 0.5, 0.75, 1.0)))+ theme_bw()
write.table(aggregate_CRS_merge_hg38_repCpGminus.20kb_sorted_counted2, "aggregate_CRS_merge_hg38_repCpGminus.20kb_sorted_counted2.txt" , quote = FALSE, append = FALSE, sep = "\t")
ggsave("CRS_merge_hg38_repCpGminus.20kb.intersected.scatter.svg", width=15*1.25, height=15*1.25, units="cm", dpi=96)



#---------------------------------Merge Gene by Comma---------------------------------
#Novel DMRs
head(aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted_cutoff10)
window20kbcutoff10array <- aggregate_CRS_merge_hg38_repCpGminus.20kb_1_sorted_counted_cutoff10
head(window20kbcutoff10array)
window20kbcutoff10array <- window20kbcutoff10array[,c(4,5,1:3)]
head(window20kbcutoff10array)
tail(window20kbcutoff10array)
colnames(window20kbcutoff10array) <- c("20Kbbins","Freq","ControlPlacenta","KHDC3LBiCHM","AnCHM")
window20kbcutoff10array["R.C"] <- data.frame(cbind(window20kbcutoff10array$KHDC3LBiCHM / window20kbcutoff10array$ControlPlacenta))
window20kbcutoff10array["S.C"] <- data.frame(cbind(window20kbcutoff10array$AnCHM / window20kbcutoff10array$ControlPlacenta))
head(window20kbcutoff10array)
window20kbcutoff10array_Hypo0.65 <- window20kbcutoff10array[window20kbcutoff10array$R.C<=0.65,]
head(window20kbcutoff10array_Hypo0.65)
dim(window20kbcutoff10array_Hypo0.65)

window20kbcutoff10array_Hyper2 <- window20kbcutoff10array[window20kbcutoff10array$R.C>=2.00,]
head(window20kbcutoff10array_Hyper2)
dim(window20kbcutoff10array_Hyper2)

window20kbcutoff10array_Hyper1.5 <- window20kbcutoff10array[window20kbcutoff10array$R.C>=1.50,]
head(window20kbcutoff10array_Hyper1.5)
dim(window20kbcutoff10array_Hyper1.5)
window20kbcutoff10array_Hyper1.5 <- window20kbcutoff10array_Hyper1.5[,c(1,3:7,2)]
write.table(window20kbcutoff10array_Hyper1.5, "window20kbcutoff10array_Hyper1.5.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = F, col.names = F)

#grep NA CR_more2.0_Hyper_GeneNOTgene.txt | sort -k1,1 -k2,2n > CR_more2.0_Hyper_GeneNOTgene_.txt
#grep NA CR_more1.5_Hyper_GeneNOTgene.txt | sort -k1,1 -k2,2n > CR_more1.5_Hyper_GeneNOTgene_.txt

#cat CR_more1.5_Hyper_GeneNOTgene_.txt CR_more2.0_Hyper_GeneNOTgene_.txt > CR_more1.5and2.0_Hyper_GeneNOTgene_.txt

#Monk_human_ICRs-oocytegDMRs_+_Hanna_Placenta_gDMRs_hg38.bed 
bedtools intersect -wa -wb -a /media/ankitv/Archivio2/ankit/Mole/bisulfite-seq/bismark_cov/dmrs/hg38/Monk_human_ICRs-oocytegDMRs_+_Hanna_Placenta_gDMRs_hg38.bed -b CR_more1.5and2.0_Hyper_GeneNOTgene_.txt > ./hyperfeature/CR_more1.5and2.0_Hyper_GeneNOTgene_updated_Monk_oocyte_Hanna_placenta_gDMRs.txt

#Oocyte-meth-CGIs_Okae-above70_unmeth-placenta_Bernstein-below10_GRCh38.txt
bedtools intersect -wa -wb -a /media/ankitv/Archivio2/ankit/Mole/bisulfite-seq/bismark_cov/cgis/GRCh38/Oocyte-meth-CGIs_Okae-above70_unmeth-placenta_Bernstein-below10_GRCh38.txt -b CR_more1.5and2.0_Hyper_GeneNOTgene_.txt > ./hyperfeature/CR_more1.5and2.0_Hyper_GeneNOTgene_OmUp.txt

#Hum_Placenta_WGBS_Bernstein_2010_20kB-windows_no-CGIs_GRCh38.txt
bedtools intersect -wa -wb -a /media/ankitv/Archivio2/ankit/Mole/bisulfite-seq/bismark_cov/cgis/GRCh38/Hum_Placenta_WGBS_Bernstein_2010_20kB-windows_no-CGIs_GRCh38.txt -b CR_more1.5and2.0_Hyper_GeneNOTgene_.txt > ./hyperfeature/CR_more1.5and2.0_Hyper_GeneNOTgene_20kB-windows_no-CGIs.txt

#Hum_Placenta_WGBS_Bernstein_2010_Meth_CGIs_above70_GRCh38.txt
bedtools intersect -wa -wb -a /media/ankitv/Archivio2/ankit/Mole/bisulfite-seq/bismark_cov/cgis/GRCh38/Hum_Placenta_WGBS_Bernstein_2010_Meth_CGIs_above70_GRCh38.txt -b CR_more1.5and2.0_Hyper_GeneNOTgene_.txt > ./hyperfeature/CR_more1.5and2.0_Hyper_GeneNOTgene_hum_placenta_Meth_CGIs.txt

#Hum_Placenta_WGBS_Bernstein_2010_Unmeth_CGIs_below10_GRCh38.txt
bedtools intersect -wa -wb -a /media/ankitv/Archivio2/ankit/Mole/bisulfite-seq/bismark_cov/cgis/GRCh38/Hum_Placenta_WGBS_Bernstein_2010_Unmeth_CGIs_below10_GRCh38.txt -b CR_more1.5and2.0_Hyper_GeneNOTgene_.txt > ./hyperfeature/CR_more1.5and2.0_Hyper_GeneNOTgene_hum_placenta_Unmeth_CGIs.txt

#Non-CGI_20kB-windows_unmeth-in-Okae-oocytes_meth-in-Bernstein-placenta_GRCh38.txt
bedtools intersect -wa -wb -a /media/ankitv/Archivio2/ankit/Mole/bisulfite-seq/bismark_cov/cgis/GRCh38/Non-CGI_20kB-windows_unmeth-in-Okae-oocytes_meth-in-Bernstein-placenta_GRCh38.txt -b CR_more1.5and2.0_Hyper_GeneNOTgene_.txt > ./hyperfeature/CR_more1.5and2.0_Hyper_GeneNOTgene_20kb-windows_no-CGIs_specific.txt

#PmCGIs_20kB-windows_unmeth-in-Okae-oocytes_meth-in-Bernstein-placenta_GRCh38.txt
bedtools intersect -wa -wb -a /media/ankitv/Archivio2/ankit/Mole/bisulfite-seq/bismark_cov/cgis/GRCh38/PmCGIs_20kB-windows_unmeth-in-Okae-oocytes_meth-in-Bernstein-placenta_GRCh38.txt -b CR_more1.5and2.0_Hyper_GeneNOTgene_.txt > ./hyperfeature/CR_more1.5and2.0_Hyper_GeneNOTgene_PmCGI_meth.txt

#Remove duplicates 20kb regions  due to intersection 
awk '{print $4"%"$5"%"$6}' CR_more1.5and2.0_Hyper_GeneNOTgene_20kb-windows_no-CGIs_specific.txt | sort -k1,1 -u | awk -F'%' '{print $1"\t"$2"\t"$3}' > CR_more1.5and2.0_Hyper_GeneNOTgene_20kb-windows_no-CGIs_specific.dedup.txt
awk '{print $4"%"$5"%"$6}' CR_more1.5and2.0_Hyper_GeneNOTgene_20kB-windows_no-CGIs.txt | sort -k1,1 -u | awk -F'%' '{print $1"\t"$2"\t"$3}' > CR_more1.5and2.0_Hyper_GeneNOTgene_20kB-windows_no-CGIs.dedup.txt
awk '{print $4"%"$5"%"$6}' CR_more1.5and2.0_Hyper_GeneNOTgene_hum_placenta_Meth_CGIs.txt | sort -k1,1 -u | awk -F'%' '{print $1"\t"$2"\t"$3}' > CR_more1.5and2.0_Hyper_GeneNOTgene_hum_placenta_Meth_CGIs.dedup.txt
awk '{print $4"%"$5"%"$6}' CR_more1.5and2.0_Hyper_GeneNOTgene_hum_placenta_Unmeth_CGIs.txt | sort -k1,1 -u | awk -F'%' '{print $1"\t"$2"\t"$3}' > CR_more1.5and2.0_Hyper_GeneNOTgene_hum_placenta_Unmeth_CGIs.dedup.txt
awk '{print $4"%"$5"%"$6}' CR_more1.5and2.0_Hyper_GeneNOTgene_OmUp.txt | sort -k1,1 -u | awk -F'%' '{print $1"\t"$2"\t"$3}' > CR_more1.5and2.0_Hyper_GeneNOTgene_OmUp.dedup.txt
awk '{print $4"%"$5"%"$6}' CR_more1.5and2.0_Hyper_GeneNOTgene_PmCGI_meth.txt | sort -k1,1 -u | awk -F'%' '{print $1"\t"$2"\t"$3}' > CR_more1.5and2.0_Hyper_GeneNOTgene_PmCGI_meth.dedup.txt
awk '{print $4"%"$5"%"$6}' CR_more1.5and2.0_Hyper_GeneNOTgene_updated_Monk_oocyte_Hanna_placenta_gDMRs.txt | sort -k1,1 -u | awk -F'%' '{print $1"\t"$2"\t"$3}' > CR_more1.5and2.0_Hyper_GeneNOTgene_updated_Monk_oocyte_Hanna_placenta_gDMRs.dedup.txt