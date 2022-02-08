#Extract positions from original files
awk '{print $2"\t"$3"\t"$4}' Hum_Placenta_WGBS_Bernstein_2010_20kB-windows_no-CGIs.txt > Hum_Placenta_WGBS_Bernstein_2010_20kB-windows_no-CGIs_hg19.txt
#122734 Hum_Placenta_WGBS_Bernstein_2010_20kB-windows_no-CGIs_hg19.txt

awk '{print $4"\t"$5"\t"$6}' Hum_Placenta_WGBS_Bernstein_2010_Meth_CGIs_above70.txt > Hum_Placenta_WGBS_Bernstein_2010_Meth_CGIs_above70_hg19.txt
#4077 Hum_Placenta_WGBS_Bernstein_2010_Meth_CGIs_above70_hg19.txt

awk '{print $4"\t"$5"\t"$6}' Hum_Placenta_WGBS_Bernstein_2010_Unmeth_CGIs_below10.txt > Hum_Placenta_WGBS_Bernstein_2010_Unmeth_CGIs_below10_hg19.txt
#12909 Hum_Placenta_WGBS_Bernstein_2010_Unmeth_CGIs_below10_hg19.txt

awk '{print $2"\t"$3"\t"$4}' Non-CGI_20kB-windows_unmeth-in-Okae-oocytes_meth-in-Bernstein-placenta.txt > Non-CGI_20kB-windows_unmeth-in-Okae-oocytes_meth-in-Bernstein-placenta_hg19.txt
#15442 Non-CGI_20kB-windows_unmeth-in-Okae-oocytes_meth-in-Bernstein-placenta_hg19.txt

awk '{print $4"\t"$5"\t"$6}' Placenta_methylated_CGIs_above70_unmet-Okae-oocyte_below20.txt > PmCGIs_20kB-windows_unmeth-in-Okae-oocytes_meth-in-Bernstein-placenta_hg19.txt
#3619 PmCGIs_20kB-windows_unmeth-in-Okae-oocytes_meth-in-Bernstein-placenta_hg19.txt


#Move name
mv Hum_Placenta_WGBS_Bernstein_2010_20kB-windows_no-CGIs_hg19.txt Hum_Placenta_WGBS_Bernstein_2010_20kB-windows_no-CGIs_GRCh37.txt
mv Hum_Placenta_WGBS_Bernstein_2010_Meth_CGIs_above70_hg19.txt Hum_Placenta_WGBS_Bernstein_2010_Meth_CGIs_above70_GRCh37.txt
mv Hum_Placenta_WGBS_Bernstein_2010_Unmeth_CGIs_below10_hg19.txt Hum_Placenta_WGBS_Bernstein_2010_Unmeth_CGIs_below10_GRCh37.txt
mv Non-CGI_20kB-windows_unmeth-in-Okae-oocytes_meth-in-Bernstein-placenta_hg19.txt Non-CGI_20kB-windows_unmeth-in-Okae-oocytes_meth-in-Bernstein-placenta_GRCh37.txt
mv PmCGIs_20kB-windows_unmeth-in-Okae-oocytes_meth-in-Bernstein-placenta_hg19.txt PmCGIs_20kB-windows_unmeth-in-Okae-oocytes_meth-in-Bernstein-placenta_GRCh37.txt


#Convert GRCh37 to GRCh38
python ~/tools_av/CrossMap-0.2.5/bin/CrossMap.py bed /home/ankitv/ref_av/GRCh37/GRCh37_to_GRCh38.chain.gz Hum_Placenta_WGBS_Bernstein_2010_20kB-windows_no-CGIs_GRCh37.txt > Hum_Placenta_WGBS_Bernstein_2010_20kB-windows_no-CGIs_GRCh37ToGRCh38.chain
python ~/tools_av/CrossMap-0.2.5/bin/CrossMap.py bed /home/ankitv/ref_av/GRCh37/GRCh37_to_GRCh38.chain.gz Hum_Placenta_WGBS_Bernstein_2010_Meth_CGIs_above70_GRCh37.txt > Hum_Placenta_WGBS_Bernstein_2010_Meth_CGIs_above70_GRCh37ToGRCh38.chain
python ~/tools_av/CrossMap-0.2.5/bin/CrossMap.py bed /home/ankitv/ref_av/GRCh37/GRCh37_to_GRCh38.chain.gz Hum_Placenta_WGBS_Bernstein_2010_Unmeth_CGIs_below10_GRCh37.txt > Hum_Placenta_WGBS_Bernstein_2010_Unmeth_CGIs_below10_GRCh37ToGRCh38.chain
python ~/tools_av/CrossMap-0.2.5/bin/CrossMap.py bed /home/ankitv/ref_av/GRCh37/GRCh37_to_GRCh38.chain.gz Non-CGI_20kB-windows_unmeth-in-Okae-oocytes_meth-in-Bernstein-placenta_GRCh37.txt > Non-CGI_20kB-windows_unmeth-in-Okae-oocytes_meth-in-Bernstein-placenta_GRCh37ToGRCh38.chain
python ~/tools_av/CrossMap-0.2.5/bin/CrossMap.py bed /home/ankitv/ref_av/GRCh37/GRCh37_to_GRCh38.chain.gz PmCGIs_20kB-windows_unmeth-in-Okae-oocytes_meth-in-Bernstein-placenta_GRCh37.txt > PmCGIs_20kB-windows_unmeth-in-Okae-oocytes_meth-in-Bernstein-placenta_GRCh37ToGRCh38.chain


#Extract GRCh38 coordinates
grep Fail Hum_Placenta_WGBS_Bernstein_2010_20kB-windows_no-CGIs_GRCh37ToGRCh38.chain -v | awk '{print "chr"$5"\t"$6"\t"$7}' | sort -k1,1 -k2,2n > Hum_Placenta_WGBS_Bernstein_2010_20kB-windows_no-CGIs_GRCh38.txt
grep Fail Hum_Placenta_WGBS_Bernstein_2010_Meth_CGIs_above70_GRCh37ToGRCh38.chain -v | awk '{print "chr"$5"\t"$6"\t"$7}' | sort -k1,1 -k2,2n > Hum_Placenta_WGBS_Bernstein_2010_Meth_CGIs_above70_GRCh38.txt
grep Fail Hum_Placenta_WGBS_Bernstein_2010_Unmeth_CGIs_below10_GRCh37ToGRCh38.chain -v | awk '{print "chr"$5"\t"$6"\t"$7}'| sort -k1,1 -k2,2n > Hum_Placenta_WGBS_Bernstein_2010_Unmeth_CGIs_below10_GRCh38.txt
grep Fail Non-CGI_20kB-windows_unmeth-in-Okae-oocytes_meth-in-Bernstein-placenta_GRCh37ToGRCh38.chain -v | awk '{print "chr"$5"\t"$6"\t"$7}' | sort -k1,1 -k2,2n > Non-CGI_20kB-windows_unmeth-in-Okae-oocytes_meth-in-Bernstein-placenta_GRCh38.txt
grep Fail PmCGIs_20kB-windows_unmeth-in-Okae-oocytes_meth-in-Bernstein-placenta_GRCh37ToGRCh38.chain -v | awk '{print "chr"$5"\t"$6"\t"$7}' | sort -k1,1 -k2,2n > PmCGIs_20kB-windows_unmeth-in-Okae-oocytes_meth-in-Bernstein-placenta_GRCh38.txt

mv *GRCh38.txt ./../GRCh38/

#Prepare 6th file
awk '{print "chr"$4"\t"$5"\t"$6}' ./../hg38/Oocyte-meth-CGIs_Okae-above70_unmeth-placenta_Bernstein-below10_GRCh38.txt > Oocyte-meth-CGIs_Okae-above70_unmeth-placenta_Bernstein-below10_GRCh38.txt
#112 Oocyte-meth-CGIs_Okae-above70_unmeth-placenta_Bernstein-below10_GRCh38.txt
wc -l Oocyte-meth-CGIs_Okae-above70_unmeth-placenta_Bernstein-below10_GRCh38.txt

wc -l Hum_Placenta_WGBS_Bernstein_2010_20kB-windows_no-CGIs_GRCh38.txt

wc -l Hum_Placenta_WGBS_Bernstein_2010_Meth_CGIs_above70_GRCh38.txt

wc -l Hum_Placenta_WGBS_Bernstein_2010_Unmeth_CGIs_below10_GRCh38.txt

wc -l Non-CGI_20kB-windows_unmeth-in-Okae-oocytes_meth-in-Bernstein-placenta_GRCh38.txt

wc -l PmCGIs_20kB-windows_unmeth-in-Okae-oocytes_meth-in-Bernstein-placenta_GRCh38.txt

wc -l Oocyte-meth-CGIs_Okae-above70_unmeth-placenta_Bernstein-below10_GRCh38.txt
112 Oocyte-meth-CGIs_Okae-above70_unmeth-placenta_Bernstein-below10_GRCh38.txt

wc -l Hum_Placenta_WGBS_Bernstein_2010_20kB-windows_no-CGIs_GRCh38.txt
135608 Hum_Placenta_WGBS_Bernstein_2010_20kB-windows_no-CGIs_GRCh38.txt

wc -l Hum_Placenta_WGBS_Bernstein_2010_Meth_CGIs_above70_GRCh38.txt
4235 Hum_Placenta_WGBS_Bernstein_2010_Meth_CGIs_above70_GRCh38.txt

wc -l Hum_Placenta_WGBS_Bernstein_2010_Unmeth_CGIs_below10_GRCh38.txt
13044 Hum_Placenta_WGBS_Bernstein_2010_Unmeth_CGIs_below10_GRCh38.txt

wc -l Non-CGI_20kB-windows_unmeth-in-Okae-oocytes_meth-in-Bernstein-placenta_GRCh38.txt
16617 Non-CGI_20kB-windows_unmeth-in-Okae-oocytes_meth-in-Bernstein-placenta_GRCh38.txt

wc -l PmCGIs_20kB-windows_unmeth-in-Okae-oocytes_meth-in-Bernstein-placenta_GRCh38.txt
3774 PmCGIs_20kB-windows_unmeth-in-Okae-oocytes_meth-in-Bernstein-placenta_GRCh38.txt
#Bedtools intersect with CRS
#Monk oocyte Placenta Hannah = mat gDMRs
bedtools intersect -wa -wb -a /media/ankitv/Archivio2/ankit/Mole/bisulfite-seq/bismark_cov/bedgraph/hg38/array_data_hg38/All/CRS_merge_hg38.txt -b /media/ankitv/Archivio2/ankit/Mole/bisulfite-seq/bismark_cov/dmrs/hg38/Monk_human_ICRs-oocytegDMRs_+_Hanna_Placenta_gDMRs_hg38.bed > CRS_merge_hg38_updated_Monk_oocyte_Hanna_placenta_gDMRs.txt
awk '{print $10"\t"$4"\t"$5"\t"$6}' CRS_merge_hg38_updated_Monk_oocyte_Hanna_placenta_gDMRs.txt > CRS_merge_hg38_updated_Monk_oocyte_Hanna_placenta_gDMRs_rearranged.txt

#Oocyte-meth-CGIs_Okae-above70_unmeth-placenta_Bernstein-below10_GRCh38.txt
bedtools intersect -wa -wb -a /media/ankitv/Archivio2/ankit/Mole/bisulfite-seq/bismark_cov/cgis/GRCh38/Oocyte-meth-CGIs_Okae-above70_unmeth-placenta_Bernstein-below10_GRCh38.txt -b /media/ankitv/Archivio2/ankit/Mole/bisulfite-seq/bismark_cov/bedgraph/hg38/array_data_hg38/All/CRS_merge_hg38.txt > CRS_merge_hg38_OmUp.txt
awk '{print $1"%"$2"%"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' CRS_merge_hg38_OmUp.txt > CRS_merge_hg38_OmUp_.txt
awk '{print $1"\t"$5"\t"$6"\t"$7}' CRS_merge_hg38_OmUp_.txt > CRS_merge_hg38_OmUp_rearranged.txt

#Hum_Placenta_WGBS_Bernstein_2010_20kB-windows_no-CGIs_GRCh38.txt
bedtools intersect -wa -wb -a /media/ankitv/Archivio2/ankit/Mole/bisulfite-seq/bismark_cov/cgis/GRCh38/Hum_Placenta_WGBS_Bernstein_2010_20kB-windows_no-CGIs_GRCh38.txt -b /media/ankitv/Archivio2/ankit/Mole/bisulfite-seq/bismark_cov/bedgraph/hg38/array_data_hg38/All/CRS_merge_hg38.txt > CRS_merge_hg38_20kB-windows_no-CGIs.txt
awk '{print $1"%"$2"%"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' CRS_merge_hg38_20kB-windows_no-CGIs.txt > CRS_merge_hg38_20kB-windows_no-CGIs_.txt
awk '{print $1"\t"$5"\t"$6"\t"$7}' CRS_merge_hg38_20kB-windows_no-CGIs_.txt > CRS_merge_hg38_20kB-windows_no-CGIs_rearranged.txt

#Hum_Placenta_WGBS_Bernstein_2010_Meth_CGIs_above70_GRCh38.txt
bedtools intersect -wa -wb -a /media/ankitv/Archivio2/ankit/Mole/bisulfite-seq/bismark_cov/cgis/GRCh38/Hum_Placenta_WGBS_Bernstein_2010_Meth_CGIs_above70_GRCh38.txt -b /media/ankitv/Archivio2/ankit/Mole/bisulfite-seq/bismark_cov/bedgraph/hg38/array_data_hg38/All/CRS_merge_hg38.txt > CRS_merge_hg38_hum_placenta_Meth_CGIs.txt
awk '{print $1"%"$2"%"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' CRS_merge_hg38_hum_placenta_Meth_CGIs.txt > CRS_merge_hg38_hum_placenta_Meth_CGIs_.txt
awk '{print $1"\t"$5"\t"$6"\t"$7}' CRS_merge_hg38_hum_placenta_Meth_CGIs_.txt > CRS_merge_hg38_hum_placenta_Meth_CGIs_rearranged.txt

#Hum_Placenta_WGBS_Bernstein_2010_Unmeth_CGIs_below10_GRCh38.txt
bedtools intersect -wa -wb -a /media/ankitv/Archivio2/ankit/Mole/bisulfite-seq/bismark_cov/cgis/GRCh38/Hum_Placenta_WGBS_Bernstein_2010_Unmeth_CGIs_below10_GRCh38.txt -b /media/ankitv/Archivio2/ankit/Mole/bisulfite-seq/bismark_cov/bedgraph/hg38/array_data_hg38/All/CRS_merge_hg38.txt > CRS_merge_hg38_hum_placenta_Unmeth_CGIs.txt
awk '{print $1"%"$2"%"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' CRS_merge_hg38_hum_placenta_Unmeth_CGIs.txt > CRS_merge_hg38_hum_placenta_Unmeth_CGIs_.txt
awk '{print $1"\t"$5"\t"$6"\t"$7}' CRS_merge_hg38_hum_placenta_Unmeth_CGIs_.txt > CRS_merge_hg38_hum_placenta_Unmeth_CGIs_rearranged.txt

#Non-CGI_20kB-windows_unmeth-in-Okae-oocytes_meth-in-Bernstein-placenta_GRCh38.txt
bedtools intersect -wa -wb -a /media/ankitv/Archivio2/ankit/Mole/bisulfite-seq/bismark_cov/cgis/GRCh38/Non-CGI_20kB-windows_unmeth-in-Okae-oocytes_meth-in-Bernstein-placenta_GRCh38.txt -b /media/ankitv/Archivio2/ankit/Mole/bisulfite-seq/bismark_cov/bedgraph/hg38/array_data_hg38/All/CRS_merge_hg38.txt > CRS_merge_hg38_20kb-windows_no-CGIs_specific.txt
awk '{print $1"%"$2"%"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' CRS_merge_hg38_20kb-windows_no-CGIs_specific.txt > CRS_merge_hg38_20kb-windows_no-CGIs_specific_.txt
awk '{print $1"\t"$5"\t"$6"\t"$7}' CRS_merge_hg38_20kb-windows_no-CGIs_specific_.txt > CRS_merge_hg38_20kB-windows_no-CGIs_specific_meth_rearranged.txt

#PmCGIs_20kB-windows_unmeth-in-Okae-oocytes_meth-in-Bernstein-placenta_GRCh38.txt
bedtools intersect -wa -wb -a /media/ankitv/Archivio2/ankit/Mole/bisulfite-seq/bismark_cov/cgis/GRCh38/PmCGIs_20kB-windows_unmeth-in-Okae-oocytes_meth-in-Bernstein-placenta_GRCh38.txt -b /media/ankitv/Archivio2/ankit/Mole/bisulfite-seq/bismark_cov/bedgraph/hg38/array_data_hg38/All/CRS_merge_hg38.txt > CRS_merge_hg38_PmCGI_meth.txt
awk '{print $1"%"$2"%"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' CRS_merge_hg38_PmCGI_meth.txt > CRS_merge_hg38_PmCGI_meth_.txt
awk '{print $1"\t"$5"\t"$6"\t"$7}' CRS_merge_hg38_PmCGI_meth_.txt > CRS_merge_hg38_PmCGI_meth_rearranged.txt


library(ggpubr)
library(ggplot2)
library(plyr)
setwd("/media/ankitv/Archivio2/ankit/Mole/trial_backup_files/Boxplots/GRCh38/")

#Mat gDMR = Monk oocyte Hannah Placenta
CRS_rearranged_mat_gDMRs <- read.table("CRS_merge_hg38_updated_Monk_oocyte_Hanna_placenta_gDMRs_rearranged.txt", header = FALSE)
head(CRS_rearranged_mat_gDMRs)
aggregate_CRS_rearranged_mat_gDMRs = aggregate(CRS_rearranged_mat_gDMRs[,2:4],by=list(CRS_rearranged_mat_gDMRs$V1), mean)
head(aggregate_CRS_rearranged_mat_gDMRs)
aggregate_CRS_rearranged_mat_gDMRs_sorted = aggregate_CRS_rearranged_mat_gDMRs[order(aggregate_CRS_rearranged_mat_gDMRs$Group.1),]
head(aggregate_CRS_rearranged_mat_gDMRs_sorted)
rownames(aggregate_CRS_rearranged_mat_gDMRs_sorted)
rownames(aggregate_CRS_rearranged_mat_gDMRs_sorted)=aggregate_CRS_rearranged_mat_gDMRs_sorted[,1]
rownames(aggregate_CRS_rearranged_mat_gDMRs_sorted)
colnames(aggregate_CRS_rearranged_mat_gDMRs_sorted)
aggregate_CRS_rearranged_mat_gDMRs_sorted = aggregate_CRS_rearranged_mat_gDMRs_sorted[,-1]
aggregate_CRS_rearranged_mat_gDMRs_sorted = as.matrix(aggregate_CRS_rearranged_mat_gDMRs_sorted)
head(aggregate_CRS_rearranged_mat_gDMRs_sorted)
dim(aggregate_CRS_rearranged_mat_gDMRs_sorted)
Count_CRS_rearranged_mat_gDMRs <- count(CRS_rearranged_mat_gDMRs, "V1")
aggregate_CRS_rearranged_mat_gDMRs_sorted_counted <- cbind(aggregate_CRS_rearranged_mat_gDMRs_sorted, Count_CRS_rearranged_mat_gDMRs)
head(aggregate_CRS_rearranged_mat_gDMRs_sorted_counted)
dim(aggregate_CRS_rearranged_mat_gDMRs_sorted_counted)
aggregate_CRS_rearranged_mat_gDMRs_sorted_counted <- aggregate_CRS_rearranged_mat_gDMRs_sorted_counted[which(aggregate_CRS_rearranged_mat_gDMRs_sorted_counted$freq >= 3),]
#write.table(aggregate_CRS_rearranged_mat_gDMRs_sorted_counted, "aggregate_CRS_rearranged_mat_gDMRs_sorted_counted.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE)
aggregate_CRS_rearranged_mat_gDMRs_sorted_counted1 <- aggregate_CRS_rearranged_mat_gDMRs_sorted_counted[,1:3]
colnames(aggregate_CRS_rearranged_mat_gDMRs_sorted_counted1) <- c("Control_Placenta","KHDC3L_BiCHM","AnCHM")
head(aggregate_CRS_rearranged_mat_gDMRs_sorted_counted1)
dim(aggregate_CRS_rearranged_mat_gDMRs_sorted_counted1)
write.table(aggregate_CRS_rearranged_mat_gDMRs_sorted_counted1, "aggregate_CRS_rearranged_mat_gDMRs_sorted_counted1.txt" , quote = FALSE, append = FALSE, sep = "\t")
summary(aggregate_CRS_rearranged_mat_gDMRs_sorted_counted1)
summary(aggregate_CRS_rearranged_mat_gDMRs_sorted_counted1$Control_Placenta)
summary(aggregate_CRS_rearranged_mat_gDMRs_sorted_counted1$KHDC3L_BiCHM)
summary(aggregate_CRS_rearranged_mat_gDMRs_sorted_counted1$AnCHM)
wilcox.test(aggregate_CRS_rearranged_mat_gDMRs_sorted_counted1$Control_Placenta, aggregate_CRS_rearranged_mat_gDMRs_sorted_counted1$KHDC3L_BiCHM)
wilcox.test(aggregate_CRS_rearranged_mat_gDMRs_sorted_counted1$Control_Placenta, aggregate_CRS_rearranged_mat_gDMRs_sorted_counted1$AnCHM)
aggregate_CRS_rearranged_mat_gDMRs_sorted_counted2 <- stack(data.frame(aggregate_CRS_rearranged_mat_gDMRs_sorted_counted1))
head(aggregate_CRS_rearranged_mat_gDMRs_sorted_counted2)
colnames(aggregate_CRS_rearranged_mat_gDMRs_sorted_counted2) <- c("Methylation", "Datasets")
head(aggregate_CRS_rearranged_mat_gDMRs_sorted_counted2)
ggboxplot(aggregate_CRS_rearranged_mat_gDMRs_sorted_counted2, x = "Datasets", y = "Methylation",
          title = "AvgMeth_oocyte_placenta_gDMRs", ylab = "Methylation",
          palette =c("#4292A6", "#AEE2EF", "#E6FBFF"),fill = "Datasets", bxp.errorbar = TRUE, ggtheme = theme_classic()) + ylim(0,1)
#ggsave("Boxplot_CRS_merge_hg38_updated_Monk_oocyte_Hanna_placenta_gDMRs_filt.svg", width=12*1.25, height=15*1.25, units="cm", dpi=96)

#Meth CGIs
CRS_rearranged_Meth_CGIs <- read.table("CRS_merge_hg38_hum_placenta_Meth_CGIs_rearranged.txt", header = FALSE)
head(CRS_rearranged_Meth_CGIs)
aggregate_CRS_rearranged_Meth_CGIs = aggregate(CRS_rearranged_Meth_CGIs[,2:4],by=list(CRS_rearranged_Meth_CGIs$V1), mean)
head(aggregate_CRS_rearranged_Meth_CGIs)
aggregate_CRS_rearranged_Meth_CGIs_sorted = aggregate_CRS_rearranged_Meth_CGIs[order(aggregate_CRS_rearranged_Meth_CGIs$Group.1),]
head(aggregate_CRS_rearranged_Meth_CGIs_sorted)
rownames(aggregate_CRS_rearranged_Meth_CGIs_sorted)
rownames(aggregate_CRS_rearranged_Meth_CGIs_sorted)=aggregate_CRS_rearranged_Meth_CGIs_sorted[,1]
rownames(aggregate_CRS_rearranged_Meth_CGIs_sorted)
colnames(aggregate_CRS_rearranged_Meth_CGIs_sorted)
aggregate_CRS_rearranged_Meth_CGIs_sorted = aggregate_CRS_rearranged_Meth_CGIs_sorted[,-1]
aggregate_CRS_rearranged_Meth_CGIs_sorted = as.matrix(aggregate_CRS_rearranged_Meth_CGIs_sorted)
head(aggregate_CRS_rearranged_Meth_CGIs_sorted)
dim(aggregate_CRS_rearranged_Meth_CGIs_sorted)
Count_CRS_rearranged_Meth_CGIs <- count(CRS_rearranged_Meth_CGIs, "V1")
aggregate_CRS_rearranged_Meth_CGIs_sorted_counted <- cbind(aggregate_CRS_rearranged_Meth_CGIs_sorted, Count_CRS_rearranged_Meth_CGIs)
head(aggregate_CRS_rearranged_Meth_CGIs_sorted_counted)
dim(aggregate_CRS_rearranged_Meth_CGIs_sorted_counted)
aggregate_CRS_rearranged_Meth_CGIs_sorted_counted <- aggregate_CRS_rearranged_Meth_CGIs_sorted_counted[which(aggregate_CRS_rearranged_Meth_CGIs_sorted_counted$freq >= 3),]
#write.table(aggregate_CRS_rearranged_Meth_CGIs_sorted_counted, "aggregate_CRS_rearranged_Meth_CGIs_sorted_counted.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE)
aggregate_CRS_rearranged_Meth_CGIs_sorted_counted1 <- aggregate_CRS_rearranged_Meth_CGIs_sorted_counted[,1:3]
colnames(aggregate_CRS_rearranged_Meth_CGIs_sorted_counted1) <- c("Control_Placenta","KHDC3L_BiCHM","AnCHM")
head(aggregate_CRS_rearranged_Meth_CGIs_sorted_counted1)
dim(aggregate_CRS_rearranged_Meth_CGIs_sorted_counted1)
write.table(aggregate_CRS_rearranged_Meth_CGIs_sorted_counted1, "aggregate_CRS_rearranged_Meth_CGIs_sorted_counted1.txt" , quote = FALSE, append = FALSE, sep = "\t")
summary(aggregate_CRS_rearranged_Meth_CGIs_sorted_counted1)
summary(aggregate_CRS_rearranged_Meth_CGIs_sorted_counted1$Control_Placenta)
summary(aggregate_CRS_rearranged_Meth_CGIs_sorted_counted1$KHDC3L_BiCHM)
summary(aggregate_CRS_rearranged_Meth_CGIs_sorted_counted1$AnCHM)
wilcox.test(aggregate_CRS_rearranged_Meth_CGIs_sorted_counted1$Control_Placenta, aggregate_CRS_rearranged_Meth_CGIs_sorted_counted1$KHDC3L_BiCHM)
wilcox.test(aggregate_CRS_rearranged_Meth_CGIs_sorted_counted1$Control_Placenta, aggregate_CRS_rearranged_Meth_CGIs_sorted_counted1$AnCHM)
#Take only Control_Placenta and KHDC3L_BiCHM
aggregate_CRS_rearranged_Meth_CGIs_sorted_counted1 <- aggregate_CRS_rearranged_Meth_CGIs_sorted_counted1[,1:2]
aggregate_CRS_rearranged_Meth_CGIs_sorted_counted2 <- stack(data.frame(aggregate_CRS_rearranged_Meth_CGIs_sorted_counted1))
head(aggregate_CRS_rearranged_Meth_CGIs_sorted_counted2)
colnames(aggregate_CRS_rearranged_Meth_CGIs_sorted_counted2) <- c("Methylation", "Datasets")
head(aggregate_CRS_rearranged_Meth_CGIs_sorted_counted2)
ggboxplot(aggregate_CRS_rearranged_Meth_CGIs_sorted_counted2, x = "Datasets", y = "Methylation",
          title = "CR_meth_CGIs", ylab = "Methylation",
          palette =c("#4292A6", "#AEE2EF"),fill = "Datasets", bxp.errorbar = TRUE, ggtheme = theme_classic()) + ylim(0,1)
#ggsave("Boxplot_CR_meth_CGIs_2_filt.svg", width=12*1.25, height=15*1.25, units="cm", dpi=96)

#Unmeth CGIs
CRS_rearranged_Unmeth_CGIs <- read.table("CRS_merge_hg38_hum_placenta_Unmeth_CGIs_rearranged.txt", header = FALSE)
head(CRS_rearranged_Unmeth_CGIs)
aggregate_CRS_rearranged_Unmeth_CGIs = aggregate(CRS_rearranged_Unmeth_CGIs[,2:4],by=list(CRS_rearranged_Unmeth_CGIs$V1), mean)
head(aggregate_CRS_rearranged_Unmeth_CGIs)
aggregate_CRS_rearranged_Unmeth_CGIs_sorted = aggregate_CRS_rearranged_Unmeth_CGIs[order(aggregate_CRS_rearranged_Unmeth_CGIs$Group.1),]
head(aggregate_CRS_rearranged_Unmeth_CGIs_sorted)
rownames(aggregate_CRS_rearranged_Unmeth_CGIs_sorted)
rownames(aggregate_CRS_rearranged_Unmeth_CGIs_sorted)=aggregate_CRS_rearranged_Unmeth_CGIs_sorted[,1]
rownames(aggregate_CRS_rearranged_Unmeth_CGIs_sorted)
colnames(aggregate_CRS_rearranged_Unmeth_CGIs_sorted)
aggregate_CRS_rearranged_Unmeth_CGIs_sorted = aggregate_CRS_rearranged_Unmeth_CGIs_sorted[,-1]
aggregate_CRS_rearranged_Unmeth_CGIs_sorted = as.matrix(aggregate_CRS_rearranged_Unmeth_CGIs_sorted)
head(aggregate_CRS_rearranged_Unmeth_CGIs_sorted)
dim(aggregate_CRS_rearranged_Unmeth_CGIs_sorted)
Count_CRS_rearranged_Unmeth_CGIs <- count(CRS_rearranged_Unmeth_CGIs, "V1")
aggregate_CRS_rearranged_Unmeth_CGIs_sorted_counted <- cbind(aggregate_CRS_rearranged_Unmeth_CGIs_sorted, Count_CRS_rearranged_Unmeth_CGIs)
head(aggregate_CRS_rearranged_Unmeth_CGIs_sorted_counted)
dim(aggregate_CRS_rearranged_Unmeth_CGIs_sorted_counted)
aggregate_CRS_rearranged_Unmeth_CGIs_sorted_counted <- aggregate_CRS_rearranged_Unmeth_CGIs_sorted_counted[which(aggregate_CRS_rearranged_Unmeth_CGIs_sorted_counted$freq >= 3),]
#write.table(aggregate_CRS_rearranged_Unmeth_CGIs_sorted_counted, "aggregate_CRS_rearranged_Unmeth_CGIs_sorted_counted.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE)
aggregate_CRS_rearranged_Unmeth_CGIs_sorted_counted1 <- aggregate_CRS_rearranged_Unmeth_CGIs_sorted_counted[,1:3]
colnames(aggregate_CRS_rearranged_Unmeth_CGIs_sorted_counted1) <- c("Control_Placenta","KHDC3L_BiCHM","AnCHM")
head(aggregate_CRS_rearranged_Unmeth_CGIs_sorted_counted1)
dim(aggregate_CRS_rearranged_Unmeth_CGIs_sorted_counted1)
write.table(aggregate_CRS_rearranged_Unmeth_CGIs_sorted_counted1, "aggregate_CRS_rearranged_Unmeth_CGIs_sorted_counted1.txt" , quote = FALSE, append = FALSE, sep = "\t")
summary(aggregate_CRS_rearranged_Unmeth_CGIs_sorted_counted1)
summary(aggregate_CRS_rearranged_Unmeth_CGIs_sorted_counted1$Control_Placenta)
summary(aggregate_CRS_rearranged_Unmeth_CGIs_sorted_counted1$KHDC3L_BiCHM)
summary(aggregate_CRS_rearranged_Unmeth_CGIs_sorted_counted1$AnCHM)
wilcox.test(aggregate_CRS_rearranged_Unmeth_CGIs_sorted_counted1$Control_Placenta, aggregate_CRS_rearranged_Unmeth_CGIs_sorted_counted1$KHDC3L_BiCHM)
wilcox.test(aggregate_CRS_rearranged_Unmeth_CGIs_sorted_counted1$Control_Placenta, aggregate_CRS_rearranged_Unmeth_CGIs_sorted_counted1$AnCHM)
#Take only Control_Placenta and KHDC3L_BiCHM
aggregate_CRS_rearranged_Unmeth_CGIs_sorted_counted1 <- aggregate_CRS_rearranged_Unmeth_CGIs_sorted_counted1[,1:2]
aggregate_CRS_rearranged_Unmeth_CGIs_sorted_counted2 <- stack(data.frame(aggregate_CRS_rearranged_Unmeth_CGIs_sorted_counted1))
head(aggregate_CRS_rearranged_Unmeth_CGIs_sorted_counted2)
colnames(aggregate_CRS_rearranged_Unmeth_CGIs_sorted_counted2) <- c("Methylation", "Datasets")
head(aggregate_CRS_rearranged_Unmeth_CGIs_sorted_counted2)
ggboxplot(aggregate_CRS_rearranged_Unmeth_CGIs_sorted_counted2, x = "Datasets", y = "Methylation",
          title = "CR_Unmeth_CGIs", ylab = "Methylation",
          palette =c("#4292A6", "#AEE2EF"),fill = "Datasets", bxp.errorbar = TRUE, ggtheme = theme_classic()) + ylim(0,1)
#ggsave("Boxplot_CR_Unmeth_CGIs_2_filt.svg", width=12*1.25, height=15*1.25, units="cm", dpi=96)

############################        Prepare Non-CGIs    ##############################################
CRS_rearranged_NonCGIs <- read.table("CRS_merge_hg38_20kB-windows_no-CGIs_rearranged.txt", header = FALSE)
head(CRS_rearranged_NonCGIs)
aggregate_CRS_rearranged_NonCGIs = aggregate(CRS_rearranged_NonCGIs[,2:4],by=list(CRS_rearranged_NonCGIs$V1), mean)
head(aggregate_CRS_rearranged_NonCGIs)
aggregate_CRS_rearranged_NonCGIs_sorted = aggregate_CRS_rearranged_NonCGIs[order(aggregate_CRS_rearranged_NonCGIs$Group.1),]
head(aggregate_CRS_rearranged_NonCGIs_sorted)
rownames(aggregate_CRS_rearranged_NonCGIs_sorted)
rownames(aggregate_CRS_rearranged_NonCGIs_sorted)=aggregate_CRS_rearranged_NonCGIs_sorted[,1]
rownames(aggregate_CRS_rearranged_NonCGIs_sorted)
colnames(aggregate_CRS_rearranged_NonCGIs_sorted)
aggregate_CRS_rearranged_NonCGIs_sorted = aggregate_CRS_rearranged_NonCGIs_sorted[,-1]
aggregate_CRS_rearranged_NonCGIs_sorted = as.matrix(aggregate_CRS_rearranged_NonCGIs_sorted)
head(aggregate_CRS_rearranged_NonCGIs_sorted)
dim(aggregate_CRS_rearranged_NonCGIs_sorted)
Count_CRS_rearranged_NonCGIs <- count(CRS_rearranged_NonCGIs, "V1")
aggregate_CRS_rearranged_NonCGIs_sorted_counted <- cbind(aggregate_CRS_rearranged_NonCGIs_sorted, Count_CRS_rearranged_NonCGIs)
head(aggregate_CRS_rearranged_NonCGIs_sorted_counted)
dim(aggregate_CRS_rearranged_NonCGIs_sorted_counted)
aggregate_CRS_rearranged_NonCGIs_sorted_counted <- aggregate_CRS_rearranged_NonCGIs_sorted_counted[which(aggregate_CRS_rearranged_NonCGIs_sorted_counted$freq >= 3),]
#write.table(aggregate_CRS_rearranged_NonCGIs_sorted_counted, "aggregate_CRS_rearranged_NonCGIs_sorted_counted.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE)
aggregate_CRS_rearranged_NonCGIs_sorted_counted1 <- aggregate_CRS_rearranged_NonCGIs_sorted_counted[,1:3]
colnames(aggregate_CRS_rearranged_NonCGIs_sorted_counted1) <- c("Control_Placenta","KHDC3L_BiCHM","AnCHM")
head(aggregate_CRS_rearranged_NonCGIs_sorted_counted1)
dim(aggregate_CRS_rearranged_NonCGIs_sorted_counted1)
write.table(aggregate_CRS_rearranged_NonCGIs_sorted_counted1, "aggregate_CRS_rearranged_NonCGIs_sorted_counted1.txt" , quote = FALSE, append = FALSE, sep = "\t")
summary(aggregate_CRS_rearranged_NonCGIs_sorted_counted1)
summary(aggregate_CRS_rearranged_NonCGIs_sorted_counted1$Control_Placenta)
summary(aggregate_CRS_rearranged_NonCGIs_sorted_counted1$KHDC3L_BiCHM)
summary(aggregate_CRS_rearranged_NonCGIs_sorted_counted1$AnCHM)
wilcox.test(aggregate_CRS_rearranged_NonCGIs_sorted_counted1$Control_Placenta, aggregate_CRS_rearranged_NonCGIs_sorted_counted1$KHDC3L_BiCHM)
wilcox.test(aggregate_CRS_rearranged_NonCGIs_sorted_counted1$Control_Placenta, aggregate_CRS_rearranged_NonCGIs_sorted_counted1$AnCHM)
#Take only Control_Placenta and KHDC3L_BiCHM
aggregate_CRS_rearranged_NonCGIs_sorted_counted1 <- aggregate_CRS_rearranged_NonCGIs_sorted_counted1[,1:2]
aggregate_CRS_rearranged_NonCGIs_sorted_counted2 <- stack(data.frame(aggregate_CRS_rearranged_NonCGIs_sorted_counted1))
head(aggregate_CRS_rearranged_NonCGIs_sorted_counted2)
colnames(aggregate_CRS_rearranged_NonCGIs_sorted_counted2) <- c("Methylation", "Datasets")
head(aggregate_CRS_rearranged_NonCGIs_sorted_counted2)
ggboxplot(aggregate_CRS_rearranged_NonCGIs_sorted_counted2, x = "Datasets", y = "Methylation",
          title = "CR_NonCGIs", ylab = "Methylation",
          palette =c("#4292A6", "#AEE2EF"),fill = "Datasets", bxp.errorbar = TRUE, ggtheme = theme_classic()) + ylim(0,1)
#ggsave("Boxplot_CR_NonCGIs_2_filt.svg", width=12*1.25, height=15*1.25, units="cm", dpi=96)

############################        Prepare Non-CGIs:Specific   ##############################################
CRS_rearranged_NonCGIs_sp <- read.table("CRS_merge_hg38_20kB-windows_no-CGIs_specific_meth_rearranged.txt", header = FALSE)
head(CRS_rearranged_NonCGIs_sp)
aggregate_CRS_rearranged_NonCGIs_sp = aggregate(CRS_rearranged_NonCGIs_sp[,2:4],by=list(CRS_rearranged_NonCGIs_sp$V1), mean)
head(aggregate_CRS_rearranged_NonCGIs_sp)
aggregate_CRS_rearranged_NonCGIs_sp_sorted = aggregate_CRS_rearranged_NonCGIs_sp[order(aggregate_CRS_rearranged_NonCGIs_sp$Group.1),]
head(aggregate_CRS_rearranged_NonCGIs_sp_sorted)
rownames(aggregate_CRS_rearranged_NonCGIs_sp_sorted)
rownames(aggregate_CRS_rearranged_NonCGIs_sp_sorted)=aggregate_CRS_rearranged_NonCGIs_sp_sorted[,1]
rownames(aggregate_CRS_rearranged_NonCGIs_sp_sorted)
colnames(aggregate_CRS_rearranged_NonCGIs_sp_sorted)
aggregate_CRS_rearranged_NonCGIs_sp_sorted = aggregate_CRS_rearranged_NonCGIs_sp_sorted[,-1]
aggregate_CRS_rearranged_NonCGIs_sp_sorted = as.matrix(aggregate_CRS_rearranged_NonCGIs_sp_sorted)
head(aggregate_CRS_rearranged_NonCGIs_sp_sorted)
dim(aggregate_CRS_rearranged_NonCGIs_sp_sorted)
Count_CRS_rearranged_NonCGIs_sp <- count(CRS_rearranged_NonCGIs_sp, "V1")
aggregate_CRS_rearranged_NonCGIs_sp_sorted_counted <- cbind(aggregate_CRS_rearranged_NonCGIs_sp_sorted, Count_CRS_rearranged_NonCGIs_sp)
head(aggregate_CRS_rearranged_NonCGIs_sp_sorted_counted)
dim(aggregate_CRS_rearranged_NonCGIs_sp_sorted_counted)
aggregate_CRS_rearranged_NonCGIs_sp_sorted_counted <- aggregate_CRS_rearranged_NonCGIs_sp_sorted_counted[which(aggregate_CRS_rearranged_NonCGIs_sp_sorted_counted$freq >= 3),]
#write.table(aggregate_CRS_rearranged_NonCGIs_sp_sorted_counted, "aggregate_CRS_rearranged_NonCGIs_sp_sorted_counted.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE)
aggregate_CRS_rearranged_NonCGIs_sp_sorted_counted1 <- aggregate_CRS_rearranged_NonCGIs_sp_sorted_counted[,1:3]
colnames(aggregate_CRS_rearranged_NonCGIs_sp_sorted_counted1) <- c("Control_Placenta","KHDC3L_BiCHM","AnCHM")
head(aggregate_CRS_rearranged_NonCGIs_sp_sorted_counted1)
dim(aggregate_CRS_rearranged_NonCGIs_sp_sorted_counted1)
write.table(aggregate_CRS_rearranged_NonCGIs_sp_sorted_counted1, "aggregate_CRS_rearranged_NonCGIs_sp_sorted_counted1.txt" , quote = FALSE, append = FALSE, sep = "\t")
summary(aggregate_CRS_rearranged_NonCGIs_sp_sorted_counted1)
summary(aggregate_CRS_rearranged_NonCGIs_sp_sorted_counted1$Control_Placenta)
summary(aggregate_CRS_rearranged_NonCGIs_sp_sorted_counted1$KHDC3L_BiCHM)
summary(aggregate_CRS_rearranged_NonCGIs_sp_sorted_counted1$AnCHM)
wilcox.test(aggregate_CRS_rearranged_NonCGIs_sp_sorted_counted1$Control_Placenta, aggregate_CRS_rearranged_NonCGIs_sp_sorted_counted1$KHDC3L_BiCHM)
wilcox.test(aggregate_CRS_rearranged_NonCGIs_sp_sorted_counted1$Control_Placenta, aggregate_CRS_rearranged_NonCGIs_sp_sorted_counted1$AnCHM)
#Take only Control_Placenta and KHDC3L_BiCHM
aggregate_CRS_rearranged_NonCGIs_sp_sorted_counted1 <- aggregate_CRS_rearranged_NonCGIs_sp_sorted_counted1[,1:2]
aggregate_CRS_rearranged_NonCGIs_sp_sorted_counted2 <- stack(data.frame(aggregate_CRS_rearranged_NonCGIs_sp_sorted_counted1))
head(aggregate_CRS_rearranged_NonCGIs_sp_sorted_counted2)
colnames(aggregate_CRS_rearranged_NonCGIs_sp_sorted_counted2) <- c("Methylation", "Datasets")
head(aggregate_CRS_rearranged_NonCGIs_sp_sorted_counted2)
ggboxplot(aggregate_CRS_rearranged_NonCGIs_sp_sorted_counted2, x = "Datasets", y = "Methylation",
          title = "CR_NonCGIs_sp", ylab = "Methylation",
          palette =c("#4292A6", "#AEE2EF"),fill = "Datasets", bxp.errorbar = TRUE, ggtheme = theme_classic()) + ylim(0,1)
#ggsave("Boxplot_CR_NonCGIs_sp_2_filt.svg", width=12*1.25, height=15*1.25, units="cm", dpi=96)


############################        Prepare PmCGIs   ##############################################
CRS_rearranged_PmCGIs <- read.table("CRS_merge_hg38_PmCGI_meth_rearranged.txt", header = FALSE)
head(CRS_rearranged_PmCGIs)
aggregate_CRS_rearranged_PmCGIs = aggregate(CRS_rearranged_PmCGIs[,2:4],by=list(CRS_rearranged_PmCGIs$V1), mean)
head(aggregate_CRS_rearranged_PmCGIs)
aggregate_CRS_rearranged_PmCGIs_sorted = aggregate_CRS_rearranged_PmCGIs[order(aggregate_CRS_rearranged_PmCGIs$Group.1),]
head(aggregate_CRS_rearranged_PmCGIs_sorted)
rownames(aggregate_CRS_rearranged_PmCGIs_sorted)
rownames(aggregate_CRS_rearranged_PmCGIs_sorted)=aggregate_CRS_rearranged_PmCGIs_sorted[,1]
rownames(aggregate_CRS_rearranged_PmCGIs_sorted)
colnames(aggregate_CRS_rearranged_PmCGIs_sorted)
aggregate_CRS_rearranged_PmCGIs_sorted = aggregate_CRS_rearranged_PmCGIs_sorted[,-1]
aggregate_CRS_rearranged_PmCGIs_sorted = as.matrix(aggregate_CRS_rearranged_PmCGIs_sorted)
head(aggregate_CRS_rearranged_PmCGIs_sorted)
dim(aggregate_CRS_rearranged_PmCGIs_sorted)
Count_CRS_rearranged_PmCGIs <- count(CRS_rearranged_PmCGIs, "V1")
aggregate_CRS_rearranged_PmCGIs_sorted_counted <- cbind(aggregate_CRS_rearranged_PmCGIs_sorted, Count_CRS_rearranged_PmCGIs)
head(aggregate_CRS_rearranged_PmCGIs_sorted_counted)
dim(aggregate_CRS_rearranged_PmCGIs_sorted_counted)
aggregate_CRS_rearranged_PmCGIs_sorted_counted <- aggregate_CRS_rearranged_PmCGIs_sorted_counted[which(aggregate_CRS_rearranged_PmCGIs_sorted_counted$freq >= 3),]
#write.table(aggregate_CRS_rearranged_PmCGIs_sorted_counted, "aggregate_CRS_rearranged_PmCGIs_sorted_counted.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE)
aggregate_CRS_rearranged_PmCGIs_sorted_counted1 <- aggregate_CRS_rearranged_PmCGIs_sorted_counted[,1:3]
colnames(aggregate_CRS_rearranged_PmCGIs_sorted_counted1) <- c("Control_Placenta","KHDC3L_BiCHM","AnCHM")
head(aggregate_CRS_rearranged_PmCGIs_sorted_counted1)
dim(aggregate_CRS_rearranged_PmCGIs_sorted_counted1)
write.table(aggregate_CRS_rearranged_PmCGIs_sorted_counted1, "aggregate_CRS_rearranged_PmCGIs_sorted_counted1.txt" , quote = FALSE, append = FALSE, sep = "\t")
summary(aggregate_CRS_rearranged_PmCGIs_sorted_counted1)
summary(aggregate_CRS_rearranged_PmCGIs_sorted_counted1$Control_Placenta)
summary(aggregate_CRS_rearranged_PmCGIs_sorted_counted1$KHDC3L_BiCHM)
summary(aggregate_CRS_rearranged_PmCGIs_sorted_counted1$AnCHM)
wilcox.test(aggregate_CRS_rearranged_PmCGIs_sorted_counted1$Control_Placenta, aggregate_CRS_rearranged_PmCGIs_sorted_counted1$KHDC3L_BiCHM)
wilcox.test(aggregate_CRS_rearranged_PmCGIs_sorted_counted1$Control_Placenta, aggregate_CRS_rearranged_PmCGIs_sorted_counted1$AnCHM)
#Take only Control_Placenta and KHDC3L_BiCHM
aggregate_CRS_rearranged_PmCGIs_sorted_counted1 <- aggregate_CRS_rearranged_PmCGIs_sorted_counted1[,1:2]
aggregate_CRS_rearranged_PmCGIs_sorted_counted2 <- stack(data.frame(aggregate_CRS_rearranged_PmCGIs_sorted_counted1))
head(aggregate_CRS_rearranged_PmCGIs_sorted_counted2)
colnames(aggregate_CRS_rearranged_PmCGIs_sorted_counted2) <- c("Methylation", "Datasets")
head(aggregate_CRS_rearranged_PmCGIs_sorted_counted2)
ggboxplot(aggregate_CRS_rearranged_PmCGIs_sorted_counted2, x = "Datasets", y = "Methylation",
          title = "CR_PmCGIs", ylab = "Methylation",
          palette =c("#4292A6", "#AEE2EF"),fill = "Datasets", bxp.errorbar = TRUE, ggtheme = theme_classic()) + ylim(0,1)
#ggsave("Boxplot_CR_PmCGIs_2_filt.svg", width=12*1.25, height=15*1.25, units="cm", dpi=96)

############################        Prepare OmUp   ##############################################
CRS_rearranged_OmUp <- read.table("CRS_merge_hg38_OmUp_rearranged.txt", header = FALSE)
head(CRS_rearranged_OmUp)
aggregate_CRS_rearranged_OmUp = aggregate(CRS_rearranged_OmUp[,2:4],by=list(CRS_rearranged_OmUp$V1), mean)
head(aggregate_CRS_rearranged_OmUp)
aggregate_CRS_rearranged_OmUp_sorted = aggregate_CRS_rearranged_OmUp[order(aggregate_CRS_rearranged_OmUp$Group.1),]
head(aggregate_CRS_rearranged_OmUp_sorted)
rownames(aggregate_CRS_rearranged_OmUp_sorted)
rownames(aggregate_CRS_rearranged_OmUp_sorted)=aggregate_CRS_rearranged_OmUp_sorted[,1]
rownames(aggregate_CRS_rearranged_OmUp_sorted)
colnames(aggregate_CRS_rearranged_OmUp_sorted)
aggregate_CRS_rearranged_OmUp_sorted = aggregate_CRS_rearranged_OmUp_sorted[,-1]
aggregate_CRS_rearranged_OmUp_sorted = as.matrix(aggregate_CRS_rearranged_OmUp_sorted)
head(aggregate_CRS_rearranged_OmUp_sorted)
dim(aggregate_CRS_rearranged_OmUp_sorted)
Count_CRS_rearranged_OmUp <- count(CRS_rearranged_OmUp, "V1")
aggregate_CRS_rearranged_OmUp_sorted_counted <- cbind(aggregate_CRS_rearranged_OmUp_sorted, Count_CRS_rearranged_OmUp)
head(aggregate_CRS_rearranged_OmUp_sorted_counted)
dim(aggregate_CRS_rearranged_OmUp_sorted_counted)
aggregate_CRS_rearranged_OmUp_sorted_counted <- aggregate_CRS_rearranged_OmUp_sorted_counted[which(aggregate_CRS_rearranged_OmUp_sorted_counted$freq >= 3),]
#write.table(aggregate_CRS_rearranged_OmUp_sorted_counted, "aggregate_CRS_rearranged_OmUp_sorted_counted.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE)
aggregate_CRS_rearranged_OmUp_sorted_counted1 <- aggregate_CRS_rearranged_OmUp_sorted_counted[,1:3]
colnames(aggregate_CRS_rearranged_OmUp_sorted_counted1) <- c("Control_Placenta","KHDC3L_BiCHM","AnCHM")
head(aggregate_CRS_rearranged_OmUp_sorted_counted1)
dim(aggregate_CRS_rearranged_OmUp_sorted_counted1)
write.table(aggregate_CRS_rearranged_OmUp_sorted_counted1, "aggregate_CRS_rearranged_OmUp_sorted_counted1.txt" , quote = FALSE, append = FALSE, sep = "\t")
summary(aggregate_CRS_rearranged_OmUp_sorted_counted1)
summary(aggregate_CRS_rearranged_OmUp_sorted_counted1$Control_Placenta)
summary(aggregate_CRS_rearranged_OmUp_sorted_counted1$KHDC3L_BiCHM)
summary(aggregate_CRS_rearranged_OmUp_sorted_counted1$AnCHM)
wilcox.test(aggregate_CRS_rearranged_OmUp_sorted_counted1$Control_Placenta, aggregate_CRS_rearranged_OmUp_sorted_counted1$KHDC3L_BiCHM)
wilcox.test(aggregate_CRS_rearranged_OmUp_sorted_counted1$Control_Placenta, aggregate_CRS_rearranged_OmUp_sorted_counted1$AnCHM)
#Take only Control_Placenta and KHDC3L_BiCHM
aggregate_CRS_rearranged_OmUp_sorted_counted1 <- aggregate_CRS_rearranged_OmUp_sorted_counted1[,1:2]
aggregate_CRS_rearranged_OmUp_sorted_counted2 <- stack(data.frame(aggregate_CRS_rearranged_OmUp_sorted_counted1))
head(aggregate_CRS_rearranged_OmUp_sorted_counted2)
colnames(aggregate_CRS_rearranged_OmUp_sorted_counted2) <- c("Methylation", "Datasets")
head(aggregate_CRS_rearranged_OmUp_sorted_counted2)
ggboxplot(aggregate_CRS_rearranged_OmUp_sorted_counted2, x = "Datasets", y = "Methylation",
          title = "CR_OmUp", ylab = "Methylation",
          palette =c("#4292A6", "#AEE2EF"),fill = "Datasets", bxp.errorbar = TRUE, ggtheme = theme_classic()) + ylim(0,1)
#ggsave("Boxplot_CR_OmUp_2_filt.svg", width=12*1.25, height=15*1.25, units="cm", dpi=96)

#Statistics
summary(aggregate_CRS_rearranged_mat_gDMRs_sorted_counted1$Control_Placenta)
summary(aggregate_CRS_rearranged_mat_gDMRs_sorted_counted1$KHDC3L_BiCHM)
summary(aggregate_CRS_rearranged_mat_gDMRs_sorted_counted1$AnCHM)
summary(aggregate_CRS_rearranged_Meth_CGIs_sorted_counted1$Control_Placenta)
summary(aggregate_CRS_rearranged_Meth_CGIs_sorted_counted1$KHDC3L_BiCHM)
summary(aggregate_CRS_rearranged_Meth_CGIs_sorted_counted1$AnCHM)
summary(aggregate_CRS_rearranged_Unmeth_CGIs_sorted_counted1$Control_Placenta)
summary(aggregate_CRS_rearranged_Unmeth_CGIs_sorted_counted1$KHDC3L_BiCHM)
summary(aggregate_CRS_rearranged_Unmeth_CGIs_sorted_counted1$AnCHM)
summary(aggregate_CRS_rearranged_NonCGIs_sorted_counted1$Control_Placenta)
summary(aggregate_CRS_rearranged_NonCGIs_sorted_counted1$KHDC3L_BiCHM)
summary(aggregate_CRS_rearranged_NonCGIs_sorted_counted1$AnCHM)
summary(aggregate_CRS_rearranged_NonCGIs_sp_sorted_counted1$Control_Placenta)
summary(aggregate_CRS_rearranged_NonCGIs_sp_sorted_counted1$KHDC3L_BiCHM)
summary(aggregate_CRS_rearranged_NonCGIs_sp_sorted_counted1$AnCHM)
summary(aggregate_CRS_rearranged_PmCGIs_sorted_counted1$Control_Placenta)
summary(aggregate_CRS_rearranged_PmCGIs_sorted_counted1$KHDC3L_BiCHM)
summary(aggregate_CRS_rearranged_PmCGIs_sorted_counted1$AnCHM)
summary(aggregate_CRS_rearranged_OmUp_sorted_counted1$Control_Placenta)
summary(aggregate_CRS_rearranged_OmUp_sorted_counted1$KHDC3L_BiCHM)
summary(aggregate_CRS_rearranged_OmUp_sorted_counted1$AnCHM)

summary(aggregate_CRS_rearranged_mat_gDMRs_sorted_counted1)
summary(aggregate_CRS_rearranged_Meth_CGIs_sorted_counted1)
summary(aggregate_CRS_rearranged_Unmeth_CGIs_sorted_counted1)
summary(aggregate_CRS_rearranged_NonCGIs_sorted_counted1)
summary(aggregate_CRS_rearranged_NonCGIs_sp_sorted_counted1)
summary(aggregate_CRS_rearranged_PmCGIs_sorted_counted1)
summary(aggregate_CRS_rearranged_OmUp_sorted_counted1)


basicStats(aggregate_CRS_rearranged_mat_gDMRs_sorted_counted1$Control_Placenta)[c("1. Quartile","Median", "3. Quartile"),]
basicStats(aggregate_CRS_rearranged_mat_gDMRs_sorted_counted1$KHDC3L_BiCHM)[c("1. Quartile","Median", "3. Quartile"),]
basicStats(aggregate_CRS_rearranged_mat_gDMRs_sorted_counted1$AnCHM)[c("1. Quartile","Median", "3. Quartile"),]

basicStats(aggregate_CRS_rearranged_Meth_CGIs_sorted_counted1$Control_Placenta)[c("1. Quartile","Median", "3. Quartile"),]
basicStats(aggregate_CRS_rearranged_Meth_CGIs_sorted_counted1$KHDC3L_BiCHM)[c("1. Quartile","Median", "3. Quartile"),]

basicStats(aggregate_CRS_rearranged_Unmeth_CGIs_sorted_counted1$Control_Placenta)[c("1. Quartile","Median", "3. Quartile"),]
basicStats(aggregate_CRS_rearranged_Unmeth_CGIs_sorted_counted1$KHDC3L_BiCHM)[c("1. Quartile","Median", "3. Quartile"),]

basicStats(aggregate_CRS_rearranged_NonCGIs_sorted_counted1$Control_Placenta)[c("1. Quartile","Median", "3. Quartile"),]
basicStats(aggregate_CRS_rearranged_NonCGIs_sorted_counted1$KHDC3L_BiCHM)[c("1. Quartile","Median", "3. Quartile"),]

basicStats(aggregate_CRS_rearranged_NonCGIs_sp_sorted_counted1$Control_Placenta)[c("1. Quartile","Median", "3. Quartile"),]
basicStats(aggregate_CRS_rearranged_NonCGIs_sp_sorted_counted1$KHDC3L_BiCHM)[c("1. Quartile","Median", "3. Quartile"),]

basicStats(aggregate_CRS_rearranged_PmCGIs_sorted_counted1$Control_Placenta)[c("1. Quartile","Median", "3. Quartile"),]
basicStats(aggregate_CRS_rearranged_PmCGIs_sorted_counted1$KHDC3L_BiCHM)[c("1. Quartile","Median", "3. Quartile"),]

basicStats(aggregate_CRS_rearranged_OmUp_sorted_counted1$Control_Placenta)[c("1. Quartile","Median", "3. Quartile"),]
basicStats(aggregate_CRS_rearranged_OmUp_sorted_counted1$KHDC3L_BiCHM)[c("1. Quartile","Median", "3. Quartile"),]


dim(aggregate_CRS_rearranged_mat_gDMRs_sorted_counted1)
dim(aggregate_CRS_rearranged_Meth_CGIs_sorted_counted1)
dim(aggregate_CRS_rearranged_Unmeth_CGIs_sorted_counted1)
dim(aggregate_CRS_rearranged_NonCGIs_sorted_counted1)
dim(aggregate_CRS_rearranged_NonCGIs_sp_sorted_counted1)
dim(aggregate_CRS_rearranged_PmCGIs_sorted_counted1)
dim(aggregate_CRS_rearranged_OmUp_sorted_counted1)


wilcox.test(aggregate_CRS_rearranged_mat_gDMRs_sorted_counted1$Control_Placenta, aggregate_CRS_rearranged_mat_gDMRs_sorted_counted1$KHDC3L_BiCHM)$p.value
wilcox.test(aggregate_CRS_rearranged_mat_gDMRs_sorted_counted1$Control_Placenta, aggregate_CRS_rearranged_mat_gDMRs_sorted_counted1$AnCHM)$p.value
wilcox.test(aggregate_CRS_rearranged_Meth_CGIs_sorted_counted1$Control_Placenta, aggregate_CRS_rearranged_Meth_CGIs_sorted_counted1$KHDC3L_BiCHM)$p.value
wilcox.test(aggregate_CRS_rearranged_Unmeth_CGIs_sorted_counted1$Control_Placenta, aggregate_CRS_rearranged_Unmeth_CGIs_sorted_counted1$KHDC3L_BiCHM)$p.value
wilcox.test(aggregate_CRS_rearranged_NonCGIs_sorted_counted1$Control_Placenta, aggregate_CRS_rearranged_NonCGIs_sorted_counted1$KHDC3L_BiCHM)$p.value
wilcox.test(aggregate_CRS_rearranged_NonCGIs_sp_sorted_counted1$Control_Placenta, aggregate_CRS_rearranged_NonCGIs_sp_sorted_counted1$KHDC3L_BiCHM)$p.value
wilcox.test(aggregate_CRS_rearranged_PmCGIs_sorted_counted1$Control_Placenta, aggregate_CRS_rearranged_PmCGIs_sorted_counted1$KHDC3L_BiCHM)$p.value
wilcox.test(aggregate_CRS_rearranged_OmUp_sorted_counted1$Control_Placenta, aggregate_CRS_rearranged_OmUp_sorted_counted1$KHDC3L_BiCHM)$p.value

wilcox.test(aggregate_CRS_rearranged_mat_gDMRs_sorted_counted1$Control_Placenta, aggregate_CRS_rearranged_mat_gDMRs_sorted_counted1$KHDC3L_BiCHM, paired = TRUE)$p.value
wilcox.test(aggregate_CRS_rearranged_mat_gDMRs_sorted_counted1$Control_Placenta, aggregate_CRS_rearranged_mat_gDMRs_sorted_counted1$AnCHM, paired = TRUE)$p.value
wilcox.test(aggregate_CRS_rearranged_Meth_CGIs_sorted_counted1$Control_Placenta, aggregate_CRS_rearranged_Meth_CGIs_sorted_counted1$KHDC3L_BiCHM, paired = TRUE)$p.value
wilcox.test(aggregate_CRS_rearranged_Unmeth_CGIs_sorted_counted1$Control_Placenta, aggregate_CRS_rearranged_Unmeth_CGIs_sorted_counted1$KHDC3L_BiCHM, paired = TRUE)$p.value
wilcox.test(aggregate_CRS_rearranged_NonCGIs_sorted_counted1$Control_Placenta, aggregate_CRS_rearranged_NonCGIs_sorted_counted1$KHDC3L_BiCHM, paired = TRUE)$p.value


wilcox.test(aggregate_CRS_rearranged_mat_gDMRs_sorted_counted1$Control_Placenta, aggregate_CRS_rearranged_mat_gDMRs_sorted_counted1$KHDC3L_BiCHM)
wilcox.test(aggregate_CRS_rearranged_mat_gDMRs_sorted_counted1$Control_Placenta, aggregate_CRS_rearranged_mat_gDMRs_sorted_counted1$AnCHM)
wilcox.test(aggregate_CRS_rearranged_Meth_CGIs_sorted_counted1$Control_Placenta, aggregate_CRS_rearranged_Meth_CGIs_sorted_counted1$KHDC3L_BiCHM)
wilcox.test(aggregate_CRS_rearranged_Unmeth_CGIs_sorted_counted1$Control_Placenta, aggregate_CRS_rearranged_Unmeth_CGIs_sorted_counted1$KHDC3L_BiCHM)
wilcox.test(aggregate_CRS_rearranged_NonCGIs_sorted_counted1$Control_Placenta, aggregate_CRS_rearranged_NonCGIs_sorted_counted1$KHDC3L_BiCHM)
wilcox.test(aggregate_CRS_rearranged_NonCGIs_sp_sorted_counted1$Control_Placenta, aggregate_CRS_rearranged_NonCGIs_sp_sorted_counted1$KHDC3L_BiCHM)
wilcox.test(aggregate_CRS_rearranged_PmCGIs_sorted_counted1$Control_Placenta, aggregate_CRS_rearranged_PmCGIs_sorted_counted1$KHDC3L_BiCHM)
wilcox.test(aggregate_CRS_rearranged_OmUp_sorted_counted1$Control_Placenta, aggregate_CRS_rearranged_OmUp_sorted_counted1$KHDC3L_BiCHM)

wilcox.test(aggregate_CRS_rearranged_Meth_CGIs_sorted_counted1$Control_Placenta, aggregate_CRS_rearranged_Meth_CGIs_sorted_counted1$AnCHM)
wilcox.test(aggregate_CRS_rearranged_Unmeth_CGIs_sorted_counted1$Control_Placenta, aggregate_CRS_rearranged_Unmeth_CGIs_sorted_counted1$AnCHM)
wilcox.test(aggregate_CRS_rearranged_NonCGIs_sorted_counted1$Control_Placenta, aggregate_CRS_rearranged_NonCGIs_sorted_counted1$AnCHM)
wilcox.test(aggregate_CRS_rearranged_NonCGIs_sp_sorted_counted1$Control_Placenta, aggregate_CRS_rearranged_NonCGIs_sp_sorted_counted1$AnCHM)
wilcox.test(aggregate_CRS_rearranged_PmCGIs_sorted_counted1$Control_Placenta, aggregate_CRS_rearranged_PmCGIs_sorted_counted1$AnCHM)
wilcox.test(aggregate_CRS_rearranged_OmUp_sorted_counted1$Control_Placenta, aggregate_CRS_rearranged_OmUp_sorted_counted1$AnCHM)


#Genes
awk '{print $1"\t"$2"\t"$3}' gencode.v29.annotation.gene.gtf > gencode.v29.annotation.gene_hg38.txt

bedtools intersect -wa -wb -a gencode.v29.annotation.gene_hg38.txt -b /media/ankitv/Archivio2/ankit/Mole/bisulfite-seq/bismark_cov/bedgraph/hg38/array_data_hg38/All/CRS_merge_hg38.txt > CRS_merge_hg38_gene.txt
awk '{print $1"%"$2"%"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' CRS_merge_hg38_gene.txt > CRS_merge_hg38_gene_.txt
awk '{print $1"\t"$5"\t"$6"\t"$7}' CRS_merge_hg38_gene_.txt > CRS_merge_hg38_gene_rearranged.txt

#Genes
CRS_merge_hg38_gene <- read.table("CRS_merge_hg38_gene_rearranged.txt", header = FALSE)
head(CRS_merge_hg38_gene)
aggregate_CRS_merge_hg38_gene = aggregate(CRS_merge_hg38_gene[,2:4],by=list(CRS_merge_hg38_gene$V1), mean)
head(aggregate_CRS_merge_hg38_gene)
aggregate_CRS_merge_hg38_gene_sorted = aggregate_CRS_merge_hg38_gene[order(aggregate_CRS_merge_hg38_gene$Group.1),]
head(aggregate_CRS_merge_hg38_gene_sorted)
rownames(aggregate_CRS_merge_hg38_gene_sorted)
rownames(aggregate_CRS_merge_hg38_gene_sorted)=aggregate_CRS_merge_hg38_gene_sorted[,1]
rownames(aggregate_CRS_merge_hg38_gene_sorted)
colnames(aggregate_CRS_merge_hg38_gene_sorted)
aggregate_CRS_merge_hg38_gene_sorted = aggregate_CRS_merge_hg38_gene_sorted[,-1]
aggregate_CRS_merge_hg38_gene_sorted = as.matrix(aggregate_CRS_merge_hg38_gene_sorted)
head(aggregate_CRS_merge_hg38_gene_sorted)
dim(aggregate_CRS_merge_hg38_gene_sorted)
Count_CRS_merge_hg38_gene <- count(CRS_merge_hg38_gene, "V1")
aggregate_CRS_merge_hg38_gene_sorted_counted <- cbind(aggregate_CRS_merge_hg38_gene_sorted, Count_CRS_merge_hg38_gene)
head(aggregate_CRS_merge_hg38_gene_sorted_counted)
dim(aggregate_CRS_merge_hg38_gene_sorted_counted)
aggregate_CRS_merge_hg38_gene_sorted_counted <- aggregate_CRS_merge_hg38_gene_sorted_counted[which(aggregate_CRS_merge_hg38_gene_sorted_counted$freq >= 3),]
#write.table(aggregate_CRS_merge_hg38_gene_sorted_counted, "aggregate_CRS_merge_hg38_gene_sorted_counted.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE)
aggregate_CRS_merge_hg38_gene_sorted_counted1 <- aggregate_CRS_merge_hg38_gene_sorted_counted[,1:3]
colnames(aggregate_CRS_merge_hg38_gene_sorted_counted1) <- c("Control_Placenta","KHDC3L_BiCHM","AnCHM")
head(aggregate_CRS_merge_hg38_gene_sorted_counted1)
dim(aggregate_CRS_merge_hg38_gene_sorted_counted1)
write.table(aggregate_CRS_merge_hg38_gene_sorted_counted1, "aggregate_CRS_merge_hg38_gene_sorted_counted1.txt" , quote = FALSE, append = FALSE, sep = "\t")
summary(aggregate_CRS_merge_hg38_gene_sorted_counted1)
summary(aggregate_CRS_merge_hg38_gene_sorted_counted1$Control_Placenta)
summary(aggregate_CRS_merge_hg38_gene_sorted_counted1$KHDC3L_BiCHM)
summary(aggregate_CRS_merge_hg38_gene_sorted_counted1$AnCHM)
wilcox.test(aggregate_CRS_merge_hg38_gene_sorted_counted1$Control_Placenta, aggregate_CRS_merge_hg38_gene_sorted_counted1$KHDC3L_BiCHM)
wilcox.test(aggregate_CRS_merge_hg38_gene_sorted_counted1$Control_Placenta, aggregate_CRS_merge_hg38_gene_sorted_counted1$AnCHM)
aggregate_CRS_merge_hg38_gene_sorted_counted2 <- stack(data.frame(aggregate_CRS_merge_hg38_gene_sorted_counted1))
head(aggregate_CRS_merge_hg38_gene_sorted_counted2)
colnames(aggregate_CRS_merge_hg38_gene_sorted_counted2) <- c("Methylation", "Datasets")
head(aggregate_CRS_merge_hg38_gene_sorted_counted2)
ggboxplot(aggregate_CRS_merge_hg38_gene_sorted_counted2, x = "Datasets", y = "Methylation",
          title = "AvgMethGenes", ylab = "Methylation",
          palette =c("#4292A6", "#AEE2EF", "#E6FBFF"),fill = "Datasets", bxp.errorbar = TRUE, ggtheme = theme_classic()) + ylim(0,1)
#ggsave("Boxplot_CRS_merge_hg38_updated_genes_filt.svg", width=12*1.25, height=15*1.25, units="cm", dpi=96)


