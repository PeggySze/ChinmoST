## Run Cut&Tag data processing pipeline
- The Cut&Tag data processing pipeline consists of the following parts :
	- quanlity control and trimming of raw reads by fastp
	- mapping trimmed reads by bowtie2
	- remove duplicates by picard
	- filter and keep the uniquely mapped reads
	- filter out reads in Blacklist Regions
	- shift reads and create bigwig by deeptools
	- call peaks by MACS2
- script : 
	- Cut&tag_processing.sh
	- run_Cut&tag_pipeline.sh


## Filter peaks and calculate Fraction of reads in peaks (FRiP)
- keep peaks in autosomes,chrX,chrY
	- script : filter_peaks.sh
- calculate Fraction of reads in peaks (FRiP) value
	- script : calculate_FRiP.sh
```R
## Load required packages
library(RColorBrewer)


## Figure S6B - number of MACS2 filtered peaks in each sample
peak_files <- list.files("~/ChinmoST/output/chinmo_cut_tag/MACS2",pattern="*filtered_peaks.narrowPeak",full.names=TRUE)
peak_samples <- sapply(1:length(peak_files),function(i){
  unlist(strsplit(unlist(strsplit(peak_files[i],"/"))[6],"_filtered"))[1]
  })
tissues <-  sapply(1:length(peak_files),function(i){
  unlist(strsplit(unlist(strsplit(peak_files[i],"/"))[6],"_"))[1]
  })
replicates <- sapply(1:length(peak_files),function(i){
  unlist(strsplit(unlist(strsplit(peak_files[i],"/"))[6],"_"))[3]
  })
days <- sapply(1:length(peak_files),function(i){
  unlist(strsplit(unlist(strsplit(peak_files[i],"/"))[6],"_"))[2]
  })
peakN_df <- c()
for (i in 1:length(peak_files)){
  peak_df <- read.table(peak_files[i])
  peakN_df <- data.frame(sample=peak_samples[i],tissue=tissues[i],replicate=replicates[i],day=days[i],peakN=nrow(peak_df)) %>% rbind(peakN_df,.)
}
peakN_df$tissue <- factor(peakN_df$tissue,levels=c("WT","ChinmoST"))
peakN_df$day <- factor(peakN_df$day,levels=c("D57","D911"),labels=c("D5-7","D9-11"))
newpalette <- brewer.pal(8,"Set2")[c(1,3)]
pdf("~/ChinmoST/output/chinmo_cut_tag/MACS2/D57_D911_filtered_peaks_number.pdf")
ggplot(data=peakN_df,aes(x=day,y=peakN,fill=tissue,group=replicate)) +
  geom_bar(stat="identity",width=0.4,position=position_dodge(0.6))+
  facet_grid(tissue ~ .,switch = "y")+
  theme_bw() +
  labs(x="Day",y="Number of peaks",fill="Condition",title="Number of peaks") + 
  scale_fill_manual(values=newpalette)+
  theme(plot.title = element_text(hjust = 0.5,size=14),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    panel.background=element_rect(fill='transparent', color='black',linetype="solid"),strip.background=element_rect(colour="black",fill="#f0f0f0",linetype="solid"),strip.text = element_text(face="bold",size=12),legend.position = "bottom",axis.title = element_text(size=12),axis.text = element_text(color="black",size=10))
dev.off()


## Figure S6C - Fraction of reads in MACS2 called peaks
df <- read.table("~/ChinmoST/bin/calculate_FRiP.o*",header = FALSE, fill = TRUE,sep=" ",stringsAsFactors=FALSE)
df <- df[,c(2,8)]
df$V8 <- gsub("%$","",df$V8)
df$Day <- sapply(1:nrow(df),function(i){
	unlist(strsplit(df$V2[i],"_"))[2]
	})
df$tissue <- sapply(1:nrow(df),function(i){
	unlist(strsplit(df$V2[i],"_"))[1]
	})
df$replicates <- sapply(1:nrow(df),function(i){
	unlist(strsplit(df$V2[i],"_"))[3]
	})
df$Day <- factor(df$Day,levels=c("D57","D911"),labels=c("D5-7","D9-11"))
df$tissue <- factor(df$tissue,levels=c("WT","ChinmoST"))
df$replicates <- factor(df$replicates,levels=c("rep1","rep2"))
df$V8 <- as.numeric(df$V8)
newpalette <- brewer.pal(8,"Set2")[c(1,3)]
pdf("~/ChinmoST/output/chinmo_cut_tag/MACS2/D57_D911_MACS2_all_peaks_FRiP.pdf")
ggplot(data=df,aes(x=Day,y=V8,fill=tissue,group=replicates)) +
  geom_bar(stat="identity",width=0.4,position=position_dodge(0.6))+
  facet_grid(tissue ~ .,switch = "y")+
  theme_bw() +
  labs(x="Day",y="FRiP(%)",title="Fraction of reads in peaks",fill="Condition") + 
  scale_fill_manual(values=newpalette)+
  theme(plot.title = element_text(hjust = 0.5,size=14),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    panel.background=element_rect(fill='transparent', color='black',linetype="solid"),strip.background=element_rect(colour="black",fill="#f0f0f0",linetype="solid"),strip.text = element_text(face="bold",size=12),legend.position = "bottom",axis.title = element_text(size=12),axis.text = element_text(color="black",size=10))
dev.off()
```


## Assess reproducibility between biological replicates
- Assess reproducibility by Irreproducibility Discovery Rate (IDR) framework and retain reproducible peaks passing IDR cutoff of 0.05
	- script : filtered_peaks_IDR.sh
- Merge reproducible peaks in different time points for the same condition (WT & ChinmoST testes)
```shell
module load bedtools/2.30.0
cd ~/ChinmoST/output/chinmo_cut_tag/IDR
cat Reproducible_WT_D57.idr Reproducible_WT_D911.idr | sort -k1,1 -k2,2n | bedtools merge -i - > merged_D57_D911_WT_reproducible_peaks.bed
cat Reproducible_ChinmoST_D57.idr Reproducible_ChinmoST_D911.idr | sort -k1,1 -k2,2n | bedtools merge -i - > merged_D57_D911_ChinmoST_reproducible_peaks.bed
```
```R
## Load required packages
library(RColorBrewer)
library(GenomicRanges)
library(ChIPseeker)
library(org.Dm.eg.db)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene

## Figure 5B - WT reproducible peaks annotation
files <- list.files("~/ChinmoST/output/chinmo_cut_tag/IDR",pattern="^Reproducible_WT_D*.idr$",full.names=TRUE)
idr_df <- c()
for (i in 1:length(files)){
  idr_df <- read.table(files[i],header=FALSE) %>% rbind(idr_df,.)
}
peaks.gr <- GRanges(idr_df$V1, IRanges(start = idr_df$V2+1, end = idr_df$V3), strand = "*")
peakAnno <- annotatePeak(peaks.gr,tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Dm.eg.db")
p <- plotAnnoBar(peakAnno)
plot_df <- p$data
plot_df$Feature <- factor(plot_df$Feature,levels=c("Promoter (<=1kb)","Promoter (1-2kb)","Promoter (2-3kb)","5' UTR","3' UTR","1st Exon","Other Exon","1st Intron","Other Intron","Downstream (<=300)","Distal Intergenic"))
newpalette <- brewer.pal(12,"Set3")[1:11]
pdf("~/ChinmoST/output/chinmo_cut_tag/IDR/WT_D57_D911_peaks_annotation_piechart.pdf")
ggplot(data=plot_df,aes(x='Feature',y=Frequency,fill=Feature)) +
  geom_bar(stat = 'identity', position = 'stack',color="black",size=0.02)  + 
  coord_polar(theta = 'y') + 
  labs(x = '', y = '', title = '') +
  scale_fill_manual(values = newpalette)+
  theme_bw() + 
  theme(axis.text = element_blank(),axis.ticks = element_blank(),panel.grid=element_blank())
dev.off()

## Figure S6E - number of merged reproducible peaks in WT and ChinmoST testes, respectively
merged_D57_D911_WT_peaks <- read.table("~/ChinmoST/output/chinmo_cut_tag/IDR/merged_D57_D911_WT_reproducible_peaks.bed",header=FALSE)
merged_D57_D911_Mu_peaks <- read.table("~/ChinmoST/output/chinmo_cut_tag/IDR/merged_D57_D911_ChinmoST_reproducible_peaks.bed",header=FALSE)
df <- data.frame(tissue=c("WT","Mutant"),PeakN=c(nrow(merged_D57_D911_WT_peaks),nrow(merged_D57_D911_Mu_peaks)))
df$tissue <- factor(df$tissue,levels=c("WT","Mutant"))
pdf("~/ChinmoST/output/chinmo_cut_tag/IDR/WT_Mu_D57_D911_merged_reproducible_peaks_number.pdf",height=4,width=6)
ggplot(data=df,aes(x=tissue,y=PeakN,fill=tissue))+
  geom_bar(stat="identity",width=0.3)+
  theme_bw() +
  labs(x="Condition",y="Number of peaks",title="Number of merged reproducible peaks",fill="Condition") + 
  scale_fill_manual(values=newpalette)+
  theme(plot.title = element_text(hjust = 0.5,size=14),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    panel.background=element_rect(fill='transparent', color='black',linetype="solid"),strip.background=element_rect(colour="black",fill="#f0f0f0",linetype="solid"),strip.text = element_text(face="bold",size=12),legend.position = "bottom",axis.title = element_text(size=12),axis.text = element_text(color="black",size=10))
dev.off()
```


## Chinmo binding sites motif analysis
- script : 
```R


```


## Footprinting analysis by HINT-ATAC
- script : 
```R


```


## Differential peak analysis by DiffBind
```R
## Load required packages
library(DiffBind)
library(RColorBrewer)
library(circlize)
library(ggplot2)
library(dplyr)
library(stringr)
library(xlsx)
library(ComplexHeatmap)

## Reading in the peaksets
bam_files <- list.files("~/ChinmoST/output/chinmo_cut_tag/bam",pattern="*_sorted_rmDup_mapped_rmbl_shift.bam$",full.names=TRUE)
bam_files <- bam_files[c(grep("ChinmoST",bam_files),grep("WT",bam_files))]
peak_files <- list.files("~/ChinmoST/output/chinmo_cut_tag/MACS2",pattern="*_filtered_peaks.narrowPeak",full.names=TRUE)
peak_files <- peak_files[c(grep("ChinmoST",peak_files),grep("WT",peak_files))]
bam_samples <- sapply(1:length(bam_files),function(i){
  sample <- unlist(strsplit(unlist(strsplit(bam_files[i],"/"))[6],"_sorted"))[1]
  sample
  })
peak_samples <- sapply(1:length(peak_files),function(i){
  sample <- unlist(strsplit(unlist(strsplit(peak_files[i],"/"))[6],"_filtered"))[1]
  sample
  })
all(bam_samples==peak_samples)
samples <- bam_samples
tissues <- sapply(1:length(samples),function(i){
	unlist(strsplit(samples[i],"_"))[1]
	})
replicates <- c(1:4,1:4)
df <- data.frame(SampleID=samples,Tissue=tissues,Factor="chinmo",Replicate=replicates,bamReads=bam_files,Peaks=peak_files,PeakCaller="narrow")
dbaObj <- dba(sampleSheet=df)

## calculate a binding matrix with scores based on read counts for every sample
# bUseSummarizeOverlaps: to use a more standard counting procedure than the built-in one by default
dbaObj <- dba.count(dbaObj, bUseSummarizeOverlaps=TRUE)

## Figure S6A - Correlation between WT and Mutant samples
pdf("~/ChinmoST/output/chinmo_cut_tag/diffbind/WT_Mu_count_scores_correlation_heatmap.pdf",height=8,width=8)
p <- dba.plotHeatmap(dbaObj,cexRow=0.8,cexCol=0.8)
dev.off()

days <- sapply(1:length(samples),function(i){
  unlist(strsplit(samples[i],"_"))[2]
  })
sample_types <- sapply(1:length(samples),function(i){
  unlist(strsplit(samples[i],"_"))[1]
  })
ha <- HeatmapAnnotation(
  df = data.frame(
      Day=days,
      Tissue=sample_types,
      row.names=samples),
  col = list(
    Day=c("D57"=brewer.pal(9,"Reds")[5],"D911"=brewer.pal(9,"Reds")[7]),
    Tissue=c("WT"=brewer.pal(8,"Set2")[1],"ChinmoST"=brewer.pal(8,"Set2")[3])),
  annotation_legend_param = list(Day = list(nrow = 1,title_position = "topcenter"),Tissue=list(nrow = 1,title_position = "topcenter"))
  )
pdf("~/ChinmoST/output/chinmo_cut_tag/diffbind/WT_Mu_count_scores_correlation_heatmap.pdf",height=8,width=7)
h <- Heatmap(p,top_annotation=ha,name="Pearson correlation",row_names_gp = gpar(fontsize = 9),column_names_gp = gpar(fontsize = 9),col=colorRamp2(breaks=c(0,0.2,0.4,0.6,0.8,1),colors=c("#FFFFFF",brewer.pal(9,"Greens")[c(2,4,6,8,9)])),heatmap_legend_param = list(direction = "horizontal",legend_width = unit(4, "cm"), title_position = "topcenter"))
draw(h,heatmap_legend_side = "bottom",annotation_legend_side = "bottom")
dev.off()

## Establishing a model design and contrast
dbaObj <- dba.contrast(dbaObj,categories=DBA_TISSUE,minMembers = 2)

## Performing the differential enrichment analysis
dbaObj <-  dba.analyze(dbaObj, method=DBA_ALL_METHODS)
dba.show(dbaObj, bContrasts=TRUE)

## Retrieving the differentially bound sites (DBSs), adjusted p-values (FDR) <0.05,abs(FC)>=log2(1.5)
Mu_vs_WT_DBSs <- dba.report(dbaObj, method=DBA_EDGER, contrast = 1,fold=log2(1.5))

# add annotation
library(ChIPseeker)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(org.Dm.eg.db)
txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
peaks_df <- as.data.frame(dba.peakset(dbaObj, bRetrieve=TRUE))
peaks_df$seqnames <- factor(peaks_df$seqnames,levels=c("chr2L","chr2R","chr3L","chr3R","chr4","chrX","chrY"))
peaks_df <- peaks_df %>% arrange(seqnames,start)
consensus.gr <- GRanges(
  seqnames=peaks_df$seqnames, 
  ranges=IRanges(peaks_df$start, peaks_df$end)
  ) 
anno <- annotatePeak(consensus.gr,tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Dm.eg.db")
anno <- as.data.frame(anno@anno)

Mu_vs_WT_DBSs_df <- as.data.frame(dba.report(dbaObj, method=DBA_EDGER, contrast = 1,fold=log2(1.5)))
Mu_vs_WT_DBSs_df$DBS_type <- ""
Mu_vs_WT_DBSs_df$DBS_type[which(Mu_vs_WT_DBSs_df$Fold>=log2(1.5)&Mu_vs_WT_DBSs_df$FDR<0.05)] <- "Mu up DBS"
Mu_vs_WT_DBSs_df$DBS_type[which(Mu_vs_WT_DBSs_df$Fold<=-log2(1.5)&Mu_vs_WT_DBSs_df$FDR<0.05)] <- "Mu down DBS"
Mu_vs_WT_DBSs_df$DBS_type <- ifelse(Mu_vs_WT_DBSs_df$DBS_type=="","Nonsignificant",Mu_vs_WT_DBSs_df$DBS_type)
Mu_vs_WT_DBSs_df <- left_join(Mu_vs_WT_DBSs_df,anno,by=c("seqnames","start","end")) %>% arrange(desc(Fold),FDR)
write.xlsx(Mu_vs_WT_DBSs_df,"~/ChinmoST/output/chinmo_cut_tag/diffbind/D57_D911_Mu_vs_WT_DBSs.xlsx",row.names=FALSE)

## Retrieving edgeR full results
# Conc: Concentration - mean read concentration over all the samples (the default calculation uses log2 normalized ChIP read counts with control read counts subtracted)
Mu_vs_WT_res <- as.data.frame(dba.report(dbaObj, method=DBA_EDGER, contrast = 1,th=1))
Mu_vs_WT_res$DBS_type <- ""
Mu_vs_WT_res$DBS_type[which(Mu_vs_WT_res$Fold>=log2(1.5)&Mu_vs_WT_res$FDR<0.05)] <- "Mu up DBS"
Mu_vs_WT_res$DBS_type[which(Mu_vs_WT_res$Fold<=-log2(1.5)&Mu_vs_WT_res$FDR<0.05)] <- "Mu down DBS"
Mu_vs_WT_res$DBS_type <- ifelse(Mu_vs_WT_res$DBS_type=="","Nonsignificant",Mu_vs_WT_res$DBS_type)
Mu_vs_WT_res <- left_join(Mu_vs_WT_res,anno,by=c("seqnames","start","end"))

# downregulated DBSs in mutant testes
Mu_down_DBS_bed <- Mu_vs_WT_res %>%
  dplyr::filter(DBS_type=="Mu down DBS") %>%
  dplyr::mutate(bed_start=start-1,bed_end=end) %>%
  dplyr::select(seqnames,bed_start,bed_end)
write.table(Mu_down_DBS_bed,"~/ChinmoST/output/chinmo_cut_tag/diffbind/D57_D911_Mu_down_DBSs.bed",row.names=FALSE,sep="\t",quote=FALSE,col.name=FALSE)

# upregulated DBSs in mutant testes
Mu_up_DBS_bed <- Mu_vs_WT_res %>%
  dplyr::filter(DBS_type=="Mu up DBS") %>%
  dplyr::mutate(bed_start=start-1,bed_end=end) %>%
  dplyr::select(seqnames,bed_start,bed_end)
write.table(Mu_up_DBS_bedd,"~/ChinmoST/output/chinmo_cut_tag/diffbind/D57_D911_Mu_up_DBSs.bed",row.names=FALSE,sep="\t",quote=FALSE,col.name=FALSE)

## get Mutant vs WT DBSs summits 
summit_Mu_vs_WT_dbaObj <- dba.count(dbaObj, summits=TRUE, score=DBA_SCORE_SUMMIT_POS)
Mu_vs_WT_summits <- as.data.frame(dba.peakset(summit_Mu_vs_WT_dbaObj, bRetrieve=TRUE))
Mu_vs_WT_DBSs_summits <- left_join(Mu_vs_WT_DBSs_df[,c(1:3,9,12)],Mu_vs_WT_summits,by=c("seqnames","start","end"))
Mu_vs_WT_DBSs_summits$summit <- round(apply(Mu_vs_WT_DBSs_summits[,8:15],1,mean))
Mu_up_DBSs_summits <- Mu_vs_WT_DBSs_summits %>%
    dplyr::filter(DBS_type=="Mu up DBS") %>%
    dplyr::arrange(desc(Fold))
Mu_up_DBSs_summits_bed <- Mu_vs_WT_DBSs_summits %>%
    dplyr::filter(DBS_type=="Mu up DBS") %>%
    dplyr::arrange(desc(Fold)) %>%
    dplyr::mutate(summit_start=summit-1) %>%
    dplyr::select(seqnames,summit_start,summit)
Mu_down_DBSs_summits <- Mu_vs_WT_DBSs_summits %>%
    dplyr::filter(DBS_type=="Mu down DBS") %>%
    dplyr::arrange(Fold)
Mu_down_DBSs_summits_bed <- Mu_vs_WT_DBSs_summits %>%
    dplyr::filter(DBS_type=="Mu down DBS") %>%
    dplyr::arrange(Fold) %>%
    dplyr::mutate(summit_start=summit-1) %>%
    dplyr::select(seqnames,summit_start,summit)

## Figure 5E,F - Enriched Heatmaps over Mutant vs WT DBSs 

```


## Binding and Expression Target Analysis (BETA)
- script : 
```R

```