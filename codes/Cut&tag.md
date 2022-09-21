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


```


## Binding and Expression Target Analysis (BETA)
- script : 
```R

```