#!/bin/bash

##PBS configure
# PBS -N calculate_FRiP
# PBS -j oe
# PSB -q batch
# PBS -S /bin/sh
# PBS -l nodes=1:ppn=10

## Source module environment and load tools
source /etc/profile.d/modules.sh
module load samtools/1.9
module load bedtools/2.30.0

## Set variables
bam_dir=~/ChinmoST/output/chinmo_cut_tag/bam/
peaks_dir=~/ChinmoST/output/chinmo_cut_tag/MACS2/
samples=(ChinmoST_D57_rep1 ChinmoST_D57_rep2 ChinmoST_D911_rep1 ChinmoST_D911_rep2 WT_D57_rep1 WT_D57_rep2 WT_D911_rep1 WT_D911_rep2)

for sample in ${samples[@]};
do
	## Convert BAM files to BED files
	bedtools bamtobed -i ${bam_dir}${sample}_sorted_rmDup_mapped_rmbl_shift.bam \
	> ${bam_dir}${sample}_sorted_rmDup_mapped_rmbl_shift.bed

	total_reads=$(wc -l ${bam_dir}${sample}_sorted_rmDup_mapped_rmbl_shift.bed |awk '{print $1}')

	## Fraction of reads in MACS2 peaks (FRiP)
	reads_in_MACS2_peaks=$(bedtools intersect -a ${bam_dir}${sample}_sorted_rmDup_mapped_rmbl_shift.bed \
        -b ${peaks_dir}${sample}_peaks.narrowPeak -u |wc -l|awk '{print $1}')

	echo "==> ${sample} MACS2 FRiP value (all peaks):" $(bc <<< "scale=2;100*$reads_in_MACS2_peaks/$total_reads")'%'
done

## Unload tools
module unload samtools/1.9
module unload bedtools/2.30.0