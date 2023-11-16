#!/bin/bash

##PBS configure
# PBS -N filter_peaks
# PBS -j oe
# PSB -q batch
# PBS -S /bin/sh
# PBS -l nodes=1:ppn=10

## Set variables
peaks_dir=~/ChinmoST/output/chinmo_cut_tag/MACS2/
samples=(ChinmoST_D57_rep1 ChinmoST_D57_rep2 ChinmoST_D911_rep1 ChinmoST_D911_rep2 WT_D57_rep1 WT_D57_rep2 WT_D911_rep1 WT_D911_rep2)

for sample in ${samples[@]};
do
	# keep peaks in autosomes and chrX, chrY
	awk '$1=="chr2L" || $1=="chr2R" || $1=="chr3L" || $1=="chr3R" || $1=="chr4" || $1=="chrX" || $1=="chrY"' ${peaks_dir}${sample}_peaks.narrowPeak > ${peaks_dir}${sample}_filtered_peaks.narrowPeak
done