#!/bin/bash

##PBS configure
# PBS -N filtered_peaks_IDR
# PBS -j oe
# PSB -q batch
# PBS -S /bin/sh
# PBS -l nodes=1:ppn=10

echo "`date +%Y/%m/%d_%H:%M:%S`"  ## Record the date and time
uname -sa  ## Information about the operating system
set -ex  ## Log everything,and quit if there is any error
START=$(date +%s.%N)

## Source module environment and load tools
source /etc/profile.d/modules.sh
module load idr/2.0.3

## Set variables
peaks_dir=~/ChinmoST/output/chinmo_cut_tag/MACS2/
output_dir=~/ChinmoST/output/chinmo_cut_tag/IDR/
samples=(ChinmoST_D57_rep1 ChinmoST_D57_rep2 ChinmoST_D911_rep1 ChinmoST_D911_rep2 WT_D57_rep1 WT_D57_rep2 WT_D911_rep1 WT_D911_rep2)
tissues=(ChinmoST_D57 ChinmoST_D911 WT_D57 WT_D911)

## Sort narrowPeak files by -log10(p-value)
for sample in ${samples[@]};
do
	sort -k8,8nr ${peaks_dir}${sample}_filtered_peaks.narrowPeak \
	> ${output_dir}${sample}_sorted_filtered_peaks.narrowPeak
done

## Access reproducibility of replicates by IDR
for tissue in ${tissues[@]};
do
	idr --samples ${output_dir}${tissue}_rep1_sorted_filtered_peaks.narrowPeak ${output_dir}${tissue}_rep2_sorted_filtered_peaks.narrowPeak \
	--input-file-type narrowPeak \
	--rank p.value \
	--output-file ${output_dir}${tissue}.idr \
	--plot \
	--log-output-file ${output_dir}${tissue}_idr.log
done


## Retain reproducible peaks passing IDR cutoff of 0.05
for tissue in ${tissues[@]};
do
	awk '{if($5 >= 540) print $0}' ${output_dir}${tissue}.idr > ${output_dir}Reproducible_${tissue}.idr
done


## Unload tools
module unload idr/2.0.3

END=$(date +%s.%N)
Duration=$(echo "$END - $START" | bc)
echo "`date +%Y/%m/%d_%H:%M:%S` Run completed"
