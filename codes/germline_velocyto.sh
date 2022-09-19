#!/bin/bash

## PBS configure
# PBS -N germline_velocyto
# PBS -j oe
# PSB -q batch
# PBS -S /bin/sh
# PBS -l nodes=1:ppn=20

echo "`date +%Y/%m/%d_%H:%M:%S`"  ## Record the date and time
uname -sa  ## Information about the operating system
START=$(date +%s.%N)

## Activate conda python3 environment

## Source module envrionment and load tools
source /etc/profile.d/modules.sh
module load samtools/1.9

## Set variables
cellranger_dir=~/ChinmoST/output/scRNA/cellranger/
output_dir=~/ChinmoST/output/scRNA/Velocyto/
samples=(WT_D35 WT_D68 WT_D911 ChinmoST_D35 ChinmoST_D68 ChinmoST_D911)


for sample in ${samples[@]};
do
	if [ ! -d ${output_dir}${sample} ]; then
		mkdir ${output_dir}${sample}
	fi

	velocyto run \
	-b ${cellranger_dir}${sample}/outs/filtered_feature_bc_matrix/barcodes.tsv \
	-m ~/DB/dm6/annotation/dm6_rmsk.gtf \
	-o ${output_dir}${sample} \
	${cellranger_dir}${sample}/outs/possorted_genome_bam.bam \
	~/DB/dm6/annotation/dmel-all-r6.33.gtf
done

## Unload tools
module unload samtools/1.9

END=$(date +%s.%N)
Duration=$(echo "$END - $START" | bc)
echo "`date +%Y/%m/%d_%H:%M:%S` Run completed"
echo $Duration

