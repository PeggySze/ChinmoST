#!/bin/bash

## PBS configure
# PBS -N testis_scRNA_cellranger
# PBS -j oe
# PSB -q batch
# PBS -S /bin/sh
# PBS -l nodes=1:ppn=20

echo "`date +%Y/%m/%d_%H:%M:%S`"  ## Record the date and time
uname -sa  ## Information about the operating system

## Source module envrionment and load tools
module load cellranger/4.0.0

## Set variables
input_dir=~/ChinmoST/input/data/scRNA/
output_dir=~/ChinmoST/output/scRNA/cellranger/
reference_dir=~/DB/10X_reference/
samples_01=(WT_D35 ChinmoST_D35)
samples_02=(WT_D68 ChinmoST_D68 WT_D911 ChinmoST_D911)

## generate single cell feature counts by cellranger count
cd ${output_dir}
for sample in ${samples_01[@]};
do
	cellranger count --id=${sample} \
	--fastqs=${input_dir}${sample}/ \
	--sample=${sample} \
	--transcriptome=${reference_dir}dm6/ \
	--localcores=8 \
	--localmem=64
done

# Force pipeline to use 5000 cells, bypassing cell detection
for sample in ${samples_02[@]};
do
	cellranger count --id=${sample} \
	--fastqs=${input_dir}${sample}/ \
	--sample=${sample} \
	--transcriptome=${reference_dir}dm6/ \
	--localcores=8 \
	--localmem=64 \
	--force-cells=5000
done

## Unload tools
module unload cellranger/4.0.0

END=$(date +%s.%N)
Duration=$(echo "$END - $START" | bc)
echo "`date +%Y/%m/%d_%H:%M:%S` Run completed"
echo $Duration

