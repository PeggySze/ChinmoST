#!/bin/bash

## set variables
input_dir=~/ChinmoST/input/data/chinmo_cut_tag/
samples=(ChinmoST_D57_rep1 ChinmoST_D57_rep2 ChinmoST_D911_rep1 ChinmoST_D911_rep2 WT_D57_rep1 WT_D57_rep2 WT_D911_rep1 WT_D911_rep2)

for sample in ${samples[@]};
do
	qsub -N ${sample} -l nodes=1:ppn=5 \
	-F "${input_dir}${sample}_R1.fq.gz ${input_dir}${sample}_R2.fq.gz" \
	~/ChinmoST/bin/Cut&tag_processing.sh
done
