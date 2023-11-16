#!/bin/bash
# This script takes two fastq files (paired-end) of CUT&Tag data,runs fastp (QC and trimming),bowtie2 (mapping trimmed reads),deeptools (shifting reads and creating bigwig) and MACS2 (calling peaks).
# USAGE: 
# sh Cut_Tag.sh <name of fastq file 1> <name of fastq file 2>
# qsub -N <samplename> -l nodes=1:ppn=3 -F "<name of fastq file 1> <name of fastq file 2>" Cut_Tag.sh
# take consideration of names of fastq files to modify ${samplename}
# take consideration of analyzed species to modify genome and blacklist file
# change ${output_dir} when analyzing different projects

##PBS configure
#PSB -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=3

echo "`date +%Y/%m/%d_%H:%M:%S`"  ## Record the date and time
set -x  ## Log everything
START=$(date +%s.%N)

# initialize a variable with an intuitive name to store the name of the input fastq file
fq_1=$1
fq_2=$2

# grab base of filename for naming outputs
samplename=`basename ${fq_1} _R1.fq.gz`
echo "Sample name is $samplename"  

# specify the number of cores to use
cores=8

# directory with the bowtie2 genome index
genome=~/DB/dm6/bowtie2_index/dm6
blacklist=~/DB/dm6/blacklist/dm6-blacklist.v2.bed

# make all of the output directories
# The -p option means no error if existing, make parent directories as needed
output_dir=~/ChinmoST/output/chinmo_cut_tag/
mkdir -p ${output_dir}fastp
mkdir -p ${output_dir}bam
mkdir -p ${output_dir}bigwig
mkdir -p ${output_dir}MACS2

# set up output directories
fastp_out=${output_dir}fastp/
bam_out=${output_dir}bam/
peak_out=${output_dir}MACS2/
bigwig_out=${output_dir}bigwig/

## Source module environment and load tools
source /etc/profile.d/modules.sh
module load bowtie2/2.4.3
module load samtools/1.9
module load picard/2.21.1
module load bedtools/2.30.0

## Quality control and read trimming by fastp
echo "Starting QC and trimming for $samplename"
~/software/fastp -i ${fq_1} \
	-I ${fq_2} \
	-o ${fastp_out}trimmed_${samplename}_R1.fastq.gz \
	-O ${fastp_out}trimmed_${samplename}_R2.fastq.gz \
	-h ${fastp_out}${samplename}_fastp.html \
	-j ${fastp_out}${samplename}_fastp.json


## Map reads to reference genome by bowtie2
echo "Starting mapping for $samplename"
bowtie2 --end-to-end --very-sensitive \
	--no-mixed --no-discordant \
	--phred33 -I 10 -X 700  -p 8 \
	-x ${genome} \
	-1 ${fastp_out}trimmed_${samplename}_R1.fastq.gz \
	-2 ${fastp_out}trimmed_${samplename}_R2.fastq.gz \
	-S ${bam_out}${samplename}.sam &> ${bam_out}${samplename}_bowtie2_summary.txt

## Compress sam files to bam files and sort bam file by samtools
samtools view -S -b ${bam_out}${samplename}.sam \
	| samtools sort -@ 4 -O bam -o ${bam_out}${samplename}_sorted.bam
rm ${bam_out}${samplename}.sam


## Remove duplicates by using picard
echo "Remove duplicates of $samplename"
java -jar picard.jar MarkDuplicates \
	I=${bam_out}${samplename}_sorted.bam \
	O=${bam_out}${samplename}_sorted_rmDup.bam \
	M=${bam_out}${samplename}_rmDup_metrics.txt \
	REMOVE_DUPLICATES=true


## Sorting BAM files by genomic coordinates
mv ${bam_out}${samplename}_sorted_rmDup.bam ${bam_out}${samplename}_rmDup.bam 
samtools sort -@ 4 -O bam -o ${bam_out}${samplename}_sorted_rmDup.bam ${bam_out}${samplename}_rmDup.bam
rm ${bam_out}${samplename}_rmDup.bam


## Filter and keep the uniquely mapped reads
# filtered out multimappers by specifying `[XS] == null`
# filtered out unmapped reads by specifying in the filter `not unmapped`
echo "Keep the uniquely mapped reads of $samplename"
sambamba view -h -t 2 -f bam -F "[XS] == null and not unmapped" \
	${bam_out}${samplename}_sorted_rmDup.bam > ${bam_out}${samplename}_sorted_rmDup_mapped.bam


## index BAM files
cd ${bam_out}
samtools index ${bam_out}${samplename}_sorted_rmDup_mapped.bam

## Filter out reads in Blacklist Regions
echo "Filter out reads in Blacklist Regions for $samplename"
bedtools intersect -a ${bam_out}${samplename}_sorted_rmDup_mapped.bam \
	-b ${blacklist} -v > ${bam_out}${samplename}_sorted_rmDup_mapped_rmbl.bam

## index BAM files
samtools index ${bam_out}${samplename}_sorted_rmDup_mapped_rmbl.bam

## shift BAM files by deeptools alignmentSieve
echo "Shifting reads of $samplename"
alignmentSieve -b ${bam_out}${samplename}_sorted_rmDup_mapped_rmbl.bam \
	-p ${cores} --ATACshift \
	-o ${bam_out}${samplename}_sorted_rmDup_mapped_rmbl_shift.bam

## generate bigwig files
echo "Generating bigwig file for $samplename"
mv ${bam_out}${samplename}_sorted_rmDup_mapped_rmbl_shift.bam ${bam_out}${samplename}_rmDup_mapped_rmbl_shift.bam
samtools sort -@ 4 -O bam -o ${bam_out}${samplename}_sorted_rmDup_mapped_rmbl_shift.bam \
	${bam_out}${samplename}_rmDup_mapped_rmbl_shift.bam
rm ${bam_out}${samplename}_rmDup_mapped_rmbl_shift.bam

samtools index ${bam_out}${samplename}_sorted_rmDup_mapped_rmbl_shift.bam

bamCoverage --bam ${bam_out}${samplename}_sorted_rmDup_mapped_rmbl_shift.bam \
	--numberOfProcessors ${cores} --binSize 10 \
	--normalizeUsing RPKM \
	-o ${bigwig_out}${samplename}_RPKM_normalized.bw


## Peak calling by MACS2
echo "Calling peaks for $samplename"
/public/home/shipy3/miniconda3/conda_software/bin/macs2 callpeak \
	-t ${bam_out}${samplename}_sorted_rmDup_mapped_rmbl_shift.bam \
	-f BAMPE -g dm --keep-dup all \
	-n ${samplename} \
	--outdir ${peak_out} \
	&>${peak_out}${samplename}_MACS2Peaks_summary.txt

## Unload tools
module unload bowtie2/2.4.3
module unload samtools/1.9
module unload picard/2.21.1
module unload bedtools/2.30.0

END=$(date +%s.%N)
Duration=$(echo "$END - $START" | bc)
echo "`date +%Y/%m/%d_%H:%M:%S` Run completed"
echo $Duration
