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


## 
