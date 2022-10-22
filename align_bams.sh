#!/bin/bash

samples=$(ls /array/psivakumar/MSA_RNASeq/merged_fastq/*fastq.gz)
names=$(ls /array/psivakumar/MSA_RNASeq/merged_fastq/*fastq.gz | cut -d"/" -f6 | cut -d"." -f1)
ids=$(ls /array/psivakumar/MSA_RNASeq/merged_fastq/*fastq.gz | cut -d"/" -f6 | cut -d"_" -f2 | sort -u)
arr_s=($samples)
arr_n=($names)
arr_i=($ids)	
length=$(($(echo $ids | wc -w) -1))

echo $length
for i in `seq 0 $length`
do
	mates=$(ls /array/psivakumar/MSA_RNASeq/merged_fastq/*${arr_i[i]}*)
	arr_m=($mates)
	/data/kronos/NGS_Software/STAR-2.7.0c/bin/Linux_x86_64_static/STAR \
	--readFilesIn ${arr_m[0]} ${arr_m[1]} \
	--readFilesCommand zcat \
	--genomeLoad LoadAndKeep \
	--genomeDir /data/kronos/NGS_Reference/STAR_genomes/human_hg38/ \
	--outFileNamePrefix /array/psivakumar/MSA_RNASeq/bams/${arr_i[i]}_ \
	--outSAMtype BAM SortedByCoordinate \
	--outSAMunmapped Within \
	--limitBAMsortRAM 50000000000
	samtools index /array/psivakumar/MSA_RNASeq/bams/${arr_i[i]}_Aligned.sortedByCoord.out.bam
done
