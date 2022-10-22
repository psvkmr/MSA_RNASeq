#!/bin/bash

samples=$(ls -d /mnt/MSA-data-2/hhoulden-jvandrovcova-results/MSA/*/*)
names=$(ls -d /mnt/MSA-data-2/hhoulden-jvandrovcova-results/MSA/*/* | cut -d"/" -f7)
arr_s=($samples)
arr_n=($names)
length=$(echo $samples | wc -w)

for i in `seq 0 $(expr ${length} - 1)`
do
	reads=$(ls ${arr_s["$i"]} | cut -d"_" -f4 | sort -u)
	arr_r=($reads)
	length_r=$(echo $reads | wc -w)
	for j in `seq 0 $(expr ${length_r} - 1)`
	do 
		all_reads=$(ls ${arr_s["$i"]}/*${arr_r["$j"]}*) 
		cat $all_reads > /array/psivakumar/MSA_RNASeq/merged_fastq/${arr_n["$i"]}_${arr_r["$j"]}.fastq.gz
	done
done
