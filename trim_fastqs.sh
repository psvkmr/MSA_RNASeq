#!/bin/bash

samples=$(find /mnt/MSA-data-2/hhoulden-jvandrovcova-results/MSA -name "*.fastq.gz")
arr_s=($samples)
length=$(echo $samples | wc -w)

for i in `seq 0 $(expr ${length} - 1)`
do
	trim...
done
