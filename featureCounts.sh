#!/bin/bash

#$1=bam_folder "/" appended
#$2=output matrix prefix

bam_folder=$1
bams=$(ls $1*.bam)

/data/kronos/NGS_Software/subread-1.6.4-Linux-x86_64/bin/featureCounts \
-a /data/kronos/NGS_Reference/STAR_genomes/Homo_sapiens.GRCh38.95.gtf \
-F GTF \
-g gene_id \
-p \
-s 1 \
--ignoreDup \
-o $2 \
$bams
