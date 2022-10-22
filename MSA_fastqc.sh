#!/bin/bash

find /mnt/MSA-data-2/hhoulden-jvandrovcova-results/MSA -name "*.fastq.gz" | xargs /data/kronos/Genetics_Software/fastqc_v0.11.8/FastQC/fastqc --outdir=/array/psivakumar/MSA_RNASeq/fastqc

