#!/bin/bash

find /array/psivakumar/MSA_RNASeq/merged_fastq -name "*.fastq.gz" | xargs /data/kronos/Genetics_Software/fastqc_v0.11.8/FastQC/fastqc --outdir=/array/psivakumar/MSA_RNASeq/fastqc/merged_fastqc
