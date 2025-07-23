#!/bin/bash

exec > >(tee log/extract_reads.log) 2>&1

for num in 52 54 60; do
    echo "Processing sample $num"
    minimap2 -t 16 -ax map-ont \
        old_run/barcode${num}_gemy.fasta barcode${num}_cassava_unmapped.fastq > barcode${num}.mapped.sam
    samtools view -b -F 4 barcode${num}.mapped.sam > barcode${num}.mapped.bam
    samtools fastq barcode${num}.mapped.bam > old_run/barcode${num}.gemy.fastq
    rm barcode${num}.mapped.sam barcode${num}.mapped.bam
    echo "Extracted reads for sample $num"
done
echo "All samples processed successfully."

