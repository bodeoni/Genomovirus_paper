#!/bin/bash

mkdir -p output/v1/log output/v1/seed output/v1/polished output/v1/final

exec > >(tee output/v1/log/output.log) 2>&1

for num in 52 54 60; do
    echo ">>> Selecting seed for barcode ${num}"

    fastq=fastq/sample${num}.fastq
    
    echo "[*] Activating nanopore environment"
    source /home/bodeoni/miniconda3/bin/activate base

    echo "=================================="
    echo "[*] Filtering long reads >2.2kb"
    echo "=================================="
    NanoFilt --maxlength 2200 -q 15 $fastq > output/v1/seed/sample_${num}.trimmed.fastq

    echo "=================================="
    echo "[*] Extracting longest read as seed"
    echo "=================================="
    
    seqkit sort -l output/v1/seed/sample_${num}.trimmed.fastq | tail -n 4 | seqkit fq2fa > output/v1/seed/sample_${num}.seed.fasta
    seqkit replace -p '^.*$' output/v1/seed/sample_${num}.seed.fasta -r "sample_${num}_seed" > output/v1/seed/sample_${num}.seed.renamed.fasta
    
    echo "=================================="
    echo "[*] Mapping reads to seed (round 1)"
    echo "=================================="
    minimap2 -t 16 -x map-ont output/v1/seed/sample_${num}.seed.renamed.fasta $fastq > output/v1/seed/sample_${num}.seed1.paf
    
    echo "=================================="
    echo "[*] Polishing with racon (round 1)"
    echo "=================================="
    #source /opt/tljh/user/bin/activate racon
    racon -t 16 $fastq output/v1/seed/sample_${num}.seed1.paf output/v1/seed/sample_${num}.seed.renamed.fasta > output/v1/polished/sample_${num}.racon1.fasta

    echo "=================================="
    echo "[*] Mapping reads to polished seed (round 2)"
    echo "=================================="
    #source /opt/tljh/user/bin/activate nanopore
    minimap2 -t 16 -x map-ont output/v1/polished/sample_${num}.racon1.fasta $fastq > output/v1/seed/sample_${num}.seed2.paf

    echo "=================================="
    echo "[*] Polishing with racon (round 2)"
    echo "=================================="
    #source /opt/tljh/user/bin/activate racon
    racon -t 16 $fastq output/v1/seed/sample_${num}.seed2.paf output/v1/polished/sample_${num}.racon1.fasta > output/v1/polished/sample_${num}.racon2.fasta

    echo "=================================="
    echo "[*] Final polishing with Medaka"
    echo "=================================="
    source /home/bodeoni/miniconda3/bin/activate medaka
    medaka_consensus -fx -i $fastq -d output/v1/polished/sample_${num}.racon2.fasta -o output/v1/final/sample_${num} -t 2
    
    #rename file
    cp output/v1/final/sample_${num}/consensus.fasta output/v1/final/sample_${num}_final.fasta

    echo "✅ Finished sample ${num}"
done

echo "All samples complete."