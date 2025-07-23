#!/bin/bash

mkdir -p output/v2/log output/v2/seed output/v2/polished output/v2/final

exec > >(tee output/v2/log/seed_correct.log) 2>&1

for num in 52 54 60; do
    echo ">>> Selecting seed for barcode ${num}"

    fastq=fastq/sample${num}.fastq
    
    echo "[*] Activating nanopore environment"
    source /home/bodeoni/miniconda3/bin/activate base

    echo "=================================="
    echo "[*] Filtering long reads >2.4kb"
    echo "=================================="
    NanoFilt --maxlength 2200 -q 15 $fastq > output/v2/seed/sample_${num}.trimmed.fastq

    echo "=================================="
    echo "[*] Extracting longest read as seed"
    echo "=================================="
    seqkit sort -l output/v2/seed/sample_${num}.trimmed.fastq | tail -n 4 | seqkit fq2fa > output/v2/seed/sample_${num}.seed.fasta
    seqkit replace -p '^.*$' output/v2/seed/sample_${num}.seed.fasta -r "sample_${num}_seed" > output/v2/seed/sample_${num}.seed.renamed.fasta

    echo "=================================="
    echo "[*] Mapping reads to seed (round 1)"
    echo "=================================="
    minimap2 -t 16 -x map-ont output/v2/seed/sample_${num}.seed.renamed.fasta output/v2/seed/sample_${num}.trimmed.fastq > output/v2/seed/sample_${num}.seed1.paf

    echo "=================================="
    echo "[*] Polishing with racon (round 1)"
    echo "=================================="
    #source /opt/tljh/user/bin/activate racon
    racon -t 16 output/v2/seed/sample_${num}.trimmed.fastq output/v2/seed/sample_${num}.seed1.paf output/v2/seed/sample_${num}.seed.renamed.fasta > output/v2/polished/sample_${num}.racon1.fasta

    echo "=================================="
    echo "[*] Mapping reads to polished seed (round 2)"
    echo "=================================="
    #source /opt/tljh/user/bin/activate nanopore
    minimap2 -t 16 -x map-ont output/v2/polished/sample_${num}.racon1.fasta output/v2/seed/sample_${num}.trimmed.fastq > output/v2/seed/sample_${num}.seed2.paf

    echo "=================================="
    echo "[*] Polishing with racon (round 2)"
    echo "=================================="
    #source /opt/tljh/user/bin/activate racon
    racon -t 16 output/v2/seed/sample_${num}.trimmed.fastq output/v2/seed/sample_${num}.seed2.paf output/v2/polished/sample_${num}.racon1.fasta > output/v2/polished/sample_${num}.racon2.fasta

    echo "=================================="
    echo "[*] Final polishing with Medaka"
    echo "=================================="
    source /home/bodeoni/miniconda3/bin/activate medaka
    medaka_consensus -i output/v2/seed/sample_${num}.trimmed.fastq -d output/v2/polished/sample_${num}.racon2.fasta -o output/v2/final/sample_${num} -t 4
    
    #rename file
    cp output/v2/final/sample_${num}/consensus.fasta output/v2/final/sample_${num}_final.fasta

    echo "✅ Finished sample ${num}"
done

echo "All samples complete."
