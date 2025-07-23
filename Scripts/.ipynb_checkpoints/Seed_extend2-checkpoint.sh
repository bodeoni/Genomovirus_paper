#!/bin/bash

mkdir -p log seed polished final

exec > >(tee log/seed_correct.log) 2>&1

for num in 52 54 60; do
    echo ">>> Selecting seed for barcode ${num}"

    fastq=fastq/sample${num}.fastq
    
    echo "[*] Activating nanopore environment"
    source /opt/tljh/user/bin/activate nanopore

    echo "[*] Filtering long reads >2.4kb"
    NanoFilt --maxlength 2200 -q 15 $fastq > seed/sample_${num}.trimmed.fastq

    echo "[*] Extracting longest read as seed"
    seqkit sort -l seed/sample_${num}.trimmed.fastq | tail -n 4 | seqkit fq2fa > seed/sample_${num}.seed.fasta
    seqkit replace -p '^.*$' seed/sample_${num}.seed.fasta -r "sample_${num}_seed" > seed/sample_${num}.seed.renamed.fasta

    echo "[*] Mapping reads to seed (round 1)"
    minimap2 -t 16 -x map-ont seed/sample_${num}.seed.renamed.fasta seed/sample_${num}.trimmed.fastq > seed/sample_${num}.seed1.paf

    echo "[*] Polishing with racon (round 1)"
    source /opt/tljh/user/bin/activate racon
    racon -t 16 seed/sample_${num}.trimmed.fastq seed/sample_${num}.seed1.paf seed/sample_${num}.seed.renamed.fasta > polished/sample_${num}.racon1.fasta

    echo "[*] Mapping reads to polished seed (round 2)"
    source /opt/tljh/user/bin/activate nanopore
    minimap2 -t 16 -x map-ont polished/sample_${num}.racon1.fasta seed/sample_${num}.trimmed.fastq > seed/sample_${num}.seed2.paf

    echo "[*] Polishing with racon (round 2)"
    source /opt/tljh/user/bin/activate racon
    racon -t 16 seed/sample_${num}.trimmed.fastq seed/sample_${num}.seed2.paf polished/sample_${num}.racon1.fasta > polished/sample_${num}.racon2.fasta

    echo "[*] Final polishing with Medaka"
    source /opt/tljh/user/bin/activate medaka
    medaka_consensus -f -i seed/sample_${num}.trimmed.fastq -d polished/sample_${num}.racon2.fasta -o final/sample_${num} -t 16
    
    #rename file
    cp final/sample_${num}/consensus.fasta final/sample_${num}_final.fasta

    echo "✅ Finished sample ${num}"
done

echo "All samples complete."
