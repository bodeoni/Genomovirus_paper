#!bin/bash

exec > >(tee log/assembly.log) 2>&1

# Flye assembly
source /opt/tljh/user/bin/activate flye-env

for num in 52 54 60; do 
    echo "Starting Flye assembly for barcode ${num}"
    mkdir -p assembly/flye/barcode${num}
    flye --nano-raw fastq/sample_*_${num}.fastq \
        --out-dir assembly/flye/barcode${num} \
        --threads 16 --genome-size 2.5k
    echo "Flye assembly for barcode ${num} completed"
done

# Canu assembly
source /opt/tljh/user/bin/activate assembly

for num in 52 54 60; do
    echo "Starting Canu assembly for barcode ${num}"
    mkdir -p assembly/canu/barcode${num}
    canu -p barcode${num} -d assembly/canu/barcode${num} \
        genomeSize=2.5k -nanopore fastq/sample_*_${num}.fastq \
        maxThreads=16
    echo "Canu assembly for barcode ${num} completed"
done

# raven assembly
source /opt/tljh/user/bin/activate raven

for num in 52 54 60; do
    echo "Starting Raven assembly for barcode ${num}"
    mkdir -p assembly/raven/barcode${num}
    raven -t 16 fastq/sample_*_${num}.fastq old_run/barcode${num}.gemy.fastq > assembly/raven/barcode${num}/assembly.fasta
    echo "Raven assembly for barcode ${num} completed"
done
