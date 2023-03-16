#!/bin/bash

HOST_FQ='sim_illumina/mosquito/mosquito_genome_1M'
PATHO_FQ='sim_illumina/mosquito/virus_1M'
OUTDIR='sim_illumina/mosquito/'

for FR in 0.005 0.01 0.05 0.1 0.2 0.3 0.4
do
    for R in 1 2
    do
        python3 spike_reads.py \
            -t "${HOST_FQ}_R$R.fastq.gz" \
            -p "${PATHO_FQ}_R$R.fastq.gz" \
            -o "$OUTDIR/host_spiked_${FR}_R$R.fastq.gz" \
            -f $FR \
            -v
    done
done