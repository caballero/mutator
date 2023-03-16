#!/bin/bash

# Usage: simulate_illumina_mosquito.sh <FASTA_DIR> <OUT_DIR>

DIR=$1
OUT=$2

for NR in 2M
do
    for FA in $DIR/*.fasta
    do
        echo $FA $DP
        iss \
            generate \
            -g $FA \
            -z \
            -m MiSeq \
            -n $NR \
            -o $OUT/$(basename $FA .fasta)_$NR
    done
done
