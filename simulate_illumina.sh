#!/bin/bash

# Usage: simulate_illumina.sh <FASTA_DIR> <OUT_DIR>

DIR=$1
OUT=$2

for DP in 1 3 5 7 10 20 40
do
    for FA in $DIR/*.fasta
    do
        echo $FA $DP
        rm ./cov/cov.txt
        for SQ in $(grep ">" $FA | perl -pe 's/>//; s/\s+.*//')
        do
            echo -e "$SQ\t$DP" >> ./cov/cov.txt
        done
        iss \
            generate \
            -p 4 \
            -g $FA \
            -z \
            -m MiSeq \
            -D ./cov/cov.txt \
            -o $OUT/$(basename $FA .fasta)_${DP}X
    done
done
