#!/bin/bash

# Usage: simulate_nanopore.sh <DIR_NAME>

DIR=$1

for DP in 1 3 5 7 10 15 20 25 30 35 40 45 50
do
    for FA in genomes/${DIR}_mut/*.fasta
    do
        echo $FA $DP
        ./pbsim --prefix sim_reads/${DIR}_${DP}X/$(basename $FA .fasta)_${DP}X \
            --depth $DP \
            --difference-ratio 23:31:46 \
            --hmm_model ~/Code/pbsim2/data/R103.model \
            --length-max 20000 \
            --accuracy-mean 0.95 \
            $FA
    done
done
