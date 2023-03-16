#!/bin/bash

# Usage: simulate_nanopore_mosquito.sh 

DIR=genomes/mosquito/

for DP in 1000000
do
    for FA in ${DIR}/deng*.fasta
    do
        echo $FA $DP
        ./pbsim --prefix sim_reads/mosquito/$(basename $FA .fasta)_${DP}X \
            --depth $DP \
            --difference-ratio 23:31:46 \
            --hmm_model ~/Code/pbsim2/data/R103.model \
            --length-max 20000 \
            --accuracy-mean 0.95 \
            $FA
    done
done
