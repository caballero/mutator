#!/usr/bin/env python3

# mutator.py - a script to mutate sequences (single strand genomes)
# (C) 2023

import argparse
import logging
from random import randint

def load_sequence(file):
    logging.debug(f"Loading sequence from {file}")
    name = ''
    seq  = ''
    with open(file, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                name = line.rstrip()
            else:
                seq = seq + line.rstrip().upper()
    logging.debug("Read " + str(len(seq)) + " bases")
    return name, seq

def mutate(base):
    if base == 'A':
        opt = ['C', 'G', 'T', 'T']
    elif base == 'T':
        opt = ['C', 'G', 'A', 'A']
    elif base == 'C':
        opt = ['A', 'T', 'G', 'G']
    elif base == 'G':
        opt = ['A', 'T', 'C', 'C']
    else:
        opt = ['A', 'T', 'C', 'G']
    return opt[ randint(0,3) ]

def add_mutations(seq, mut):
    logging.debug(f"Adding {mut} mutations")
    aseq = list(seq)
    for n in range(mut):
        pos = randint(0, len(aseq) -1)
        old = aseq[pos]
        new = mutate(base=old)
        logging.debug(f"  iter {n}: changing {pos} ({old} -> {new})")
        aseq[pos] = new
    return ''.join(aseq)

def write_seq(name, seq, file):
    with open(file, "wt") as fh:
        fh.write(name + "\n" + seq + "\n")
        

def main():
    rate = 0.01
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fasta", help="Path to Fasta file", required=True)
    parser.add_argument("-o", "--output", help="Mutated Fasta output file", required=True)
    parser.add_argument("-r", "--rate", type=float, help=f"Mutation rate [0.0 - 1.0], default = {rate}", required=False)
    parser.add_argument("-v", "--verbose", help="Verbose mode", action="store_true")
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    
    if args.rate:
        rate = args.rate

    seq_name, seq = load_sequence(file=args.fasta)
    mutations = int(rate * len(seq))
    logging.debug(f"Total mutations = {mutations} with rate = {rate}")
    seq_mut = add_mutations(seq=seq, mut=mutations)
    write_seq(name=seq_name, seq=seq_mut, file=args.output)

    logging.debug("Completed")

if __name__ == "__main__":
    main()
