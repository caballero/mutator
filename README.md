# Mutator

Simple script to add mutations in a genomic sequence (single chromosome, i.e. virus or bacteria).

## Requirement

- Python3

## Usage

```bash
usage: mutator.py [-h] -f FASTA -o OUTPUT [-r RATE] [-v]

options:
  -h, --help            show this help message and exit
  -f FASTA, --fasta FASTA
                        Path to Fasta file
  -o OUTPUT, --output OUTPUT
                        Mutated Fasta output file
  -r RATE, --rate RATE  Mutation rate [0.0 - 1.0], default = 0.01
  -v, --verbose         Verbose mode
```

(C) 2023
