#!/usr/bin/env python3

# spike_reads.py - a script to spike nanopore reads in a fastq file

import argparse
import logging
import gzip
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def main():
    fraction = 0.01
    host_reads = 1000000
    pathogen_reads = 1000000
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--host", help="Path to host Fastq file", required=True)
    parser.add_argument("-p", "--pathogen", help="Path to pathogen Fastq file", required=True)
    parser.add_argument("-o", "--output", help="Fastq output file", required=True)
    parser.add_argument("-f", "--fraction", type=float, help=f"Fraction to spike [0.0 - 1.0], default = {fraction}", required=False)
    parser.add_argument("--host-reads", type=int, help=f"Total reads for host fastq, default = {host_reads}", required=False)
    parser.add_argument("--pathogen-reads", type=int, help=f"Total reads for pathogen fastq, default = {pathogen_reads}", required=False)
    parser.add_argument("-v", "--verbose", help="Verbose mode", action="store_true")
    args = parser.parse_args()
 
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    
    if args.fraction:
        fraction = args.fraction
    if args.host_reads:
        host_reads = args.host_reads
    if args.pathogen_reads:
        pathogen_reads = args.pathogen_reads
    
    host_fastq = args.host
    pathogen_fastq = args.pathogen
    out_fastq = args.output

    total_host_reads = int((1 - fraction) * host_reads)
    total_pathogen_reads = host_reads - total_host_reads

    logging.debug(f"Host file: {host_fastq}")
    logging.debug(f"Pathogen file: {pathogen_fastq}")
    logging.debug(f"Output file: {out_fastq}")
    logging.debug(f"Host reads: {host_reads}")
    logging.debug(f"Pathogen reads: {pathogen_reads}")
    logging.debug(f"Reads fraction: {fraction}")
    logging.debug(f"Total host reads: {total_host_reads}")
    logging.debug(f"Total pathogen reads: {total_pathogen_reads}")

    with gzip.open(out_fastq, "wt") as output_fh:
        counter = 0
        logging.debug("Adding host reads")
        with gzip.open(host_fastq, "rt") as host_fh:
            for record in SeqIO.parse(host_fh, "fastq"):
                new_id = record.id.replace("S", "H")
                new_record = SeqRecord(
                            record.seq,
                            id=new_id,
                            name=new_id,
                            description='',
                            letter_annotations=record.letter_annotations
                        )

                SeqIO.write(new_record, output_fh, "fastq")
                counter += 1
                if counter >= total_host_reads:
                    break
        logging.debug(f" {counter} reads added")

        logging.debug("Adding pathogen reads")
        with gzip.open(pathogen_fastq, "rt") as pathogen_fh:
            for record in SeqIO.parse(pathogen_fh, "fastq"):
                new_id = record.id.replace("S", "P")
                new_record = SeqRecord(
                            record.seq,
                            id=new_id,
                            name=new_id,
                            description='',
                            letter_annotations=record.letter_annotations
                        )

                SeqIO.write(new_record, output_fh, "fastq")
                counter += 1
                if counter >= host_reads:
                    break
        logging.debug(f" {counter} reads added")

    logging.debug("Completed")

if __name__ == "__main__":
    main()
