#!/usr/bin/env python3

import argparse
import random
import json
import logging
import shutil
import sys
import gzip
from typing import Generator, List, Dict
from pathlib import Path

import numpy as np
from scipy.stats import skewnorm
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from utils import setup_logging, format_output, log_args

OUTPUT_DIR = Path('.')

def begin_log() -> None:
    """
    Sets up logging and saves initial information to the log file.
    """
    setup_logging("Genome Chunking and Barcode Assignment", file_path=OUTPUT_DIR / 'log.txt')
    logging.info(f"Saving all output, including a more detailed version of this log, to {OUTPUT_DIR}/")
    raw_command_line = " ".join(sys.argv)
    logging.debug(f"Command run: {raw_command_line}")

def parse_args() -> argparse.Namespace:
    """
    Parse the command line arguments to generate an in silico library of barcoded plasmids.
    """
    parser = argparse.ArgumentParser(description='Generate an in silico library of barcoded plasmids')
    parser.add_argument('--unique-barcodes', type=int, default=int(1e6), help='Number of unique barcodes to generate (default: 1e6)')
    parser.add_argument('--barcode-length', type=int, default=20, help='Length of barcodes (default: 20)')
    parser.add_argument('--library-size', type=int, default=int(1e5), help='Final number of barcoded plasmids (default: 1e5)')
    parser.add_argument('--coverage', type=int,default=10, help='Sequencing coverage for mapping the library (default: 10). This is *NOT* the number of passes an individual amplicon gets in the sequencer; it is the number of identical amplicons that get sent to be sequenced. This generates a distribution (scale 1/10 of the median)around the default. For absolute coverage instead, use the flag --absolute_coverage.')
    parser.add_argument('--absolute-coverage', action='store_true', help='Use absolute coverage instead of a coverage distribution (default: False)')
    parser.add_argument('--map-reads', action='store_true', help='Include mapping of simulated reads (default: False)')
    parser.add_argument('--oligos', type=Path, help=f'Path to oligo FASTA file')
    parser.add_argument('--random-seed', type=int, default=None, help='Number to set random seed to (default: None)')
    parser.add_argument('fasta', type=Path, help='Input genome sequence in FASTA format (REQUIRED)')
    args = parser.parse_args()
    log_args(vars(args))
    return args

def get_genome(fasta_path: Path) -> SeqIO.SeqRecord:
    """
    Parses a genome sequence from a given FASTA file.

    This function supports both plain and gzipped FASTA files and extracts the
    genome sequence as a SeqRecord object. It ensures that the input file
    contains only a single contig; otherwise, it raises an error.

    Args:
        fasta_path (Path): Path to the input genome sequence in FASTA format.
                           Can be a plain text or gzipped file.

    Returns:
        SeqIO.SeqRecord: The genome sequence as a SeqRecord object.
    """

    try:
        if fasta_path.suffix == ".gz":
            with gzip.open(fasta_path, "rt") as handle:
                records = list(SeqIO.parse(handle, "fasta"))
        else:
            records = list(SeqIO.parse(fasta_path, "fasta"))
        if len(records) != 1:
            raise ValueError("More than one contig detected. Only single-contig genomes are currently supported.")
        return records[0]
    except ValueError as e:
        logging.exception(f"Invalid genome file: {e}")
        sys.exit(1)
    except Exception as e:
        logging.exception(f"Error parsing genome file: {e}")
        sys.exit(1)

def generate_random_integer(median: int, scale: float, skewness: int) -> Generator[int, None, None]:
    """
    Generator to get appropriate distribution of sizes, e.g. for inserts

    Args:
        median (int): The median value of the distribution.
        scale (int): The scale parameter (~standard deviation) of the distribution.
        skewness (int): The skewness parameter of the distribution.
    """
    while True:
        random_number = skewnorm.rvs(a=skewness, loc=median, scale=scale, size=1)[0]
        yield int(random_number)

def design_plasmids(unique_barcodes, library_size, barcode_length, genome_length):

    def quantify_barcode_redundancy(barcodes: np.ndarray) -> None:
        """
        Prints barcode redundancy summary using NumPy.
        """
        unique_elements, counts = np.unique(barcodes, return_counts=True)
        redundancy_distribution = np.bincount(counts)
        logging.info("Barcode Redundancy Distribution:")
        for copies, num_barcodes in enumerate(redundancy_distribution):
            if num_barcodes > 0:
                logging.info(f'Number of barcodes with {copies} copies: {num_barcodes}')

    def random_dna(length: int) -> str:
        return ''.join(random.choices('ATCG', k=length))

    def convert_numbers_to_sequences(sampled_barcodes: np.ndarray, barcode_length: int) -> List[str]:
        """
        Converts barcode numbers into nucleotide sequences.
        """
        barcode_dict = {num: random_dna(barcode_length) for num in set(sampled_barcodes)}
        return [barcode_dict[num] for num in sampled_barcodes]

    def chunk_genome(genome_length: int, sampled_barcodes: np.ndarray) -> List[Dict[str, object]]:
        chunk_size_generator = generate_random_integer(median = 2500, scale = 2000, skewness = 10)
        plasmid_list = []
        for ix, barcode in enumerate(sampled_barcodes):
            chunk_size = next(chunk_size_generator)
            direction = random.choice(['+', '-'])
            start = random.randint(0, genome_length - 1)
            end = (start + chunk_size) % genome_length
            plasmid_list.append({"id": ix+1, "barcode_sequence": barcode, "direction": direction, "start": start, "end": end})
        return plasmid_list

    barcode_number_list = np.arange(unique_barcodes)
    sampled_barcode_indexes = np.random.choice(barcode_number_list, library_size, replace=True)
    quantify_barcode_redundancy(sampled_barcode_indexes)
    sampled_barcode_sequences = convert_numbers_to_sequences(sampled_barcode_indexes, barcode_length)
    plasmid_data = chunk_genome(genome_length, sampled_barcode_sequences)
    return plasmid_data

def export_plasmids(plasmids: List[Dict[str, object]], output_path: Path) -> None:
    with open(output_path, 'w') as f:
        json.dump(plasmids, f, indent=4)
    logging.info(f"Exported {len(plasmids)} plasmid records to {output_path.relative_to(OUTPUT_DIR)}")

def generate_pcr_sequences(plasmids: List[Dict[str, object]], genome: SeqRecord) -> List[str]:
    """
    Generates PCR sequences using defined flanking regions and barcodes.
    Outputs have constant fwd/rev indexes for further processing for now.
    Order is ix5-left-insert-mid-barcode-right-ix3
    """
    ix5 = "GTTCTTATCTTTGCAGTCTC"
    left = "TGTTGACAATTAATCATCCGGCTCGTATAATGTGTGGAATTGTGAGCGGATAACAATTTCAGAATTCAC"
    mid = "GTGAGCCTCGGTACCAAATTCCAGAAAAGAGGCCTCCCGAAAGGGGGGCCTTTTTTCGTTTTGGTCCGCCTCGTGAGAGCTGGTCGACCTGCAGCGTACG"
    right = "AGAGACCTCGTGGACATCTATCAGAGACTATCAGTTTTTTTGATTTCTTCCCTTGCCTTGTCAATCCTTGCTTGCAGCTCCGGGGTTATCATCAAATCTTCACGACCAACTTTTACCAAAGCGTAAATCTC"
    ix3 = "GAGATTTACGCTTTGGTAAAAGTTGG"

    genome_sequence = genome.seq

    sequences = []

    for plasmid in plasmids:
        start, end = plasmid['start'], plasmid['end']
        insert_seq = genome_sequence[start:end] if start < end else genome_sequence[start:] + genome_sequence[end]
        insert_seq = insert_seq if plasmid['direction'] == '+' else insert_seq.reverse_complement()
        sequences.append(ix5 + left + insert_seq + mid + plasmid['barcode_sequence'] + right + ix3)

    return sequences

def export_pcr_fasta(sequences: List[str], median_count: int, output_path: Path, absolute_coverage: bool = None) -> None:
    """
    Exports PCR sequences to a FASTA file.
    """
    scale = 0 if absolute_coverage else median_count/10
    num_digits = len(str(len(sequences)))
    count_generator = generate_random_integer(median = median_count, scale = scale, skewness = 0)
    with open(output_path, 'w') as f:
        for ix, seq in enumerate(sequences):
            seq_name = f"amplicon-{ix+1:0{num_digits}d}"
            seq_record = SeqRecord(Seq(seq), id=seq_name, name=seq_name, description='')
            for count in range(next(count_generator)):
                SeqIO.write(seq_record, f, 'fasta')
    logging.info(f"Exported {len(sequences)} PCR sequences, amplified ~{median_count}-fold")

def main() -> None:
    begin_log()
    args = parse_args()

    if args.random_seed is not None:
        random.seed(args.random_seed)

    genome = get_genome(args.fasta)

    plasmids = design_plasmids(args.unique_barcodes, args.library_size, args.barcode_length, len(genome))
    export_plasmids(plasmids=plasmids, output_path=OUTPUT_DIR / 'plasmids.json')

    pcrs = generate_pcr_sequences(plasmids, genome)
    export_pcr_fasta(
        pcrs,
        median_count=args.coverage,
        output_path=OUTPUT_DIR / 'pcrs.fasta',
        absolute_coverage=args.absolute_coverage
    )

if __name__ == '__main__':
    main()
