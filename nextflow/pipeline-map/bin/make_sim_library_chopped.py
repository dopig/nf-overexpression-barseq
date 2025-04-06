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

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Generate an in silico library of barcoded plasmids')
    parser.add_argument('--unique-barcodes', type=int, default=int(1e6), help='Number of unique barcodes to generate (default: 1e6)')
    parser.add_argument('--barcode-length', type=int, default=20, help='Length of barcodes (default: 20)')
    parser.add_argument('--library-size', type=int, default=int(1e5), help='Final number of barcoded plasmids (default: 1e5)')
    parser.add_argument('--coverage', type=int,default=10, help='Sequencing coverage for mapping the library (default: 10). This is *NOT* the number of passes an individual amplicon gets in the sequencer; it is the number of identical amplicons that get sent to be sequenced. This generates a distribution (scale 1/10 of the median)around the default. For absolute coverage instead, use the flag --absolute_coverage.')
    parser.add_argument('--absolute-coverage', action='store_true', help='Use absolute coverage instead of a coverage distribution (default: False)')
    parser.add_argument('--map-reads', action='store_true', help='Include mapping of simulated reads (default: False)')
    parser.add_argument('--oligos', help=f'Path to oligo FASTA file')
    parser.add_argument('fasta', help='Input genome sequence in FASTA format (REQUIRED)')
    return parser.parse_args()

# def setup_reference_files(ref_output_dir: Path, args) -> Dict[str, Path]:
#     ref_output_dir.mkdir(parents=True, exist_ok=True)
#     genome_path = Path(args.fasta)
#     gff_path = genome_path.with_suffix(".gff")
#     feature_table_path = gff_path.parent / (gff_path.stem.rsplit('_',1)[0] + '_feature_table.txt')
#     oligos_path = Path(args.oligos)

#     abs_path_dict = {
#         "ref": ref_output_dir,
#         "fasta": ref_output_dir / genome_path.name,
#         "gff": ref_output_dir / gff_path.name,
#         "feature_table": ref_output_dir / feature_table_path.name,
#         "oligos": ref_output_dir / oligos_path.name,
#     }
#     try:
#         shutil.copyfile(genome_path, abs_path_dict["fasta"])
#         shutil.copyfile(gff_path, abs_path_dict["gff"])
#         shutil.copyfile(feature_table_path, abs_path_dict["feature_table"])
#         shutil.copyfile(oligos_path, abs_path_dict["oligos"])
#     except FileNotFoundError as e:
#         logging.exception(f"File not found: {e.filename}")
#         sys.exit(1)

#     return abs_path_dict

def get_genome(fasta_file: str) -> SeqIO.SeqRecord:
    try:
        if fasta_file.endswith(".gz"):
            with gzip.open(fasta_file, "rt") as handle:
                records = list(SeqIO.parse(handle, "fasta"))
        else:
            records = list(SeqIO.parse(fasta_file, "fasta"))
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
        # if direction == 1:
        #     insert_seq = genome_sequence[start:end] if start < end else genome_sequence[start:] + genome_sequence[:end]
        # else:
        #     insert_seq = genome_sequence[end:start][::-1] if end < start else (genome_sequence[:start] + genome_sequence[end:])[::-1]
        plasmid_list.append({"id": ix+1, "barcode_sequence": barcode, "direction": direction, "start": start, "end": end})
    return plasmid_list

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

def export_pcr_fasta(sequences: List[str], median_count: int, file_path: Path, absolute_coverage: bool = None) -> None:
    """
    Exports PCR sequences to a FASTA file.
    """
    scale = 0 if absolute_coverage else median_count/10
    num_digits = len(str(len(sequences)))
    count_generator = generate_random_integer(median = median_count, scale = scale, skewness = 0)
    with open(file_path, 'w') as f:
        for ix, seq in enumerate(sequences):
            seq_name = f"amplicon-{ix+1:0{num_digits}d}"
            seq_record = SeqRecord(Seq(seq), id=seq_name, name=seq_name, description='')
            for count in range(next(count_generator)):
                SeqIO.write(seq_record, f, 'fasta')
    logging.info(f"Exported {len(sequences)} PCR sequences, amplified ~{median_count}-fold")



def prepare_for_mapping(abs_map_dict: Dict[str, Path], output_dir: Path, map_temp_output_dir: Path, consensus_path: Path) -> None:
    map_temp_output_dir.mkdir(parents=True, exist_ok=False)
    shutil.copy(consensus_path, map_temp_output_dir)

    # Load the JSON file
    replacements = {
        'lib_names': [str(consensus_path.stem)],
        "lib_genome_dir": ".",
        'lib_genome_filenames': [str(abs_map_dict['fasta'].name)],
        'lib_genome_gffs': [str(abs_map_dict['gff'].name)],
    }

    with open(DEFAULT_BOBA_JSON_PATH, 'r') as f:
        json_data = json.load(f)

    for key, value in replacements.items():
        json_data[key] = value

    json_data['primer_info']['oligo_db_fp'] = str(BOBASEQ_WORK_DIR / abs_map_dict['oligos'].relative_to(output_dir))

    # Save the JSON file
    with open(abs_map_dict['ref'] / 'bobaseq_config.json', 'w') as f:
        json.dump(json_data, f, indent=4)



def main() -> None:
    output_dir = Path('.') # DATA_DIR / ('output/output-'+ datetime.now().strftime("%Y-%m-%d-%H%M%S"))
    ref_output_dir = output_dir / 'ref'
    map_temp_output_dir = output_dir / 'map_temp'


    setup_logging("Genome Chunking and Barcode Assignment", file_path=output_dir / 'log.txt')
    logging.info(f"Saving all output, including a more detailed version of this log, to {output_dir}/")

    raw_command_line = " ".join(sys.argv)
    logging.debug(f"Command run: {raw_command_line}")

    args = parse_args()
    log_args(vars(args))

    # Get genome here
    genome = get_genome(args.fasta)
    barcode_number_list = np.arange(args.unique_barcodes)
    sampled_barcode_indexes = np.random.choice(barcode_number_list, args.library_size, replace=True)
    quantify_barcode_redundancy(sampled_barcode_indexes)
    sampled_barcode_sequences = convert_numbers_to_sequences(sampled_barcode_indexes, args.barcode_length)

    plasmid_data = chunk_genome(len(genome), sampled_barcode_sequences)

    plasmid_path = output_dir / 'plasmids.json'
    with open(plasmid_path, 'w') as f:
        json.dump(plasmid_data, f, indent=4)
    logging.info(f"Exported {len(plasmid_data)} plasmid records to {plasmid_path.relative_to(output_dir)}")

    pcr_path = output_dir / 'pcrs.fasta'

    pcrs = generate_pcr_sequences(plasmid_data, genome)
    export_pcr_fasta(pcrs, median_count=args.coverage, file_path=pcr_path, absolute_coverage=args.absolute_coverage)

    # bam_file_stem = 'sequenced'
    # run_pbsim(pcr_path, output_dir=output_dir, output_prefix=bam_file_stem)
    # bam_path = output_dir / (bam_file_stem + '.bam')
    # consensus_file_path = generate_read_consensus(bam_path)

    # prepare_for_mapping(abs_map_dict=ref_path_dict, output_dir=output_dir, map_temp_output_dir=map_temp_output_dir, consensus_path=consensus_file_path)

    # if args.map_reads:
    #     map_reads(output_dir)
    # else:
    #     print('Skipping mapping reads. To map reads next time, use flag "--map-reads".')
    #     logging.info(f'To map these specific reads in the futue, run "src/map_reads.sh {output_dir}"')

if __name__ == '__main__':
    main()
