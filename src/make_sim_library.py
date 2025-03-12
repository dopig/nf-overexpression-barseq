#!/usr/bin/env python3

import argparse
import random
import json
import subprocess
import gzip
import shutil
import sys
from typing import Generator, List, Dict
from pathlib import Path
from datetime import datetime

import numpy as np
from scipy.stats import skewnorm
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from utils import setup_logging, loginfo, logerror


SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR = SCRIPT_DIR.parent
DATA_DIR = ROOT_DIR / 'data'
DEFAULT_BOBA_JSON_PATH = ROOT_DIR / 'config/default_bobaseq_config.json'

BOBASEQ_WORK_DIR = Path('/work') # Directory within Docker container


def parse_args():
    parser = argparse.ArgumentParser(description='Generate an in silico library of barcoded plasmids')
    parser.add_argument('--unique_barcodes', type=int, default=int(1e6), help='Number of unique barcodes to generate (default: 1e6)')
    parser.add_argument('--barcode_length', type=int, default=20, help='Length of barcodes (default: 20)')
    parser.add_argument('--library_size', type=int, default=int(1e5), help='Final number of barcoded plasmids (default: 1e5)')
    parser.add_argument('--coverage', type=int,default=10, help='Sequencing coverage for mapping the library (default: 10). This is *NOT* the number of passes an individual amplicon gets in the sequencer; it is the number of identical amplicons that get sent to be sequenced. This generates a distribution (scale 1/10 of the median)around the default. For absolute coverage instead, use the flag --absolute_coverage.')
    parser.add_argument('--absolute_coverage', action='store_true', help='Use absolute coverage instead of a coverage distribution (default: False)')
    default_oligos_path = f'{DATA_DIR}/reference/oligos.fasta'
    parser.add_argument('--oligos', default=default_oligos_path, help=f'Path to oligo FASTA file (default: {default_oligos_path})')
    parser.add_argument('fasta', help='Input genome sequence in FASTA format (REQUIRED)')
    return parser.parse_args()

def setup_reference_files(ref_output_dir: Path, args) -> Dict[str, Path]:
    ref_output_dir.mkdir(parents=True, exist_ok=True)
    genome_path = Path(args.fasta)
    gff_path = Path(args.fasta).with_suffix(".gff")
    oligos_path = Path(args.oligos)

    abs_path_dict = {
        "ref": ref_output_dir,
        "fasta": ref_output_dir / genome_path.name,
        "gff": ref_output_dir / gff_path.name,
        "oligos": ref_output_dir / oligos_path.name,
    }
    shutil.copyfile(genome_path, abs_path_dict["fasta"])
    shutil.copyfile(gff_path, abs_path_dict["gff"])
    shutil.copyfile(oligos_path, abs_path_dict["oligos"])

    return abs_path_dict

def validate_genome(fasta_file: str) -> SeqIO.SeqRecord:
    records = list(SeqIO.parse(fasta_file, "fasta"))
    if len(records) != 1:
        loginfo("Error: More than one contig detected. Only single-contig genomes are supported.", True)
        return None
    return records[0]

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
    loginfo("Barcode Redundancy Distribution:", True)
    for copies, num_barcodes in enumerate(redundancy_distribution):
        if num_barcodes > 0:
            loginfo(f'Number of barcodes with {copies} copies: {num_barcodes}', True)

def random_dna(length: int) -> str:
    return ''.join(random.choices('ATCG', k=length))

def convert_numbers_to_sequences(sampled_barcodes: np.ndarray, barcode_length: int) -> List[str]:
    """
    Converts barcode numbers into nucleotide sequences.
    """
    barcode_dict = {num: random_dna(barcode_length) for num in set(sampled_barcodes)}
    return [barcode_dict[num] for num in sampled_barcodes]

def chunk_genome(genome_sequence: str, sampled_barcodes: np.ndarray) -> List[Dict[str, object]]:
    genome_length = len(genome_sequence)
    chunk_size_generator = generate_random_integer(median = 2500, scale = 2000, skewness = 10)
    plasmid_list = []
    for barcode in sampled_barcodes:
        start = random.randint(0, genome_length - 1)
        chunk_size = next(chunk_size_generator)
        direction = random.choice([1, -1])
        end = (start + direction * chunk_size) % genome_length
        if direction == 1:
            insert_seq = genome_sequence[start:end] if start < end else genome_sequence[start:] + genome_sequence[:end]
        else:
            insert_seq = genome_sequence[end:start][::-1] if end < start else (genome_sequence[:start] + genome_sequence[end:])[::-1]
        plasmid_list.append({"barcode_sequence": barcode, "insert_sequence": str(insert_seq.seq), "start": start, "end": end})
    return plasmid_list

def generate_pcr_sequences(plasmids: List[Dict[str, object]]) -> List[str]:
    """
    Generates PCR sequences using defined flanking regions and barcodes.
    Outputs have constant fwd/rev indexes for further processing for now.
    """
    ix5 = "GTTCTTATCTTTGCAGTCTC"
    left = "TGTTGACAATTAATCATCCGGCTCGTATAATGTGTGGAATTGTGAGCGGATAACAATTTCAGAATTCAC"
    mid = "GTGAGCCTCGGTACCAAATTCCAGAAAAGAGGCCTCCCGAAAGGGGGGCCTTTTTTCGTTTTGGTCCGCCTCGTGAGAGCTGGTCGACCTGCAGCGTACG"
    right = "AGAGACCTCGTGGACATCTATCAGAGACTATCAGTTTTTTTGATTTCTTCCCTTGCCTTGTCAATCCTTGCTTGCAGCTCCGGGGTTATCATCAAATCTTCACGACCAACTTTTACCAAAGCGTAAATCTC"
    ix3 = "GAGATTTACGCTTTGGTAAAAGTTGG"
    sequences = [ix5 + left + plasmid['insert_sequence'] + mid + plasmid['barcode_sequence'] + right + ix3 for plasmid in plasmids]
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
    loginfo(f"Exported {len(sequences)} PCR sequences, amplified ~{median_count}-fold, to {file_path}", True)

def run_pbsim(pcr_file_path: Path, output_dir: Path, output_prefix: str = None) -> None:
    """
    Runs pbsim on the PCR sequences.  Currently does 25 passes and outputs a bam and maf.gz file.
    """
    # Construct absolute path to QSHMM-RSII.model
    qshmm_model_path_abs = (DATA_DIR / "reference/pbsim/QSHMM-RSII.model").resolve()
    pcr_file_path_abs = pcr_file_path.resolve()

    command = [
        "pbsim",
        "--strategy", "templ",
        "--method", "qshmm",
        "--pass-num", "25",
        "--qshmm", str(qshmm_model_path_abs),  # Convert Path to string
        "--template", str(pcr_file_path_abs),    # Convert Path to string
    ]

    if output_prefix:
        command.extend(["--prefix", output_prefix])

    loginfo(f"Running command: {' '.join(command)}", True)

    # Run pbsim command with subprocess.run and capture output
    # Merge stderr because all the information goes to stderr and stdout is empty
    result = subprocess.run(
        command,
        cwd=output_dir,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,  # Merge stderr into stdout
    )

    # Log pbsim output
    if result.stdout:
        lines = result.stdout.splitlines()
        indented_output = "\n".join([" " * 11 + ":" * 8 + ' - ' + line for line in lines])
        loginfo(f"PBSIM OUTPUT:\n{indented_output}")

    # Check the return code
    if result.returncode != 0:
        logerror(f"PBSIM exited with return code: {result.returncode}")

    loginfo(f"PBSIM completed", True)

def generate_read_consensus(bam_path: Path) -> Path:
    consensus_path = bam_path.parent / "consensus.fastq"
    command = [
        SCRIPT_DIR/"run_pbccs.sh",
        bam_path.parent,
        bam_path.name,
        consensus_path.name
    ]

    loginfo(f"Running command: {' '.join([str(x) for x in command])}")

    result = subprocess.run(
        command,
        cwd=SCRIPT_DIR,
        capture_output=True,
        text=True
    )

    if result.returncode != 0:
        indented_output = "\n".join([" " * 11 + ":" * 8 + ' - ' + line for line in result.stderr.splitlines() if line.strip() != ''])
        message = f"PBCCS exited with return code: {result.returncode}\n{indented_output}"
        logerror(message, print_out=True)
    else:
        loginfo(f"PBCCS completed", True)
        return consensus_path

def prepare_for_mapping(abs_map_dict: Dict[str, Path], output_dir: Path, map_output_dir: Path, consensus_path: Path) -> None:
    map_output_dir.mkdir(parents=True, exist_ok=False)
    shutil.copy(consensus_path, map_output_dir)

    # Load the JSON file
    replacements = {
        'lib_names': [str(consensus_path.stem)],
        'lib_genome_dir': str(BOBASEQ_WORK_DIR / abs_map_dict['fasta'].relative_to(output_dir).parent),
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
    output_dir = DATA_DIR / ('output/output-'+ datetime.now().strftime("%Y-%m-%d-%H%M%S"))
    ref_output_dir = output_dir / 'ref'
    sim_output_dir = output_dir / 'sim'
    map_output_dir = output_dir / 'map'

    setup_logging("Genome Chunking and Barcode Assignment", file_path=output_dir / 'log.txt')
    loginfo(f"Saving all output to {output_dir}/", True)

    raw_command_line = " ".join(sys.argv)
    loginfo(f"Command run: {raw_command_line}", print_out=False)

    args = parse_args()
    loginfo(f"Running with arguments: {vars(args)}")

    ref_path_dict = setup_reference_files(ref_output_dir,args)
    print(ref_path_dict)

    genome = validate_genome(ref_path_dict["fasta"])
    barcode_number_list = np.arange(args.unique_barcodes)
    sampled_barcode_indexes = np.random.choice(barcode_number_list, args.library_size, replace=True)
    quantify_barcode_redundancy(sampled_barcode_indexes)
    sampled_barcode_sequences = convert_numbers_to_sequences(sampled_barcode_indexes, args.barcode_length)
    loginfo("Generating inserts from genome")
    plasmid_data = chunk_genome(genome, sampled_barcode_sequences)

    sim_output_dir.mkdir(parents=True, exist_ok=True)
    plasmid_path = sim_output_dir / 'plasmids.json'
    with open(plasmid_path, 'w') as f:
        json.dump(plasmid_data, f, indent=4)
    loginfo(f"Exported {len(plasmid_data)} plasmid records to {plasmid_path}", True)

    pcr_path = sim_output_dir / 'pcrs.fasta'
    pcrs = generate_pcr_sequences(plasmid_data)
    export_pcr_fasta(pcrs, median_count=args.coverage, file_path=pcr_path, absolute_coverage=args.absolute_coverage)

    bam_file_stem = 'sequenced'
    run_pbsim(pcr_path, output_dir=sim_output_dir, output_prefix=bam_file_stem)
    bam_path = sim_output_dir / (bam_file_stem + '.bam')
    consensus_file_path = generate_read_consensus(bam_path)

    prepare_for_mapping(abs_map_dict=ref_path_dict, output_dir=output_dir, map_output_dir=map_output_dir, consensus_path=consensus_file_path)


if __name__ == '__main__':
    main()
