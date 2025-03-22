#!/usr/bin/env python3


"""
This script processes plasmid and CDS feature data to simulate barseq reads. It selects "winner" CDS features,
determines which plasmids contain them, and simulates barseq reads with variations in read counts.

Input Files:
    - [working_dir]/library/plasmids.json: Plasmid information.
    - [working_dir]/ref/*.gff: GFF file containing CDS feature data.
    - [working_dir]/reference/bobaseq_barseq_samples.tsv: Sample information.
    - [working_dir]/reference/barseq4.index2: Multiplex primer index information.

Output Files:
    - [working_dir]/barseq/feature_to_plasmid.json: Map of CDS features to plasmids.
    - [working_dir]/barseq/reads.fastq: Simulated barseq reads in FASTQ format.
    - [working_dir]/log.txt: log file of the run.

Parameters:
    --working-dir: Path to the working directory.
    --samples-tsv-path: Path to the samples TSV file.
        If none provided, [working_dir]/reference/bobaseq_barseq_samples.tsv will be used.
    --multiplex-index-tsv: Path to the multiplex index TSV file.
        If none provided, [working_dir]/reference/barseq4.index2 will be used.
    --winner-count: Number of "winner" CDS features to select.
    --winner-strength: Increase in read count for "winner" features.
    --count-range: Range of read counts for simulation.
    --plusminus: Variation in read counts.
    --quality-score: Quality score for FASTQ output.

Example Usage:
    python3 simulate_barseq_reads.py --working-dir /path/to/working_dir --winner-count 30
"""




"""
Point this script (with --working-dir) to the directory made by the make_sim_library script.
Within that working directory it will use the following files produced by the make_sim_library.py script:
- [working_dir]/library/plasmids.json
- [working_dir]/ref/*.gff (needs to match only one GFF file containing only one contig)

You will also need a few more files. See the examples in data/reference for more details:
- samples_tsv_path - You can provide path as an argument or it will use a the project-level
        data/reference/bobaseq_barseq_samples.tsv; whichever is chosen, is copied to
        [working_dir]/reference/bobaseq_barseq_samples.tsv and used from there
- multiplex_index_tsv - You can provide path as an argument or it will use a the project-level
        data/reference/barseq4.index2; whichever is chosen, is copied to
        [working_dir]/reference/barseq4.index2 and used from there

This script will select random "winner" CDS features from the GFF and find which plasmids contain those features.
"Loser" CDS features are not a thing at this time. A JSON map from feature to plasmid will be created and exported.
It will then simulate barseq reads from those plasmids as follows:
- Establish a primary distribution for all barcodes: these will be in range `count_range`
- All samples begin with slight variations, within `plusminus` distance, of the primary distribution
- Non-Time0 barcodes will be shifted by `winner_strength` for the "winning" features
- The reads are exported to a fastq file
"""




import argparse
import json
import random
import logging
import sys
from pathlib import Path
from pprint import pformat
from shutil import copyfile
from typing import List, Dict, Tuple, Union

import gffutils
import pandas as pd

from utils import setup_logging, format_output, log_args


SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR = SCRIPT_DIR.parent
DATA_DIR = ROOT_DIR / 'data'
TIME0_NAME = 'Time0'


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('working_dir', help='Working directory which has library, map, and ref subdirectories')
    parser.add_argument('--samples-tsv-path', default=f'{DATA_DIR}/reference/bobaseq_barseq_samples.tsv',
        help='Path to samples TSV file. This script only needs the index and Group columns. If no file path is provided, data/reference/bobaseq_barseq_samples.tsv will be used')
    parser.add_argument('--multiplex-index-tsv', default=f'{DATA_DIR}/reference/barseq4.index2',
        help='Path to multiplex primer barseq4 TSV file. If no file path is provided, data/reference/barseq4.index2 will be used')
    parser.add_argument('--winner-count', type=int, default=20, help='How many will be selected as winners (default: 20)')
    parser.add_argument('--winner-strength', type=int, default=5000, help='How many counts stronger winners will be (default: 5000)')
    parser.add_argument('--count-range', type=int, nargs=2, default=[0, 1000], help='Count range (min, max), entered as "--count-range min max" (default: 0 1000)')
    parser.add_argument('--plusminus', type=int, default=1, help='Plus/minus value (default: 1)')
    parser.add_argument('--quality-score', type=int, default=30, help='Quality score to be used for all bases (default: 30)')
    return parser.parse_args()

def get_samples_df(samples_path: Path, index_path: Path) -> pd.DataFrame:
    """
    Read in samples TSV file and Barseq index TSV file.

    Returns a pandas dataframe with columns 'Group', 'index_name', 'index',
    and 'name'. The 'index' column is the actual index sequence.
    """
    df_samples = pd.read_csv(samples_path, sep='\t', usecols=['index', 'Group'])
    df_samples.rename(columns={'index': 'index_name'}, inplace=True)

    df_barseq4 = pd.read_csv(index_path, sep='\t')
    df_samples = pd.merge(df_samples, df_barseq4, on='index_name')

    return df_samples

def extract_cds_from_gff(gff_path: Path) -> List[Tuple[str, int, int]]:
    """
    Extract CDS features from a GFF file.

    :param gff_path: Path to a GFF file
    :return: A list of tuples, where each tuple contains the id, start, and end of a CDS feature
    """
    db = gffutils.create_db(gff_path, ":memory:")
    count = db.count_features_of_type('CDS')
    location_list = [(feat.id, feat.start, feat.end) for feat in db2.all_features(featuretype='CDS')]
    db.close()
    return records

def write_feature_to_plasmid_json(feature_to_plasmid_dict: Dict[str, List[str]], output_file_path: Path) -> None:
    """
    Write a JSON file mapping feature IDs to lists of overlapping plasmid IDs.

    :param feature_to_plasmid_dict: A dict mapping feature IDs to lists of overlapping plasmid IDs
    :param output_file_path: The path to which the JSON file should be written
    """
    output_file_path.parent.mkdir(exist_ok=True)
    with open(output_file_path, 'w') as f:
        json.dump(feature_to_plasmid_dict, f, indent=4)
    logging.debug(f'Writing feature-to-plasmid JSON for {len(feature_to_plasmid_dict)} features to {output_file_path}')

def find_overlapping_plasmids(
    plasmids: List[Dict[str, Union[int, str]]],
    cds_regions: List[Tuple[str, int, int]],
    feature_to_plasmid_json_path: Path
) -> tuple[set, Dict[str, List[str]]]:
    """
    Given a list of plasmids and a list of CDS regions, return a set of all
    plasmid IDs that overlap with any of the CDS regions and a dict mapping
    each CDS feature ID to a list of overlapping plasmid IDs.
    """
    all_matching_plasmid_ids = set()
    feat_to_plasmid_dict = dict()
    for feature_id, cds_start, cds_end in cds_regions:
        feature_matching_plasmid_ids = []
        for plasmid in plasmids:
            # Cover normal case (first line below) and then if plasmid spans the origin
            if (plasmid['start'] <= cds_start and plasmid['end'] >= cds_end) or \
               (plasmid['start'] > plasmid['end'] and (cds_end < plasmid['end'] or cds_start > plasmid['start'] )):
                feature_matching_plasmid_ids.append(plasmid['id'])
        if feature_matching_plasmid_ids:
            feat_to_plasmid_dict[feature_id] = feature_matching_plasmid_ids
            all_matching_plasmid_ids.update(feature_matching_plasmid_ids)
        else:
            feat_to_plasmid_dict[feature_id] = None

    if not feat_to_plasmid_dict:
        logging.error(f"None of the {len(plasmids)} plasmids contained any of the {len(cds_regions)} CDSs")
        sys.exit(1)

    found_feature_count = sum(1 for value in feat_to_plasmid_dict.values() if value is not None)
    logging.info(f"Found {len(all_matching_plasmid_ids)} plasmids that contained the winning features")
    logging.info(f"In all, {found_feature_count} out of the {len(cds_regions)} features hit at least one plasmid")

    write_feature_to_plasmid_json(feat_to_plasmid_dict, feature_to_plasmid_json_path)

    return all_matching_plasmid_ids

def get_gff_path(gff_dir: Path) -> Path:
    """
    Find the .gff file in the directory.

    Args:
        gff_dir: str, path to the directory containing the GFF file.

    Returns:
        Path, the path to the GFF file.

    Raises:
        ValueError: if there is not exactly one .gff file in the directory.
        FileNotFoundError: if the directory does not exist.
        Exception: if any other error occurs.
    """
    try:
        gff_files = list(Path(gff_dir).glob("*.gff"))
        if len(gff_files) != 1:
            raise ValueError(f"There should be exactly one .gff file, but found {len(gff_files)} in {gff_dir}")
        return gff_files[0]
    except ValueError as e:
        logging.error(f"Error finding GFF file: {e}") # include the exception message itself.
        sys.exit(1)
    except FileNotFoundError:
        logging.error(f"Directory {gff_dir} not found")
        sys.exit(1)
    except Exception as e:
        logging.error(f"An unexpected error occurred while finding GFF file: {e}")
        sys.exit(1)

def pick_random_genes(gff_dir: Path, num_to_pick: int) -> List[Tuple[str, int, int]]:
    """
    Pick a specified number of random CDS features from the GFF file in the given directory.

    Args:
        gff_dir: str, path to the directory containing the GFF file.
        num_to_pick: int, the number of genes to pick.

    Returns:
        list of tuples, where each tuple contains the feature ID, start position and end position of one of the randomly selected genes.

    Raises:
        ValueError: if the number of genes requested is greater than the number of CDS features in the GFF file.
        FileNotFoundError: if the directory does not exist.
        gffutils.exceptions.FeatureNotFoundError: if the GFF file does not contain a 'CDS' feature type.
        Exception: if any other error occurs.
    """

    try:
        gff_path = str(get_gff_path(gff_dir))
        db = gffutils.create_db(gff_path, ":memory:")

        count = db.count_features_of_type('CDS')
        if num_to_pick > count:
            raise ValueError(f"Requested {num_to_pick} genes, but only {count} available")

        cds_list = [(feat.id, feat.start, feat.end) for feat in db.all_features(featuretype='CDS')]
        winning_features = random.sample(cds_list, num_to_pick)
        return winning_features

    except ValueError as e:
        logging.error(f"Error picking random genes: {e}")
        sys.exit(1)
    except gffutils.exceptions.FeatureNotFoundError as e:
        logging.error(f"Error with GFF file structure: {e}")
        sys.exit(1)
    except Exception as e:
        logging.error(f"An unexpected error occurred while picking random genes: {e}")
        sys.exit(1)

def load_plasmids(plasmids_json_path: Path) -> List[Dict[str, Union[int, str]]]:
    with open(plasmids_json_path) as f:
        plasmids = json.load(f)
    return plasmids

def simulate_read_counts(
    df_samples: pd.DataFrame, plasmid_count: int, count_range: List[int],
    plusminus: int, winners: List[int], strength: int
) -> pd.Series:
    """
    For each samples, generates a list of read counts for each barcode by doing the following:

    1. Establish a primary distribution for all barcodes: these will be in range `count_range`
    2. All samples begin with slight variations, within `plusminus` distance, of the primary distribution
    3. Non-Time0 barcodes will be shifted by `winner_strength` for the "winning" features

    Args:
        df_samples: Dataframe with a 'Group' column and other columns that will be returned unchanged
        plasmid_count: Number of barcodes
        count_range: Range (inclusive) of values for the primary distribution
        plusminus: Maximum distance (in either direction) to shift the primary distribution
        winners: List of indices of the "winning" features
        strength: Amount to shift the "winning" features by

    Returns:
        Series with the same index as `df_samples` and a value for each barcode which is a list of read counts
    """

    def _vary(x, y):
        return max(x + random.randint(-plusminus, plusminus) + y, 0)

    def simple_variation(row):
        return [_vary(x, 0) for x in primary_counts]

    def winner_variation(row, winners):
        return [_vary(x, 0) if ix not in winners else _vary(x, strength) for ix, x in enumerate(primary_counts)]

    primary_counts = [random.randint(count_range[0], count_range[1]) for _ in range(plasmid_count)]

    time0rows = df_samples['Group'] == 'Time0'

    return pd.concat([
        df_samples[time0rows].apply(simple_variation, axis=1),
        df_samples[~time0rows].apply(winner_variation, winners=winners, axis=1)
    ], axis=0)

def export_to_fastq(output_file: Path, df_samples: pd.DataFrame, plasmids: List[Dict[str, Union[int, str]]], quality_score: int):
    output_file.parent.mkdir(exist_ok=True)
    with open(output_file, 'w') as f:
        read_count = 0
        for row in df_samples.itertuples():
            index=row.index2
            for ix, counts in enumerate(row.results):
                n = "".join([random.choice(['A', 'T', 'C', 'G']) for x in range(row.nN)])
                barcode = plasmids[ix]['barcode_sequence']
                read_sequence = f"{n}{index}GTCGACCTGCAGCGTACG{barcode}AGAGACCTCGTGGAC"
                quality_string = chr(quality_score + 33) * len(read_sequence)  # Phred+33 encoding
                for copies in range(counts):
                    read_count += 1
                    f.write(f"@read{read_count}\n{read_sequence}\n+\n{quality_string}\n")
    logging.info(f'Fastq written to {output_file}')

def main() -> None:
    args = parse_args()
    working_dir = Path(args.working_dir)
    json_path = working_dir / 'library' / 'plasmids.json'
    ref_dir = working_dir / 'ref'
    output_file_dir = working_dir / 'barseq'
    output_json_path = output_file_dir / 'feature_to_plasmid.json'
    output_reads_path = output_file_dir / 'reads.fastq'
    log_path = working_dir / 'log.txt'

    for directory in [working_dir, json_path.parent, ref_dir]:
        if not directory.is_dir():
            raise FileNotFoundError(f"The following does not exist or is not a directory: {directory}")

    if not log_path.exists():
        print(f"Unexpectedly, there is no log file in the expected location. Creating a new log file here: {log_path}")

    args = parse_args()
    setup_logging("Simulating Barseq Reads", file_path=log_path)
    log_args(vars(args))

    # Pick winning features (returns a list of (feature_name, start, end))
    winning_features = pick_random_genes(ref_dir, args.winner_count)

    # Get sample information
    samples_path = ref_dir / 'bobaseq_barseq_samples.tsv'
    multiplex_path = ref_dir / 'barseq4.index2'
    copyfile(args.samples_tsv_path, samples_path)
    copyfile(args.multiplex_index_tsv, multiplex_path)
    df_samples = get_samples_df(samples_path, multiplex_path)

    # Find (indices of) plasmids that contain those features
    plasmids = load_plasmids(json_path)
    matching_plasmid_ids = find_overlapping_plasmids(plasmids, winning_features, output_json_path)

    # Simulate read counts
    df_samples['results'] = simulate_read_counts(
        df_samples = df_samples,
        plasmid_count = len(plasmids),
        count_range = args.count_range,
        plusminus = args.plusminus,
        winners = matching_plasmid_ids,
        strength = args.winner_strength
    )

    # Export to fastq
    export_to_fastq(
        output_file = output_reads_path,
        df_samples = df_samples,
        plasmids = plasmids,
        quality_score = args.quality_score
    )

if __name__ == '__main__':
    main()
