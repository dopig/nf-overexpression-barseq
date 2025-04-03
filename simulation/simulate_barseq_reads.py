#!/usr/bin/env python3

"""
Simulates BarSeq reads by selecting "winner" CDS features, mapping them to plasmids,
and generating synthetic FASTQ data with varying read counts. When determining "groups"
of samples, the "desc" column is used as a groupby column so replicates look similar.

Overview:
- Selects random "winner" CDS features from a GFF file.
- Determines which plasmids contain those features.
- Assigns read counts to plasmid barcodes with variations.
- Exports simulated BarSeq reads in FASTQ format.

Required Input Files (inside --working-dir):
- library/plasmids.json: Plasmid metadata.
- ref/*.gff: GFF file containing CDS feature data.
- reference/bobaseq_barseq_samples.tsv: Sample index/group data.
- external/primers/barseq4.index2: Multiplex primer index information.

Generated Output:
- barseq/feature_to_plasmid.json: Mapping of CDS features to plasmids.
- barseq/reads.fastq: Simulated FASTQ file with BarSeq reads.
- log.txt: Run log.

Command-line Parameters:
- --working-dir: Root directory containing input files.
- --samples-tsv-path: Path to sample metadata (default: reference/bobaseq_barseq_samples.tsv).
- --multiplex-index-tsv: Path to primer index file (default: reference/barseq4.index2).
- --winner-count: Number of CDS features to mark as "winners" (default: 20).
- --winner-strength: Read count increase for "winner" features (default: 5000).
- --count-range: Min/max range of read counts per barcode (default: 0 1000).
- --plusminus: Variation in read counts (default: 1).
- --quality-score: Quality score for all bases in FASTQ output (default: 30).

Example:
    python3 simulate_barseq_reads.py --working-dir /path/to/working_dir --winner-count 30
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

from bin.utils import setup_logging, format_output, log_args


SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR = SCRIPT_DIR.parent
DATA_DIR = ROOT_DIR / 'data'
TIME0_NAME = 'Time0'


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('working_dir', help='Working directory which has library, map, and ref subdirectories')
    parser.add_argument('--samples-tsv-path', default=DATA_DIR / 'reference' / 'bobaseq_barseq_samples.tsv',
        help='Path to samples TSV file. This script needs the index, Group, and desc columns. If no file path is provided, data/reference/bobaseq_barseq_samples.tsv will be used')
    parser.add_argument('--multiplex-index-tsv', default=ROOT_DIR / 'shared' / 'external' / 'primers' / 'barseq4.index2',
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
    df_samples = pd.read_csv(samples_path, sep='\t', usecols=['index', 'Group', 'desc'])
    df_samples.rename(columns={'index': 'index_name'}, inplace=True)

    df_barseq4 = pd.read_csv(index_path, sep='\t')
    df_samples = pd.merge(df_samples, df_barseq4, on='index_name')

    return df_samples

def load_plasmids(plasmids_json_path: Path) -> List[Dict[str, Union[int, str]]]:
    with open(plasmids_json_path) as f:
        plasmids = json.load(f)
    return plasmids

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

def find_overlapping_plasmids(
    plasmids: List[Dict[str, Union[int, str]]],
    cds_regions: List[Tuple[str, int, int]]
) -> Dict[str, List[str]]:
    """
    Given a list of plasmids and a list of CDS regions, return a dict mapping
    each CDS feature ID to a set of all the plasmid IDs that contain it
    """
    feat_to_plasmid_dict = dict()
    for feature_id, cds_start, cds_end in cds_regions:
        feature_matching_plasmid_ids = []
        for plasmid in plasmids:
            # Cover normal case (first line below) and then if plasmid spans the origin
            if (plasmid['start'] <= cds_start and plasmid['end'] >= cds_end) or \
               (plasmid['start'] > plasmid['end'] and (cds_end < plasmid['end'] or cds_start > plasmid['start'] )):
                feature_matching_plasmid_ids.append(plasmid['id'])
                # print('--',feature_id, cds_start, cds_end)
                # print('-----', plasmid)
        feat_to_plasmid_dict[feature_id] = feature_matching_plasmid_ids
    return feat_to_plasmid_dict

    feat_to_plasmid_dict = defaultdict(list)  # No need to check if key exists
    for feature_id, cds_start, cds_end in cds_regions:
        for plasmid in plasmids:
            if (plasmid['start'] <= cds_start and plasmid['end'] >= cds_end) or \
            (plasmid['start'] > plasmid['end'] and (cds_end < plasmid['end'] or cds_start > plasmid['start'])):
                feat_to_plasmid_dict[feature_id].append(plasmid['id'])

def write_feature_to_plasmid_json(feature_to_plasmid_dict: Dict[str, List[str]], output_file_path: Path) -> None:
    """
    Write a JSON file mapping feature IDs to lists of overlapping plasmid IDs.

    :param feature_to_plasmid_dict: A dict mapping feature IDs to lists of overlapping plasmid IDs
    :param output_file_path: The path to which the JSON file should be written
    """
    output_file_path.parent.mkdir(exist_ok=True)
    with open(output_file_path, 'w') as f:
        json.dump(feature_to_plasmid_dict, f, indent=4)
    logging.debug(f'Writing feature-to-plasmid JSON for {len(feature_to_plasmid_dict)} groups to {output_file_path}')

def get_all_plasmid_ids(dictionary):
    """
    Return a set of all the plasmid IDs found in the dictionary.

    :param dictionary: A dictionary where some values are lists of plasmid IDs
    :return: A set of all the plasmid IDs found in the dictionary
    """
    return {num for value in dictionary.values() if isinstance(value, list) for num in value}

def get_winning_plasmid_ids(
    groups: List[str],
    plasmids: List[Dict[str, Union[int, str]]],
    ref_dir: Path,
    winner_count: int,
    output_json_path: Path
) -> Dict[str, List[int]]:
    group_to_feature_to_plasmid_ids = {}
    group_to_plasmid_id = {}
# Loop through each group (except 'Time0')
    logging.info(f"Processing {len(groups)} groups...")
    for group in groups:
        # Pick winning features for this group
        winning_features = pick_random_genes(ref_dir, winner_count)
        # Find matching plasmid ids for this group
        feature_to_plasmid_ids = find_overlapping_plasmids(plasmids, winning_features)
        found_feature_count = sum(1 for value in feature_to_plasmid_ids.values() if value != [])
        if found_feature_count == 0:
            logging.error(f"For group '{group}', none of the {len(plasmids)} plasmids contained any of the {len(winning_features)} CDSs")
            sys.exit(1)
        logging.info(f"For group '{group}', {found_feature_count} out of the {len(winning_features)} features hit at least one plasmid")
        # Store the matching plasmid ids in the dictionary
        group_to_feature_to_plasmid_ids[group] = feature_to_plasmid_ids
        group_to_plasmid_id[group] = get_all_plasmid_ids(feature_to_plasmid_ids)

    write_feature_to_plasmid_json(group_to_feature_to_plasmid_ids, output_json_path)

    return group_to_plasmid_id

def simulate_read_counts(
    df_samples: pd.DataFrame, plasmid_count: int, count_range: List[int],
    plusminus: int, winners_dict: Dict[str, List[int]], strength: int
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
        winners_dict: Dict with keys of group names and values of lists of indices of the "winning" plasmids
        strength: Amount to shift the "winning" features by

    Returns:
        Series with the same index as `df_samples` and a value for each barcode which is a list of read counts
    """

    def _vary(x, y):
        return max(x + random.randint(-plusminus, plusminus) + y, 0)

    def simple_variation(row):
        return [_vary(x, 0) for x in primary_counts]

    def winner_variation(row, winners_dict):
        # Plasmids are not 0-indexed, so need to add 1
        return [_vary(x, 0) if ix+1 not in winners_dict[row.desc] else _vary(x, strength) for ix, x in enumerate(primary_counts)]

    primary_counts = [random.randint(count_range[0], count_range[1]) for _ in range(plasmid_count)]

    logging.debug(f"Simulated primary read counts: {primary_counts[:10]}... (showing first 10)")

    time0rows = df_samples['Group'] == TIME0_NAME

    return pd.concat([
        df_samples[time0rows].apply(simple_variation, axis=1),
        df_samples[~time0rows].apply(winner_variation, winners_dict=winners_dict, axis=1)
    ], axis=0)

def export_to_fastq(output_file: Path, df_samples: pd.DataFrame, plasmids: List[Dict[str, Union[int, str]]], quality_score: int):
    output_file.parent.mkdir(exist_ok=True)
    with open(output_file, 'w') as f:
        read_count = 0
        for row in df_samples.itertuples():
            index=row.index2
            for ix, plasmid in enumerate(plasmids):
                # print(ix, '--', row.list_of_read_counts[ix], '--',plasmid)
                # quit()
                n = "".join([random.choice(['A', 'T', 'C', 'G']) for x in range(row.nN)])
                barcode = plasmid['barcode_sequence']
                read_sequence = f"{n}{index}GTCGACCTGCAGCGTACG{barcode}AGAGACCTCGTGGAC"
                quality_string = chr(quality_score + 33) * len(read_sequence)  # Phred+33 encoding
                counts = row.list_of_read_counts[ix]
            # for ix, counts in enumerate(row.list_of_read_counts):
            #     n = "".join([random.choice(['A', 'T', 'C', 'G']) for x in range(row.nN)])
            #     barcode = plasmids[ix-1]['barcode_sequence'] # Plasmid IDs start at 1
            #     read_sequence = f"{n}{index}GTCGACCTGCAGCGTACG{barcode}AGAGACCTCGTGGAC"
            #     quality_string = chr(quality_score + 33) * len(read_sequence)  # Phred+33 encoding
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

    setup_logging("Simulating Barseq Reads", file_path=log_path)
    log_args(vars(args))

    # Get sample information
    working_samples_path = ref_dir / 'bobaseq_barseq_samples.tsv'
    working_multiplex_path = ref_dir / 'barseq4.index2'
    copyfile(args.samples_tsv_path, working_samples_path)
    copyfile(args.multiplex_index_tsv, working_multiplex_path)

    df_samples = get_samples_df(working_samples_path, working_multiplex_path)
    plasmids = load_plasmids(json_path)

    unique_desc = df_samples[df_samples['Group'] != TIME0_NAME].desc.unique()

    group_to_plasmid_ids = get_winning_plasmid_ids(unique_desc, plasmids, ref_dir, args.winner_count, output_json_path)
    print(group_to_plasmid_ids)
    # Simulate read counts
    df_samples['list_of_read_counts'] = simulate_read_counts(
        df_samples = df_samples,
        plasmid_count = len(plasmids),
        count_range = args.count_range,
        plusminus = args.plusminus,
        winners_dict = group_to_plasmid_ids,
        strength = args.winner_strength
    )
    df_samples.to_csv(output_file_dir / 'barseq_samples.tsv', sep='\t', index=False)
    # Export to fastq
    export_to_fastq(
        output_file = output_reads_path,
        df_samples = df_samples,
        plasmids = plasmids,
        quality_score = args.quality_score
    )

if __name__ == '__main__':
    main()
