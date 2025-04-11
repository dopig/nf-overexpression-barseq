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
- barseq/reads/feature_to_plasmid.json: Mapping of CDS features to plasmids.
- barseq/reads/reads.fastq: Simulated FASTQ file with BarSeq reads.
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
    python3 simulate_barseq_reads.py --working-dir /path/to/out_dir --winner-count 30
"""

import argparse
import json
import random
import logging
import sys
from pathlib import Path
from shutil import copyfile
from typing import List, Dict, Tuple, Union

import gffutils
import pandas as pd

from utils import setup_logging, format_output, log_args


SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR = SCRIPT_DIR.parent.parent
DATA_DIR = ROOT_DIR / 'data'
TIME0_NAME = 'Time0'


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('--out_dir', default=".", help='Output directory (default: .)')
    parser.add_argument('--gff-path', required=True, help='Path to GFF file')
    parser.add_argument('--plasmid-json-path', required=True, help='Path to plasmid JSON file')
    parser.add_argument('--samples-tsv-path', default=DATA_DIR / 'reference' / 'bobaseq_barseq_samples.tsv',
        help='Path to samples TSV file. This script needs the index, Group, and desc columns. If no file path is provided, data/reference/bobaseq_barseq_samples.tsv will be used')
    parser.add_argument('--multiplex-index-tsv', default=ROOT_DIR / 'shared' / 'external' / 'primers' / 'barseq4.index2',
        help='Path to multiplex primer barseq4 TSV file. If no file path is provided, data/reference/barseq4.index2 will be used')
    parser.add_argument('--winner-count', type=int, default=5, help='How many will be selected as winners (default: 20)')
    parser.add_argument('--winner-strength', type=int, default=5000, help='How many counts stronger winners will be (default: 5000)')
    parser.add_argument('--count-range', type=int, nargs=2, default=[0, 1000], help='Count range (min, max), entered as "--count-range min max" (default: 0 1000)')
    parser.add_argument('--plusminus', type=int, default=1, help='Plus/minus value (default: 1)')
    parser.add_argument('--quality-score', type=int, default=30, help='Quality score to be used for all bases (default: 30)')
    parser.add_argument('--set-random-seed', default=False, action='store_true', help='Set random seed to 42 (default: False)')
    return parser.parse_args()

def get_samples_df(samples_path: Path, index_path: Path) -> pd.DataFrame:
    """
    Read in samples TSV file and Barseq index TSV file.

    Returns a pandas dataframe with columns 'Group', 'index_name', 'index',
    and 'name'. The 'index' column is the actual index sequence.
    """
    df_samples = pd.read_csv(samples_path, sep='\t', usecols=['index', 'Group', 'desc', 'lib'])
    df_samples.rename(columns={'index': 'index_name'}, inplace=True)

    df_barseq4 = pd.read_csv(index_path, sep='\t')
    df_samples = pd.merge(df_samples, df_barseq4, on='index_name')

    return df_samples

def load_plasmids(plasmids_json_path: Path) -> List[Dict[str, Union[int, str]]]:
    with open(plasmids_json_path) as f:
        plasmids = json.load(f)
    return plasmids

def find_overlapping_plasmids(plasmids: List[Dict[str, Union[int, str]]], start: int, end: int) -> List[int]:
    """
    Given list of plasmids
    """
    feature_matching_plasmid_ids = []
    for plasmid in plasmids:
        # Cover normal case (first line below) and then if plasmid spans the origin
        if (plasmid['start'] <= start and plasmid['end'] >= end) or \
            (plasmid['start'] > plasmid['end'] and (end < plasmid['end'] or start > plasmid['start'])):
            feature_matching_plasmid_ids.append(plasmid['id'])
    return feature_matching_plasmid_ids

def pick_random_genes(gff_path: Path, num_to_pick: int) -> List[gffutils.feature.Feature]:
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
        db = gffutils.create_db(gff_path, ":memory:")

        count = db.count_features_of_type('CDS')
        if (num_to_pick * 2) > count:
            raise ValueError(f"Requested {num_to_pick} genes, but only {count} genes available, and you are only permitted to pick up to half of the genome as winners.")

        crude_winning_features = random.sample(list(db.all_features(featuretype='CDS')), num_to_pick * 2)
        crude_winning_proteins = [feat for feat in crude_winning_features if 'pseudo' not in feat.attributes][:num_to_pick]
        return crude_winning_proteins

    except ValueError as e:
        logging.error(f"Error picking random genes: {e}")
        sys.exit(1)
    except gffutils.exceptions.FeatureNotFoundError as e:
        logging.error(f"Error with GFF file structure: {e}")
        sys.exit(1)
    except Exception as e:
        logging.error(f"An unexpected error occurred while picking random genes: {e}")
        sys.exit(1)

def make_feature_plasmid_dict(group: 'str',
    plasmids: List[Dict[str, Union[int, str]]],
    winning_features: List[gffutils.feature.Feature]) -> List[ Dict ]:

    winners = [
        {
            'group': group,
            'id': feat.id,
            'locus_tag': feat.attributes.get('locus_tag',[None])[0],
            'protein_id': feat.attributes.get('protein_id',[None])[0],
            'plasmid_ids': find_overlapping_plasmids(plasmids, feat.start, feat.end),
        } for feat in winning_features
    ]
    return winners

def check_plasmid_success(feature_to_plasmid_ids, group, len_plasmids, len_winning_features):
    found_feature_count = sum(1 for locus in feature_to_plasmid_ids if locus['plasmid_ids'])
    if found_feature_count == 0:
        logging.error(f"For group '{group}', none of the {len_plasmids} plasmids contained any of the {len_winning_features} CDSs")
        sys.exit(1)
    logging.info(f"For group '{group}', {found_feature_count} out of the {len_winning_features} features hit at least one plasmid")

def get_winning_plasmid_ids(
    groups: List[str],
    plasmids: List[Dict[str, Union[int, str]]],
    gff_path: Path,
    winner_count: int,
    winners_tsv_path: Path
) -> Dict[str, List[int]]:

    group_feature_plasmid = []

    logging.info(f"Processing {len(groups)} groups...")
    for group in groups:
        # Pick winning features for this group
        winning_features = pick_random_genes(gff_path, winner_count)

        feature_plasmid_dicts = make_feature_plasmid_dict(group, plasmids, winning_features)

        check_plasmid_success(feature_plasmid_dicts, group, len(plasmids), len(winning_features))

        group_feature_plasmid.extend(feature_plasmid_dicts)

    df_gfp = pd.DataFrame(group_feature_plasmid)

    # Save the detailed table for debugging
    df_gfp.to_csv(winners_tsv_path, sep='\t', index=False)
    logging.debug(f'Writing chosen winner details to {winners_tsv_path}')

    group_to_plasmid_id = df_gfp.groupby('group')['plasmid_ids'].apply(lambda x: set().union(*x)).to_dict()

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
                n = "".join([random.choice(['A', 'T', 'C', 'G']) for x in range(row.nN)])
                barcode = plasmid['barcode_sequence']
                read_sequence = f"{n}{index}GTCGACCTGCAGCGTACG{barcode}AGAGACCTCGTGGAC"
                quality_string = chr(quality_score + 33) * len(read_sequence)  # Phred+33 encoding
                counts = row.list_of_read_counts[ix]
                for copies in range(counts):
                    read_count += 1
                    f.write(f"@read{read_count}\n{read_sequence}\n+\n{quality_string}\n")
    logging.info(f'Fastq written to {output_file}')

# def make_json_lib(df_samples: pd.DataFrame, out_dir: Path) -> None:
#     """
#     Fitness scripts expect this JSON to be built
#     """
#     lib_json_path = out_dir / 'barseq' / 'reads' / 'lib.json'
#     lib_value = str(df_samples['lib'].unique()[0])
#     dir_value = (out_dir / "map" / lib_value / "05-BC_and_genes_dfs").resolve()

#     # Get the ft_path
#     feature_tables = list((out_dir / "ref").glob("*_feature_table.txt*"))
#     ft_path = feature_tables[0].resolve()

#     lib_dict = {
#         "lib": lib_value,
#         "dir": str(dir_value),
#         "feature_table_path": str(ft_path)
#     }

#     with open(lib_json_path, 'w') as f:
#         json.dump(lib_dict, f, indent=4)

#     df_lib = pd.DataFrame([lib_dict])
#     df_lib.to_csv(lib_json_path.with_suffix('.csv'), index=False)


def main() -> None:
    args = parse_args()
    out_dir = Path(args.out_dir)

    winners_tsv_path = out_dir / 'chosen_winners.tsv'
    output_reads_path = out_dir / 'reads.fastq'
    log_path = out_dir / 'log.txt'

    if args.set_random_seed:
        random.seed(42)

    if not log_path.exists():
        print(f"Unexpectedly, there is no log file in the expected location. Creating a new log file here: {log_path}")

    out_dir.mkdir(parents=True, exist_ok=True)

    setup_logging("Simulating Barseq Reads", file_path=log_path)
    raw_command_line = " ".join(sys.argv)
    logging.debug(f"Command run: {raw_command_line}")
    log_args(vars(args))

    df_samples = get_samples_df(args.samples_tsv_path, args.multiplex_index_tsv)
    plasmids = load_plasmids(args.plasmid_json_path)

    unique_desc = df_samples[df_samples['Group'] != TIME0_NAME].desc.unique()

    group_to_plasmid_ids = get_winning_plasmid_ids(unique_desc, plasmids, args.gff_path, args.winner_count, winners_tsv_path)

    # Simulate read counts
    df_samples['list_of_read_counts'] = simulate_read_counts(
        df_samples = df_samples,
        plasmid_count = len(plasmids),
        count_range = args.count_range,
        plusminus = args.plusminus,
        winners_dict = group_to_plasmid_ids,
        strength = args.winner_strength
    )
    # df_samples.to_csv(out_dir / 'barseq_samples.tsv', sep='\t', index=False)

    export_to_fastq(
        output_file = output_reads_path,
        df_samples = df_samples,
        plasmids = plasmids,
        quality_score = args.quality_score
    )

    # make_json_lib(df_samples, out_dir)

if __name__ == '__main__':
    main()
