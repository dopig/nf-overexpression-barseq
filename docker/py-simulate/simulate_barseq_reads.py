#!/usr/bin/env python3

"""
simulate_barseq_reads.py

Part of the nf-overexpression-barseq Nextflow pipeline.

This script simulates next-generation sequencing reads in FASTQ format
from a synthetically generated plasmid library. It's designed to mimic
real-world BarSeq experiments by incorporating variations in read counts.

Key functionalities:
- Selects random "winner" genes from a supplied genome file (GFF format).
- Identifies which plasmids (from a JSON file) contain these winner genes.
- Assigns variable read counts to plasmid barcodes. Plasmids with winner
  genes receive increased simulated read counts to reflect experimental
  enrichment.
- Exports the simulated reads in standard FASTQ format, including sequence
  and quality scores.

Note on sample grouping:
"Winner" genes are chosen independently for each sample group. The "desc"
column in the input samples TSV file is used to define these groups.
This column name is used for compatibility with other tools in our workflow.

Note on terminology:
"Gene" and "CDS" and "feature" are used interchangeably in the script.

Required inputs:
- Plasmid JSON file (path specified by --plasmid-json-path)
- GFF genome file (path specified by --gff-path)
- Samples TSV file (path specified by --samples-tsv-path)
- Multiplex index TSV file (path specified by --multiplex-index-tsv-path)

For a comprehensive list of all command-line arguments and their descriptions
(including optional parameters like read length, error rates, etc.),
run the script with the `--help` flag:
    python simulate_barseq_reads.py --help
"""

import argparse
import json
import random
import logging
import sys
import gzip
from pathlib import Path
from typing import List, Dict

import gffutils # type: ignore[import-untyped]
import pandas as pd

from utils import begin_log, Plasmid

TIME0_NAME = "Time0"

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out_dir", type=Path, default=".", help="Output directory (default: .)")
    parser.add_argument("--gff-path", type=Path, required=True, help="Path to GFF file")
    parser.add_argument("--plasmid-json-path", type=Path, required=True, help="Path to plasmid JSON file")
    parser.add_argument("--samples-tsv-path", type=Path, required=True, help='Path to samples TSV file. The "desc" column is used for separating into groups.')
    parser.add_argument("--multiplex-index-tsv-path", type=Path, required=True, help="Path to multiplex primer barseq4 TSV file.")
    parser.add_argument("--winner-count", type=int, default=5, help="How many will be selected as winners (default: 5)")
    parser.add_argument("--winner-strength", type=int, default=5000, help="How many counts stronger winners will be (default: 5000)")
    parser.add_argument("--count-range", type=int, nargs=2, default=[0, 1000], help='Count range (min, max), entered as "--count-range min max" (default: 0 1000)')
    parser.add_argument("--plusminus", type=int, default=1, help="Plus/minus value (default: 1)")
    parser.add_argument("--quality-score", type=int, default=30, help="Quality score to be used for all bases (default: 30)")
    parser.add_argument("--random-seed", type=int, default=None, help="Number to set random seed to (default: None)")
    return parser.parse_args()

def get_samples_df(samples_path: Path, index_path: Path) -> pd.DataFrame:
    """
    Read in samples TSV file and Barseq index TSV file.

    Returns:
        pd.DataFrame: DataFrame with columns "Group", "index_name", "index",
        and "name". The "index" column is the actual index sequence.
    """
    df_samples = pd.read_csv(
        samples_path,
        sep="\t",
        usecols=["index", "Group", "desc", "lib"],
        dtype={"index": str, "Group": str, "desc": str, "lib": str}
    ).rename(columns={"index": "multiplex_index_name"})

    df_barseq4 = pd.read_csv(
        index_path,
        sep="\t",
        usecols=["index_name", "index2", "nN"],
        dtype={"index_name": str, "index2": str, "nN": int}
    ).rename(columns={"index_name": "multiplex_index_name"})

    df_samples = pd.merge(df_samples, df_barseq4, on="multiplex_index_name")

    return df_samples

def load_plasmids(plasmids_json_path: Path) -> List[Plasmid]:
    with open(plasmids_json_path) as f:
        plasmids = json.load(f)
    return plasmids

def find_overlapping_plasmids(plasmids: List[Plasmid], start: int, end: int) -> List[int]:
    """
    Given list of plasmids, and a gene's starting and ending coordinates,
    return the IDs of the plasmids that contain the gene.

    There is some complexity here because genomes are circular and plasmids
    can span the origin (i.e. index 1 of the genome). For example, for a
    4,000,000 nucleotide genome. a 1,000 nucleotide plasmid could start at
    3,999,501 and end at 500. This needs to be accounted for.
    """
    feature_matching_plasmid_ids = []
    for plasmid in plasmids:
        # Cover typical case (first line below); in subsequent line, worry about rare cases when plasmid genomic origin
        if (plasmid["start"] <= start and plasmid["end"] >= end) or \
            (plasmid["start"] > plasmid["end"] and (end < plasmid["end"] or start > plasmid["start"])):
            feature_matching_plasmid_ids.append(plasmid["id"])
    return feature_matching_plasmid_ids

def pick_random_genes(gff_path: Path, num_to_pick: int) -> List[gffutils.feature.Feature]:
    """
    Picks a specified number of random CDS features from the supplied GFF file.
    Excludes pseudogenes.

    Args:
        gff_path: path to the GFF file.
        num_to_pick: int, the number of genes to pick.

    Returns:
        list of winning features
    """

    try:
        # Convert path to string for gffutils
        db = gffutils.create_db(str(gff_path), ":memory:")

        count = db.count_features_of_type("CDS")
        if (num_to_pick * 2) > count:
            raise ValueError(f"Requested {num_to_pick} genes, but only {count} genes available, and you are only permitted to pick up to half of the genome as winners.")

        # Pick double the amount to account for attrition in case of pseudogenes
        crude_winning_features = random.sample(list(db.all_features(featuretype="CDS")), num_to_pick * 2)
        crude_winning_proteins = [feat for feat in crude_winning_features if "pseudo" not in feat.attributes][:num_to_pick]
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

def make_feature_plasmid_dict(group: str,
    plasmids: List[Plasmid],
    winning_features: List[gffutils.feature.Feature]) -> List[Dict]:
    """
    Loops through supplied list of Features. For each feature, pulls out the
    relevants details AND runs find_overlapping_plasmids to identify
    the plasmids that contain that feature.

    Returns a list of dictionaries with the following keys:
        group: str - This is the value from the "desc" column in samples TSV
            (e.g. "LB-low salt")
        id: str - NCBI's ID for the feature (e.g. "cds-NP_050645.1")
        locus_tag: str - ID for the locus from GFF file (e.g. "Mup41")
        protein_id: str - NCBI's ID for the protein it encodes
            (e.g. "NP_050645.1")
        plasmid_ids: list - List of IDs of plasmids containing the feature
    """
    winners = [
        {
            "group": group,
            "id": feat.id,
            "locus_tag": feat.attributes.get("locus_tag",[None])[0],
            "protein_id": feat.attributes.get("protein_id",[None])[0],
            "plasmid_ids": find_overlapping_plasmids(plasmids, feat.start, feat.end),
        } for feat in winning_features
    ]
    return winners

def check_plasmid_success(
    feature_to_plasmid_ids: List[Dict],
    group: str,
    len_plasmids: int,
    len_winning_features: int
) -> None:
    """
    Within a group, ensures that at least one plasmid contains each of the
    winning features. Quits program if this is not the case.

    Args:
        feature_to_plasmid_ids: list of dictionaries, where each dictionary
            represents one of the winning features chosen for that group.
            The plasmids (if applicable) that contain the winning feature are
            listed under the 'plasmid_ids' key.
        len_plasmids: int - Number of plasmids in library.

    First it dermines found_feature_count, by checking how many of the winning
    features have at least one matching plasmid (and this the key 'plasmid_ids').
    The value of found_feature_count will range between 0 (if no feautres have
    any matching plasmids) and len_winning_features (if all features have
    matching plasmids).
    """
    found_feature_count = sum(1 for locus in feature_to_plasmid_ids if locus["plasmid_ids"])

    if found_feature_count == 0:
        logging.error(f'For group "{group}", none of the {len_plasmids} plasmids contained any of the {len_winning_features} CDSs')
        sys.exit(1)
    logging.info(f'For group "{group}", {found_feature_count} out of the {len_winning_features} features hit at least one plasmid')

def get_winning_plasmid_ids(
    groups: List[str],
    plasmids: List[Plasmid],
    gff_path: Path,
    winner_count: int,
    winners_tsv_path: Path
) -> Dict[str, List[int]]:
    """
    Given a list of groups, a list of plasmids, a GFF file path, the number of
    winning features to pick, and a path to write the detailed results to, picks
    the winning features for each group and determines the plasmids that contain
    each of them.

    Returns a dictionary with the group names as the keys, and a list of plasmid
    IDs as the values.
    """
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
    df_gfp.to_csv(winners_tsv_path, sep="\t", index=False)
    logging.debug(f"Writing chosen winner details to {winners_tsv_path}")

    group_to_plasmid_id = df_gfp.groupby("group")["plasmid_ids"].apply(lambda x: set().union(*x)).to_dict()

    return group_to_plasmid_id

def simulate_read_counts(
    df_samples: pd.DataFrame,
    plasmid_count: int,
    count_range: List[int],
    plusminus: int,
    winners_dict: Dict[str, List[int]],
    strength: int
) -> pd.Series:
    """
    For each sample, generates a list of read counts for each barcode by doing the following:

    1. Establish a primary distribution for all barcodes in range `count_range`
    2. So that samples in the same group are not identical, for each sample,
        vary by `plusminus` distance from the primary distribution
    3. Non-Time0 barcodes will be shifted by `winner_strength` for the "winning"
        features

    Args:
        df_samples: Dataframe with a "Group" column (among other columns)
        plasmid_count: Number of barcoded plasmids
        count_range: Range (inclusive) of values for the primary distribution
        plusminus: Maximum distance (in either direction) to shift the primary
            distribution
        winners_dict: Dict with keys of group names and values of lists of
            indices of the "winning" plasmids
        strength: Amount to shift the "winning" features by

    Returns:
        Series with the same index as `df_samples` and a value of list of read
            counts for each plasmid's barcode
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

    time0rows = df_samples["Group"] == TIME0_NAME

    return pd.concat([
        df_samples[time0rows].apply(simple_variation, axis=1),
        df_samples[~time0rows].apply(winner_variation, winners_dict=winners_dict, axis=1)
    ], axis=0)

def export_to_fastq_gz(
    output_file: Path,
    df_samples: pd.DataFrame,
    plasmids: List[Plasmid],
    quality_score: int
) -> None:
    """
    Exports sequencing reads to a gzipped FASTQ file.

    Constructs sequencing reads using barcode sequences from plasmids and sample
    information from a DataFrame. Each read is written to the specified gzipped
    FASTQ file with the provided quality score.

    Args:
        output_file: Path to the gzipped FASTQ file to be created.
        df_samples: DataFrame containing sample data, including read counts and
            indices.
        plasmids: List of dictionaries, each containing plasmid information such
            as barcode sequences.
        quality_score: Integer representing the base quality score for the reads.
    """
    output_file.parent.mkdir(exist_ok=True)
    read_count = 0
    with gzip.open(output_file, "wb") as gz_file:
        for _, row in df_samples.iterrows():
            index_val = row["index2"]
            for ix, plasmid in enumerate(plasmids):
                n = "".join([random.choice(["A", "T", "C", "G"]) for x in range(row["nN"])])
                barcode = plasmid["barcode_sequence"]
                read_sequence = f"{n}{index_val}GTCGACCTGCAGCGTACG{barcode}AGAGACCTCGTGGAC"
                # ASCII offset for Phred quality scores
                quality_string = chr(quality_score + 33) * len(read_sequence)
                counts = row["list_of_read_counts"][ix]
                for _ in range(counts):
                    read_count += 1
                    gz_file.write(f"@read{read_count}\n{read_sequence}\n+\n{quality_string}\n".encode("utf-8"))
    logging.info(f"Gzipped fastq written to {output_file}")

def main() -> None:
    args = parse_args()
    out_dir = Path(args.out_dir)

    winners_tsv_path = out_dir / "chosen_winners.tsv"
    output_reads_path = out_dir / "reads.fastq.gz"
    log_path = out_dir / "log.txt"

    if args.random_seed is not None:
        random.seed(args.random_seed)

    out_dir.mkdir(parents=True, exist_ok=True)

    begin_log("Simulating Barseq Reads", log_path, vars(args))

    df_samples = get_samples_df(args.samples_tsv_path, args.multiplex_index_tsv_path)
    plasmids = load_plasmids(args.plasmid_json_path)

    # Get names of groups of biological replicates
    unique_desc: List[str] = list(df_samples[df_samples["Group"] != TIME0_NAME].desc.unique())

    # Pick winning genes for each group, and identify the plasmids that contain them
    group_to_plasmid_ids = get_winning_plasmid_ids(unique_desc, plasmids, args.gff_path, args.winner_count, winners_tsv_path)

    # Simulate read counts
    df_samples["list_of_read_counts"] = simulate_read_counts(
        df_samples = df_samples,
        plasmid_count = len(plasmids),
        count_range = args.count_range,
        plusminus = args.plusminus,
        winners_dict = group_to_plasmid_ids,
        strength = args.winner_strength
    )

    export_to_fastq_gz(
        output_file = output_reads_path,
        df_samples = df_samples,
        plasmids = plasmids,
        quality_score = args.quality_score
    )

if __name__ == "__main__":
    main()
