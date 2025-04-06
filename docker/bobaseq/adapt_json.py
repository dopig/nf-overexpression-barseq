#!/usr/bin/env python3

import argparse
import json

parser = argparse.ArgumentParser(description='Modify JSON file for Bobaseq script')
parser.add_argument('json', help='Input JSON file (REQUIRED)')
parser.add_argument('--name', required=True, help='Library name (REQUIRED)')
parser.add_argument('--fasta', required=True, help='Input genome FASTA filename (REQUIRED)')
parser.add_argument('--gff', required=True, help='Input genome GFF filename (REQUIRED)')
parser.add_argument('--oligos', required=True, help='Input oligo FASTA path (REQUIRED)')
args = parser.parse_args()

with open(args.json, 'r') as f:
    json_data = json.load(f)

json_data['lib_names'] = [args.name]
json_data['lib_genome_dir'] = '.'
json_data['lib_genome_filenames'] = [args.fasta]
json_data['lib_genome_gffs'] = [args.gff]
json_data['primer_info']['oligo_db_fp'] = args.oligos

with open('bobaseq_config.json', 'w') as f:
    json.dump(json_data, f, indent=4)
