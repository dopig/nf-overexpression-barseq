#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 1 ]; then
  echo "Error: Number of arguments: $#. Please provide one assembly ID as the only argument."
  exit 1
fi

# Access the argument
assembly_id="$1"

docsum=$(esearch -db assembly -query "$assembly_id" | efetch -format docsum)

species=$(echo $docsum | xtract -pattern DocumentSummary -element SpeciesName)
echo "$species" > species.txt

ftp_url=$(echo $docsum | xtract -pattern DocumentSummary -element FtpPath_RefSeq)

if [ -z "$ftp_url" ]; then
    echo "Could not retrieve FTP URL for $species assembly $assembly_id" >&2
    exit 1
else
    echo "FTP URL identified for $species: $ftp_url"
fi

scaffold_count=$(echo $docsum | xtract -pattern Stat -if @category -equals scaffold_count -and @sequence_tag -equals all -element Stat)
if [ $scaffold_count != 1 ]; then
    echo "Script currently limited to organisms with 1 scaffold; unfortunately $assembly_id has $scaffold_count" >&2
    exit 1
fi

# Extract the full file basename, e.g. GCF_000027325.1_ASM2732v1
base_name=$(basename "$ftp_url")

# Download the files
for ext in genomic.fna.gz genomic.gff.gz feature_table.txt.gz; do
    wget --no-verbose "${ftp_url}/${base_name}_${ext}"
done

echo "Downloading complete"
