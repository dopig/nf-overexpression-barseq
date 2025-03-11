#!/bin/bash
set -eu  # Exit on error

# Ensure correct usage
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <local_data_dir> <file_in> <file_out>"
    exit 1
fi

local_data_dir="$1"
file_in="$2"
file_out="$3"
log_file="pbccs-log.txt"

mapping=""$local_data_dir":/work"
command="ccs "$file_in" "$file_out""

docker run --rm -v $mapping pbccs-container $command

echo "Processing complete: $file_out"