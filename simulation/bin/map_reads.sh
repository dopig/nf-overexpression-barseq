#!/bin/bash
set -euo pipefail

# Ensure correct usage
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <session_data_dir>"
    exit 1
fi

work_dir="$1"

# If .DS_Store file exists, remove it so it doesn't interfere
if [ -f "$work_dir/map_temp/.DS_Store" ]; then
    rm "$work_dir/map_temp/.DS_Store"
fi

command="/work/ref/bobaseq_config.json /work/map_temp /work/map 1"

echo "Running: docker run --rm -v "$work_dir":/work bobaseq-image $command"

docker run --rm -v "$work_dir":/work bobaseq-image $command

# Remove map_temp dir
rm -rf "$work_dir/map_temp"

echo "Mapping complete"