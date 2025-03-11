#!/bin/bash
set -euo pipefail

# Set the environment path
envpath="./envs/sim-overex-lib"

# Set the Python command and arguments
# The current ones are for a tiny library for testing
pycommand_args=(
    "src/make_sim_library.py"
    "PATH/TO/FASTAFILE.fna"
    "--unique_barcodes" "5000"
    "--library_size" "500"
    "--coverage" "10"
)

# Run the Python command within the Micromamba environment
echo "Running Python script in micromamba environment: $envpath"
echo "Command: python3 ${pycommand_args[@]}"

micromamba run -p "$envpath" python3 "${pycommand_args[@]}"
