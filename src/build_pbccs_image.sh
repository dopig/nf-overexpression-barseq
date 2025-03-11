#!/bin/bash
set -eu

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONTAINER_DIR="$SCRIPT_DIR/../containers/pbccs"

# Build the pbccs container
docker build -t pbccs-container $CONTAINER_DIR

echo "Container pbccs-container built successfully."