#!/bin/bash
set -eu

# Ensure correct usage
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <image-to-build (e.g. pbccs, bobaseq)>"
    exit 1
fi

image_name="$1"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DOCKERFILE_DIR="$SCRIPT_DIR/../docker/$image_name"

# Ensure the specified Dockerfile directory exists
if [ ! -d "$DOCKERFILE_DIR" ]; then
    echo "Error: Directory '$DOCKERFILE_DIR' not found."
    exit 1
fi

# Build the image
docker build -t $image_name-image $DOCKERFILE_DIR

echo "Successfully built $image_name-image docker image."
