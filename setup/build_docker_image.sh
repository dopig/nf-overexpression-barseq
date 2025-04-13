#!/bin/bash
set -euo pipefail

#!/bin/bash

# Get the directory where the script is located
SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"

# Define the relative path to the docker directory from the script's location
DOCKER_RELATIVE_TO_SCRIPT="../docker"
DOCKER_ABSOLUTE_DIR="${SCRIPT_DIR}/${DOCKER_RELATIVE_TO_SCRIPT}"

# Define the setup directory (where the script lives)
SETUP_ABSOLUTE_DIR="${SCRIPT_DIR}"
SCRIPT_NAME="$(basename "$0")"
SCRIPT_PATH="${SETUP_ABSOLUTE_DIR}/${SCRIPT_NAME}"

# Function to check if Docker is running
check_docker_running() {
  if ! command -v docker &> /dev/null; then
    echo "Error: Docker command not found. Please ensure Docker is installed and in your PATH."
    return 1
  fi

  if ! docker info > /dev/null 2>&1; then
    echo "Error: Docker daemon is not running. Please start Docker."
    return 1
  fi
  return 0
}

# Function to check if a Docker image exists
check_image_exists() {
  local image_name="$1"
  local image_id=$(docker images -q "$image_name")
  if [ -n "$image_id" ]; then
    return 0 # Image exists
  else
    return 1 # Image does not exist
  fi
}

# Check if Docker is running before proceeding
if ! check_docker_running; then
  exit 1
fi

# Set a flag to control skipping existing images (default is to rebuild)
SKIP_EXISTING=true

# Process command-line arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --force)
      SKIP_EXISTING=false
      shift
      ;;
    *)
      echo "Unknown option: $1"
      echo "Usage (to force remaking all repo images): $0 [--force]"
      exit 1
      ;;
  esac
done

if ! "$SKIP_EXISTING"; then
  echo "Note: to force recreate existing images in the future, use: build_docker_image.sh --force"
fi

# Iterate over the directories in $DOCKER_ABSOLUTE_DIR, excluding 'shared'
for dir in "${DOCKER_ABSOLUTE_DIR}"/*; do
  if [[ -d "$dir" && "$(basename "$dir")" != "shared" ]]; then
    image_name="$(basename "$dir")"
    dockerfile_path="${dir}/Dockerfile"

    if [[ -f "$dockerfile_path" ]]; then
      if "$SKIP_EXISTING" && check_image_exists "$image_name"; then
        echo "Image '$image_name' already exists. Skipping build."
      else
        echo "Building Docker image '$image_name' from '../docker/$image_name'..."
        docker build -t "$image_name" "${dir}"
        if [ $? -eq 0 ]; then
          echo "-- Successfully built Docker image '$image_name'."
        else
          echo "-- Error building Docker image '$image_name'."
          exit 1
        fi
      fi
    else
      echo "** Warning: Dockerfile not found in '$dir'. Skipping. **"
    fi
  fi
done

echo "Finished processing Docker directories."

# Make the script executable (this might already be done, but it's good to ensure)
chmod +x "${SCRIPT_PATH}"
