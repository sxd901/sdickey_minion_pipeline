#!/bin/bash

# Set log file name (e.g., log to "install_images.log" in the current directory)
LOG_FILE="install_singularity_images.log"

# Function to log messages with timestamps
log_message() {
  local message="$1"
  echo "$(date '+%Y-%m-%d %H:%M:%S') - $message" | tee -a "$LOG_FILE"
}

# Log script start
log_message "Script started."

# Default TMPDIR location, can be overridden by environment variable
TMPDIR=${TMPDIR:-~/tmp}

# Check if 'images' directory exists, create it if not
if [ ! -d "images" ]; then
  log_message "Creating 'images' directory..."
  mkdir sdickey_minion_pipeline/images/
else
  log_message "'images' directory already exists."
fi

# Ensure TMPDIR exists
mkdir -p $TMPDIR
log_message "Temporary directory set to: $TMPDIR"

# Set the necessary environment variables for Singularity
export TMP=$TMPDIR
export SINGULARITY_TMPDIR=$TMPDIR

# Define an array of image names and their corresponding Docker image URIs
declare -A images
images=( 
  ["nanofilt"]="jdelling7igfl/nanofilt:2.8.0"
  ["seqkit"]="nanozoo/seqkit:2.6.1--022e008"
  ["cutadapt"]="olistr12/cutadapt:4.9"
  ["porechop"]="biocontainers/porechop:v0.2.4dfsg-1-deb_cv1"
  ["vsearch"]="olistr12/vsearch:2.28.1"
  ["r_env"]="olistr12/r_env:0.0.7"
  ["bbmap"]="nanozoo/bbmap:38.86--9ebcbfa"
  ["blast"]="ncbi/blast:2.15.0"
  ["ncbitax2lin"]="pegi3s/ncbitax2lin:2.3.1"  # Add ncbitax2lin image entry
)

# Check if Singularity is installed and the version is correct
if ! command -v singularity &> /dev/null; then
  log_message "Singularity could not be found. Please install Singularity before running this script."
  exit 1
fi

SINGULARITY_VERSION=$(singularity --version | awk '{print $3}')
EXPECTED_VERSION="3.1.0"
if [[ "$SINGULARITY_VERSION" != "$EXPECTED_VERSION" ]]; then
  log_message "Warning: Singularity version $SINGULARITY_VERSION detected. This script was written for version $EXPECTED_VERSION."
  log_message "Proceeding with the current version, but consider updating Singularity."
fi

# Iterate over the images array and build each image
for image_name in "${!images[@]}"; do
  docker_image="${images[$image_name]}"
  log_message "Building Singularity image for $image_name from $docker_image..."
  
  # Build the image and log the output
  singularity build "images/$image_name.sif" "docker://$docker_image" >> "$LOG_FILE" 2>&1
  
  # Check if the build was successful
  if [ $? -eq 0 ]; then
    log_message "$image_name.sif built successfully!"
  else
    log_message "Failed to build $image_name.sif. Please check the log file for more details."
  fi
done

log_message "Singularity images build process completed."
