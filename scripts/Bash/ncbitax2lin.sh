#!/bin/bash

# Define logging functions
log_info() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') [INFO] $1"
}

log_error() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') [ERROR] $1"
}

log_info "Starting script execution."

# Define the target directory for taxdump and final files
TAXDUMP_DIR="/local/home/shannond/sdickey_minion_pipeline/blastdb/NCBI_tax_id_lineage_key/taxdump"
FINAL_CSV_DIR="/local/home/shannond/sdickey_minion_pipeline/blastdb/NCBI_tax_id_lineage_key"
DATE=$(date '+%Y-%m-%d')  # Dynamic date for the file names

# Step 1: Create the 'taxdump' directory if it doesn't exist
log_info "Creating taxdump directory at $TAXDUMP_DIR if it doesn't exist..."
mkdir -p "$TAXDUMP_DIR"
if [ $? -ne 0 ]; then
    log_error "Failed to create taxdump directory at $TAXDUMP_DIR"
    exit 1
fi
log_info "Taxdump directory created successfully."

# Step 2: Download the taxdump.tar.gz file into the 'taxdump' directory
log_info "Downloading taxdump.tar.gz from NCBI..."
wget -N ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz -P "$TAXDUMP_DIR"
if [ $? -ne 0 ]; then
    log_error "Failed to download taxdump.tar.gz"
    exit 1
fi
log_info "Successfully downloaded taxdump.tar.gz."

# Step 3: Extract the taxdump.tar.gz file
log_info "Extracting taxdump.tar.gz into $TAXDUMP_DIR..."
tar zxf "$TAXDUMP_DIR/taxdump.tar.gz" -C "$TAXDUMP_DIR"
if [ $? -ne 0 ]; then
    log_error "Failed to extract taxdump.tar.gz"
    exit 1
fi
log_info "Successfully extracted taxdump.tar.gz."

# Step 4: Ensure the final output directory exists
log_info "Ensuring final CSV output directory exists: $FINAL_CSV_DIR"
mkdir -p "$FINAL_CSV_DIR"
if [ $? -ne 0 ]; then
    log_error "Failed to create final CSV directory at $FINAL_CSV_DIR"
    exit 1
fi

# Step 5: Check for required files before proceeding
log_info "Checking for required files..."
if [ ! -f "$TAXDUMP_DIR/nodes.dmp" ]; then
    log_error "nodes.dmp not found in $TAXDUMP_DIR"
    exit 1
fi
if [ ! -f "$TAXDUMP_DIR/names.dmp" ]; then
    log_error "names.dmp not found in $TAXDUMP_DIR"
    exit 1
fi
log_info "Required files (nodes.dmp, names.dmp) found."

# Step 6: Set the Singularity image path
SINGULARITY_IMAGE="/local/home/shannond/sdickey_minion_pipeline/images/ncbitax2lin.sif"
CSV_FILE="${FINAL_CSV_DIR}/ncbi_lineages_${DATE}.csv.gz"  # Explicit output file path

log_info "Attempting to run ncbitax2lin via Singularity to convert taxonomy files..."
singularity exec "$SINGULARITY_IMAGE" \
    ncbitax2lin --nodes-file "$TAXDUMP_DIR/nodes.dmp" \
                --names-file "$TAXDUMP_DIR/names.dmp" \
                --output "$CSV_FILE"
if [ $? -ne 0 ]; then
    log_error "Failed to run ncbitax2lin using Singularity"
    exit 1
fi
log_info "Successfully ran ncbitax2lin using Singularity."

# Step 7: Check if the CSV file was created and decompress it
log_info "Checking if the CSV file was created at $CSV_FILE..."
if [ ! -f "$CSV_FILE" ]; then
    log_error "CSV file not found at $CSV_FILE"
    exit 1
fi

log_info "Decompressing $CSV_FILE..."
gunzip "$CSV_FILE"
if [ $? -ne 0 ]; then
    log_error "Failed to decompress $CSV_FILE"
    exit 1
fi
log_info "Successfully decompressed $CSV_FILE."

# Step 8: Recompress the CSV file (keeping the original CSV)
log_info "Recompressing $FINAL_CSV_DIR/ncbi_lineages_${DATE}.csv..."
gzip -c "${FINAL_CSV_DIR}/ncbi_lineages_${DATE}.csv" > "${FINAL_CSV_DIR}/ncbi_lineages_${DATE}.csv.gz"
if [ $? -ne 0 ]; then
    log_error "Failed to recompress $FINAL_CSV_DIR/ncbi_lineages_${DATE}.csv"
    exit 1
fi
log_info "Successfully recompressed $CSV_FILE."

log_info "Script execution completed successfully."
