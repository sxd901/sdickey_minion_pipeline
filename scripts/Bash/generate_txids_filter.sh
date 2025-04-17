#!/bin/bash
set -euo pipefail

# Add EDirect to PATH in case it was just installed
export PATH="${HOME}/edirect:${PATH}"

# ---------------------------------------------------------------------- #
# Run code for global blast filter out for only chordata_no_human.txids #
# ---------------------------------------------------------------------- #

# Define logging functions
log_info() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') [INFO] $1"
}

log_error() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') [ERROR] $1"
}

log_info "Starting the script execution."

# Set environment variables
iwd="/local/home/shannond/sdickey_minion_pipeline/blastdb"
global_blastdb="${iwd}/global_blast_db/nt_euk"
chordata_no_human_file="${iwd}/global_blast_db/taxid_db/chordata_no_human.txids"

# Check if chordata_no_human.txids exists
log_info "Checking if chordata_no_human.txids exists at $chordata_no_human_file..."
if [ ! -f "$chordata_no_human_file" ]; then
    log_info "chordata_no_human.txids does not exist. Generating now (requires internet)."

    cd "$iwd"

    # Check if EDirect tools are available
    if ! command -v esearch &>/dev/null; then
        log_info "Installing EDirect..."
        sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"
        export PATH="${HOME}/edirect:${PATH}"
    fi

    if ! command -v esearch &>/dev/null; then
        log_error "EDirect installation failed. Exiting script."
        exit 1
    fi
    log_info "EDirect installed successfully."

    # Fetch Chordata taxids
    log_info "Fetching Chordata taxids..."
    esearch -db taxonomy -query "Chordata[Subtree]" | \
    efetch -format docsum | \
    xtract -pattern Taxon -element TaxId > chordata.txids 2> chordata.log

    if [ ! -s chordata.txids ]; then
        log_error "Failed to fetch Chordata taxids. Please check chordata.log."
        exit 1
    fi
    log_info "Chordata taxids fetched successfully."

    # Fetch Hominidae (Human-related) taxids
    log_info "Fetching Hominidae (human-related) taxids..."
    esearch -db taxonomy -query "Hominidae[Subtree]" | \
    efetch -format docsum | \
    xtract -pattern Taxon -element TaxId > hominidae.txids 2> hominidae.log

    if [ ! -s hominidae.txids ]; then
        log_error "Failed to fetch Human taxids. Please check hominidae.log."
        exit 1
    fi
    log_info "Human taxids fetched successfully."

    # Show preview
    log_info "Chordata taxids (first 10 lines):"
    head -n 10 chordata.txids
    log_info "Hominidae taxids (first 10 lines):"
    head -n 10 hominidae.txids

    # Sort the files
    log_info "Sorting taxid files..."
    sort chordata.txids -o chordata.txids
    sort hominidae.txids -o hominidae.txids

    # Remove human taxids from Chordata list
    log_info "Filtering human taxids out of chordata..."
    comm -23 chordata.txids hominidae.txids > chordata_no_human.txids

    if [ ! -s chordata_no_human.txids ]; then
        log_error "Filtering failed. Output file is empty."
        exit 1
    fi
    log_info "Filtering complete."

    log_info "chordata_no_human.txids (first 10 lines):"
    head -n 10 chordata_no_human.txids

    # Move final output to target location
    log_info "Moving result to $chordata_no_human_file"
    mv chordata_no_human.txids "$chordata_no_human_file"
    if [ $? -ne 0 ]; then
        log_error "Failed to move output file."
        exit 1
    fi

    # Clean up
    log_info "Cleaning up temporary files..."
    rm -f chordata.txids hominidae.txids
    log_info "Cleanup complete."

else
    log_info "chordata_no_human.txids already exists. Skipping generation."
fi

# Return to original directory
log_info "Returning to working directory: $iwd"
cd "$iwd"

log_info "Script execution completed successfully."
