#!/bin/bash

# Logging functions
log_info() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') [INFO] $1"
}

log_error() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') [ERROR] $1"
}

log_info "Starting the BLAST database update and download process..."

# Check available databases
log_info "Checking available BLAST databases..."
update_blastdb.pl --showall [*]
if [ $? -ne 0 ]; then
    log_error "Failed to check available BLAST databases."
    exit 1
fi
log_info "Successfully checked available BLAST databases."

# With Singularity: Check available databases
log_info "Using Singularity to check available BLAST databases..."
singularity exec /local/home/shannond/sdickey_minion_pipeline/images/blast.sif update_blastdb.pl --showall [*]
if [ $? -ne 0 ]; then
    log_error "Failed to check available BLAST databases with Singularity."
    exit 1
fi
log_info "Successfully checked available BLAST databases with Singularity."

# Change to the appropriate directory for database storage
log_info "Changing directory to /local/home/shannond/sdickey_minion_pipeline/blastdb/global_blast_db..."
cd /local/home/shannond/sdickey_minion_pipeline/blastdb/global_blast_db || { log_error "Failed to change directory."; exit 1; }

# Create a new directory for the nt_euk database
log_info "Creating directory for nt_euk database..."
mkdir -p nt_euk
cd nt_euk || { log_error "Failed to change to global_db_nt_euk directory."; exit 1; }

# Download and decompress the nt_euk database
log_info "Downloading and decompressing the nt_euk database..."
update_blastdb.pl --decompress nt_euk
if [ $? -ne 0 ]; then
    log_error "Failed to download and decompress nt_euk database."
    exit 1
fi
log_info "Successfully downloaded and decompressed the nt_euk database."

# With Singularity: Download and decompress the nt_euk database
log_info "Using Singularity to download and decompress the nt_euk database..."
singularity exec /local/home/shannond/sdickey_minion_pipeline/images/blast.sif update_blastdb.pl --decompress nt_euk
if [ $? -ne 0 ]; then
    log_error "Failed to download and decompress nt_euk database with Singularity."
    exit 1
fi
log_info "Successfully downloaded and decompressed the nt_euk database with Singularity."

# Check the integrity of the nt_euk database
log_info "Checking the integrity of the nt_euk database..."
blastdbcheck -db nt_euk -dbtype nucl
if [ $? -ne 0 ]; then
    log_error "Failed to check integrity of nt_euk database."
    exit 1
fi
log_info "Successfully checked the integrity of the nt_euk database."

# With Singularity: Check the integrity of the nt_euk database
log_info "Using Singularity to check the integrity of the nt_euk database..."
singularity exec /local/home/shannond/sdickey_minion_pipeline/images/blast.sif blastdbcheck -db nt_euk -dbtype nucl
if [ $? -ne 0 ]; then
    log_error "Failed to check integrity of nt_euk database with Singularity."
    exit 1
fi
log_info "Successfully checked the integrity of the nt_euk database with Singularity."

# Compress the folder for transfer
log_info "Compressing the nt_euk_new directory for transfer..."
cd /local/home/shannond/sdickey_minion_pipeline/blastdb/global_blast_db || { log_error "Failed to change directory."; exit 1; }
tar -czvf nt_euk_new.tar.gz global_db_nt_euk
if [ $? -ne 0 ]; then
    log_error "Failed to compress the nt_euk_new directory."
    exit 1
fi
log_info "Successfully compressed the nt_euk_new directory for transfer."

# Final log message
log_info "BLAST database update and download process completed successfully."
