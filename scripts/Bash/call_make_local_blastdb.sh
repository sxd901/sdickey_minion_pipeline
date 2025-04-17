#!/bin/bash

# ------------------------------------------------------- #

# Logging functions
log_info() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') [INFO] $1"
}

log_error() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') [ERROR] $1"
}

log_info "Starting script execution..."

# Change to home directory
cd ~
log_info "Changed directory to home (~)"

# Make make_local_blastdb.sh executable
chmod +x /local/home/shannond/molra_xprize_edna_minion_pipeline-main/scripts/make_local_blastdb.sh
if [ $? -ne 0 ]; then
    log_error "Failed to make make_local_blastdb.sh executable"
    exit 1
fi
log_info "Made make_local_blastdb.sh executable"

# Fix line endings (if calling script from gives this error: /bin/bash^M: bad interpreter: No such file or directory)
log_info "Fixing line endings in make_local_blastdb.sh"
sed -i -e 's/\r$//' /local/home/shannond/molra_xprize_edna_minion_pipeline-main/scripts/make_local_blastdb.sh
if [ $? -ne 0 ]; then
    log_error "Failed to fix line endings in make_local_blastdb.sh"
    exit 1
fi
log_info "Line endings fixed successfully"

# ------------------------------------------------------- #

# Replace 'taxid=NA;' with 'taxid=;' in the FASTA file
log_info "Replacing 'taxid=NA;' with 'taxid=;' in the reference FASTA file"
sed -i 's/taxid=NA;/taxid=;/g' /local/home/shannond/molra_xprize_edna_minion_pipeline-main/data/tuckerton_run1/TidalCycle_combined_labeled.q10.mifish_linked_unlinked.dd.con.sub.fasta
if [ $? -ne 0 ]; then
    log_error "Failed to replace 'taxid=NA;' in FASTA file"
    exit 1
fi
log_info "'taxid=NA;' replaced with 'taxid=;' successfully"

# Run the make_local_blastdb script
log_info "Running make_local_blastdb.sh to create BLAST database"
/local/home/shannond/molra_xprize_edna_minion_pipeline-main/scripts/make_local_blastdb.sh \
  /local/home/shannond/molra_xprize_edna_minion_pipeline-main/data/tuckerton_run1/TidalCycle_combined_labeled.q10.mifish_linked_unlinked.dd.con.sub.fasta \
  blastdb/MiFish_Tuckerton/MiFish_Tuckerton \
  "blastdb/taxdb" "linux"
if [ $? -ne 0 ]; then
    log_error "Failed to run make_local_blastdb.sh"
    exit 1
fi
log_info "make_local_blastdb.sh executed successfully"

# Final log message
log_info "Script execution completed successfully."
