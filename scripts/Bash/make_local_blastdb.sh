#!/bin/bash

# This script will create a local blast database from a fasta file that has the taxid=; attribute for each sequence
# Note: requires internet to download ncbi taxonomy file
# If using Linux, it assumes you are using the images/blast.sif singularity image
# If using Mac, it assumes you have blast+ installed

# $1 Input fasta that is the local reference database, must have taxid=; attribute for each seq
# $2 Output blastdb dir for the local blast database e.g. blastdb/florida-mifish/florida-mifish
# $3 taxdb directory: directory where the taxdb.tar.gz file is located; If left blank as "", it will download a new one (needs internet)
# $4 os: operating system, either "linux" or "mac"

# Example call:
# ./make_local_blastdb.sh \
# blastdb/tuckerton-mifish/tuckerton-mifish \
# "taxdb" "linux"

# Assign parameters to variables
input_fasta="/local/home/shannond/molra_xprize_edna_minion_pipeline-main/data/tuckerton_run1/TidalCycle_combined_labeled.q10.mifish_linked_unlinked.dd.con.sub.fasta"
out_blastdb_dir="/local/home/shannond/sdickey_minion_pipeline/blastdb/local_blast_db/local_blast_db"
taxdb_dir=""
os="linux"

# Check if the input fasta exists
if [ ! -f "$input_fasta" ]; then
    echo "Input fasta file does not exist!"
    exit 1
fi

### Create a temporary directory
temp_dir=$(mktemp -d /tmp/blastdb.XXXXXX)

# Get filename without extension (e.g., test/test2/file.fastq -> file)
input_fasta_no_extension=$(basename "${input_fasta%.*}")

# Generate a taxid map file
awk -F "[>;= ]" '/^>/{for(i=1; i<=NF; i++) if($i == "taxid") print $2, $(i+1)}' \
"$input_fasta" > \
"${temp_dir}/${input_fasta_no_extension}.taxid.txt"

# Make sure the temporary directory exists
mkdir -p "$temp_dir"

# Run makeblastdb
if [ "$os" == "linux" ]; then
    singularity exec /local/home/shannond/molra_xprize_edna_minion_pipeline-main/images/blast.sif \
    makeblastdb -in "$input_fasta" \
    -parse_seqids \
    -taxid_map "${temp_dir}/${input_fasta_no_extension}.taxid.txt" \
    -dbtype nucl \
    -out "$out_blastdb_dir"
elif [ "$os" == "mac" ]; then
    makeblastdb -in "$input_fasta" \
    -parse_seqids \
    -taxid_map "${temp_dir}/${input_fasta_no_extension}.taxid.txt" \
    -dbtype nucl \
    -out "$out_blastdb_dir"
fi

# Handle the taxdb directory
if [ "$taxdb_dir" == "" ]; then
    echo "Downloading taxonomy file"
    # Get the taxonomy file
    wget https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz
    # Extract the tar.gz file to the directory
    tar -xvzf taxdb.tar.gz -C $(dirname "$out_blastdb_dir")
    # Optional: remove the taxdb.tar.gz file
    rm taxdb.tar.gz
else
    # Copy taxdb files from the provided directory
    cp $taxdb_dir/* $(dirname "$out_blastdb_dir")
fi
