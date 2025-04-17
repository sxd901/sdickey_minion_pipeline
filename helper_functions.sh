#!/bin/bash

# Get read count and append to file
get_counts_append() {
    local input_fastx="$1"
    local output_file="$2"
    local os="$3"

    logger "Counting number of reads in $input_fastx"
    echo "Counting number of reads in $input_fastx"

    if [ "$os" = "linux" ]; then
        input_fastx_count=$(singularity exec /local/home/shannond/sdickey_minion_pipeline/images/seqkit.sif \
        seqkit stats "$input_fastx" | \
        awk '{print $4}' | grep -o '[0-9,]*' | tr -d ',')
    elif [ "$os" = "mac" ]; then
        input_fastx_count=$(seqkit stats "$input_fastx" | \
        awk '{print $4}' | grep -o '[0-9,]*' | tr -d ',') 
    fi  

    logger "Number of reads in $input_fastx: $input_fastx_count"
    echo "Number of reads in $input_fastx: $input_fastx_count"

    echo -e "$input_fastx\t$input_fastx_count" >> "$output_file"
}

# Get read length distribution
get_read_len_distribs() {
    local input_fastx="$1"
    local output_file="$2"
    local os="$3"

    logger "Getting read length distribution for $input_fastx"
    echo "Getting read length distribution for $input_fastx"

    if [ "$os" = "linux" ]; then
        bbmap_readlen_call="singularity exec /local/home/shannond/sdickey_minion_pipeline/images/bbmap.sif readlength.sh"
    elif [ "$os" = "mac" ]; then
        bbmap_readlen_call="~/bbmap/readlength.sh"
    fi

    $bbmap_readlen_call in="$input_fastx" out="$output_file"
}

# Generate read histogram plots
generate_read_hist_plots() {
    local input_file="$1"
    local output_file="$2"
    local min_length="$3"
    local max_length="$4"
    local os="$5"

    logger "Running generate_read_hist_plots for input_file = $input_file, outputfile = $output_file, min_length = $min_length, max_length = $max_length"
    $r_env_call /local/home/shannond/sdickey_minion_pipeline/scripts/R/process_read_distribs_CL.R \
    --distrib_file="$input_file" --output_file="$output_file" \
    --hmin="$min_length" --hmax="$max_length"
}

# Subset vsearch fasta
vsearch_fasta_subset() {
    local ids_file="$1"
    local fasta_file="$2"
    local output_file="$3"
    local os="$4"

    logger "Subsetting vsearch fasta"
    sed 's/^/centroid=/' $ids_file > temp-newids.txt
    sed '/^>/ s/;/ /g' $fasta_file > temp-vs-clusters.fasta
    $seqkit_call grep --invert-match -f temp-newids.txt temp-vs-clusters.fasta -o temp-vs-clusters.sub.fasta
    sed '/^>/ s/ /;/g' temp-vs-clusters.sub.fasta > "$output_file"

    # Cleanup
    rm temp-newids.txt temp-vs-clusters.fasta temp-vs-clusters.sub.fasta
}
