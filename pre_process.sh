#!/bin/bash

log_step() {
    local step_name="$1"
    echo "$step_name - $(date)" >> "$log_file"
}

preprocess_reads() {
    log_step "Starting preprocess reads"
    
    reads_linked_unlinked="${fasta_dir}/${input_fastq_no_extension}.q${quality_score}.l${raw_min_length}.L${raw_max_length}.${primer_name}_linked_unlinked.sub.fastq"
    reads_linked_unlinked_dd="${fasta_dir}/${input_fastq_no_extension}.q${quality_score}.l${raw_min_length}.L${raw_max_length}.${primer_name}_linked_unlinked.sub.dd.fastq"
    
    log_step "Combining linked and unlinked reads"
    cat "$reads_linked_dd" "$reads_unlinked_sub" > "$reads_linked_unlinked"

    log_step "Removing duplicates with seqkit"
    $seqkit_call rmdup "$reads_linked_unlinked" -o "$reads_linked_unlinked_dd" || { echo "Error removing duplicates"; exit 1; }
}

tabulate_reads() {
    log_step "Tabulating reads"
    
    reads_fasta="${fasta_dir}/${input_fastq_no_extension}.q${quality_score}.l${raw_min_length}.L${raw_max_length}.${primer_name}_linked_unlinked.sub.dd.fasta"
    reads_fasta_att="${fasta_dir}/${input_fastq_no_extension}.q${quality_score}.l${raw_min_length}.L${raw_max_length}.${primer_name}_linked_unlinked.sub.dd.att.fasta"
    reads_tab="${vsearch_dir}/${input_fastq_no_extension}.q${quality_score}.l${raw_min_length}.L${raw_max_length}.${primer_name}_linked_unlinked.sub.dd.att.txt"

    log_step "Converting reads to FASTA format"
    $seqkit_call fq2fa "$reads_linked_unlinked_dd" -o "$reads_fasta"

    log_step "Dropping unnecessary attributes"
    awk '/^>/ {split($0, a, " "); printf("%s", a[1]); for (i=2; i<=length(a); i++) { if (match(a[i], /^barcode=|^start_time=/)) printf(" %s", a[i]); } printf("\n"); next} {print}' \
    "$reads_fasta" > "$reads_fasta_att"

    log_step "Tabulating FASTA file"
    $seqkit_call fx2tab -l "$reads_fasta_att" -o "$reads_tab"
}

# Usage: source preprocess.sh
