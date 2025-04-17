#!/bin/bash

log_step() {
    local step_name="$1"
    echo "$step_name - $(date)" >> "$log_file"
}

local_blast() {
    log_step "Starting Local BLAST"

    blast_output="${vsearch_dir}/${input_fastq_no_extension}.q${quality_score}.l${raw_min_length}.L${raw_max_length}.${primer_name}_linked_unlinked.sub.dd.att.con.sub.csv"
    log_step "Defining BLAST output: $blast_output"

    log_step "Running blastn"
    $blastn_call \
    -db "/local/home/shannond/sdickey_minion_pipeline/blastdb/local_blast_db \
    -num_threads "$n_threads" \
    -outfmt "6 delim=, std qlen slen staxids sscinames scomnames sskingdoms" \
    -max_target_seqs 50 \
    -out "$blast_output" \
    -query "$vsearch_con_sub" || { echo "blastn command failed"; exit 1; }
}

process_blast_output() {
    log_step "Processing Local BLAST output"

    output_cluster_filepath="${blast_dir}/blast_local.${primer_name}.vs-cluster-${vsearch_perc_id}.csv"
    output_taxa_filepath="${blast_dir}/blast_local.${primer_name}.vs-cluster-${vsearch_perc_id}.taxa.csv"
    clusters_matched_filepath="${blast_dir}/blast_local.${primer_name}.vs-cluster-${vsearch_perc_id}.clusters_matched.txt"

    log_step "Running R script for BLAST output processing"
    $r_env_call /local/home/shannond/sdickey_minion_pipeline/scripts/R/blast_compile_results_CL.R \
    --cluster_key "$vsearch_uc" \
    --barcode_key "$reads_tab" \
    --conseq_key "$vsearch_con_sub_tab_edit" \
    --blast_results "$blast_output" \
    --sample_key "$sample_key" \
    --primer_name "$primer_name" \
    --min_alignment_length "$min_alignment_length" \
    --local_global "local" \
    --output_cluster_filepath "$output_cluster_filepath" \
    --output_taxa_filepath "$output_taxa_filepath" \
    --clusters_matched_filepath "$clusters_matched_filepath" \
    --taxdump_filepath "$ncbi_lineages_loc" || { echo "R script failed"; exit 1; }
}

# Usage: source blast_processing.sh
