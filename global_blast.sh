#!/bin/bash

log_step() {
    local step_name="$1"
    echo "$step_name - $(date)" >> "$log_file"
}

global_blast() {
    log_step "Running Global BLAST for clusters with no local hits"

    blast_global_output="${vsearch_dir}/${input_fastq_no_extension}.q${quality_score}.l${raw_min_length}.L${raw_max_length}.${primer_name}_linked_unlinked.sub.dd.att.con.sub.no_local.global.csv"
    log_step "Global BLAST output: $blast_global_output"

    if [ -z "$global_filter_taxids" ]; then
        log_step "Running global BLAST without taxid filter"
        $blastn_call -db $(basename "$global_blastdb") \
        -num_threads "$n_threads" \
        -outfmt "6 delim=, std qlen slen staxids sscinames scomnames sskingdoms" \
        -max_target_seqs 50 \
        -out "$blast_global_output" \
        -query "${vsearch_con_sub_no_local}" || { echo "Global BLAST command failed"; exit 1; }
    else
        log_step "Running global BLAST with taxid filter"
        $blastn_call -db $(basename "$global_blastdb") \
        -taxidlist "../$(basename "$global_filter_taxids")" \
        -num_threads "$n_threads" \
        -outfmt "6 delim=, std qlen slen staxids sscinames scomnames sskingdoms" \
        -max_target_seqs 50 \
        -out "$blast_global_output" \
        -query "${vsearch_con_sub_no_local}" || { echo "Global BLAST command with taxid filter failed"; exit 1; }
    fi
}

process_global_blast_output() {
    log_step "Processing Global BLAST output"

    output_cluster_filepath_global="${blast_dir}/blast_global.${primer_name}.vs-cluster-${vsearch_perc_id}.csv"
    output_taxa_filepath_global="${blast_dir}/blast_global.${primer_name}.vs-cluster-${vsearch_perc_id}.taxa.csv"
    clusters_matched_filepath_global="${blast_dir}/blast_global.${primer_name}.vs-cluster-${vsearch_perc_id}.clusters_matched.txt"

    log_step "Running R script for global BLAST output processing"
    $r_env_call /local/home/shannond/sdickey_minion_pipeline/scripts/R/blast_compile_results_CL.R \
    --cluster_key "$vsearch_uc" \
    --barcode_key "$reads_tab" \
    --conseq_key "$vsearch_con_sub_tab_edit" \
    --blast_results "$blast_global_output" \
    --sample_key "$sample_key" \
    --primer_name "$primer_name" \
    --min_alignment_length "$min_alignment_length" \
    --local_global "global" \
    --output_cluster_filepath "$output_cluster_filepath_global" \
    --output_taxa_filepath "$output_taxa_filepath_global" \
    --clusters_matched_filepath "$clusters_matched_filepath_global" \
    --taxdump_filepath "$ncbi_lineages_loc" || { echo "R script for global blast failed"; exit 1; }
}

# Usage: source global_blast.sh
