#!/bin/bash

log_step() {
    local step_name="$1"
    echo "$step_name - $(date)" >> "$log_file"
}

cluster_reads() {
    log_step "Clustering reads using vsearch"

    vsearch_cen="${vsearch_dir}/${input_fastq_no_extension}.q${quality_score}.l${raw_min_length}.L${raw_max_length}.${primer_name}_linked_unlinked.sub.dd.att.cen.fasta"
    vsearch_con="${vsearch_dir}/${input_fastq_no_extension}.q${quality_score}.l${raw_min_length}.L${raw_max_length}.${primer_name}_linked_unlinked.sub.dd.att.con.fasta"
    vsearch_con_sub="${vsearch_dir}/${input_fastq_no_extension}.q${quality_score}.l${raw_min_length}.L${raw_max_length}.${primer_name}_linked_unlinked.sub.dd.att.con.sub.fasta"
    vsearch_uc="${vsearch_dir}/${input_fastq_no_extension}.q${quality_score}.l${raw_min_length}.L${raw_max_length}.${primer_name}_linked_unlinked.sub.dd.att.clusters.uc"

    log_step "Running vsearch clustering"
    if [ "$os" = "linux" ]; then
        vsearch_call="singularity exec /local/home/shannond/sdickey_minion_pipeline/images/vsearch.sif vsearch"
    elif [ "$os" = "mac" ]; then
        vsearch_call="vsearch"
    fi

    $vsearch_call --threads "$n_threads" --cluster_fast \
    "$reads_fasta_att" \
    --centroids "$vsearch_cen" \
    --id "$vsearch_perc_id" \
    --consout "$vsearch_con" \
    --uc "$vsearch_uc" || { echo "VSEARCH clustering failed"; exit 1; }

    log_step "Appending read count of vsearch_con"
    get_counts_append "$vsearch_con" "$read_counts" "$os"

    log_step "Removing clusters under 5 reads"
    seqkit grep -p 'seqs=([6-9]|[1-9]\d+)$' "$vsearch_con" -o "$vsearch_con_sub"

    log_step "Tabulating vsearch_con_sub"
    seqkit fx2tab -l "$vsearch_con_sub" -o "$vsearch_con_sub_tab"
}

# Usage: source clustering.sh
