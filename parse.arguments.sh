#!/bin/bash

# Parse command-line arguments
parse_arguments() {
    while [[ "$#" -gt 0 ]]; do
        case $1 in
            --input_fastq) input_fastq="$2"; shift ;;
            --output_dir) output_dir="$2"; shift ;;
            --quality_score) quality_score="$2"; shift ;;
            --raw_min_length) raw_min_length="$2"; shift ;;
            --raw_max_length) raw_max_length="$2"; shift ;;
            --primer_name) primer_name="$2"; shift ;;
            --primer_file) primer_file="$2"; shift ;;
            --cutadapt_error_rate) cutadapt_error_rate="$2"; shift ;;
            --min_length) min_length="$2"; shift ;;
            --max_length) max_length="$2"; shift ;;
            --skip_porechop) skip_porechop="$2"; shift ;;
            --vsearch_perc_id) vsearch_perc_id="$2"; shift ;;
            --local_blastdb) local_blastdb="$2"; shift ;;
            --global_blastdb) global_blastdb="$2"; shift ;;
            --global_filter_taxids) global_filter_taxids="$2"; shift ;;
            --min_alignment_length) min_alignment_length="$2"; shift ;;
            --sample_key) sample_key="$2"; shift ;;
            --positive_control_list) positive_control_list="$2"; shift ;;
            --contaminant_list) contaminant_list="$2"; shift ;;
            --tag_jump_rate) tag_jump_rate="$2"; shift ;;
            --ncbi_lineages_loc) ncbi_lineages_loc="$2"; shift ;;
            --n_threads) n_threads="$2"; shift ;;
            --os) os="$2"; shift ;;
            --help) usage; exit 0 ;;
            *) log_error "Unknown parameter passed: $1"; exit 1 ;;
        esac
        shift
    done
}
