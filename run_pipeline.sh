#!/bin/bash

# Run the main pipeline script with the required parameters
./main.sh \
--input_fastq "/local/home/shannond/sdickey_minion_pipeline/data/input_fastq/TidalCycle_combined_labeled.q10.mifish_linked_unlinked.dd.con.sub.fasta" \
--output_dir "/local/home/shannond/sdickey_minion_pipeline/completed_pipeline" \
--quality_score 10 \
--raw_min_length 200 \
--raw_max_length 350 \
--primer_name "MiFish" \
--primer_file "/local/home/shannond/sdickey_minion_pipeline/primers/mifish_primer.fasta" \
--cutadapt_error_rate 0.2 \
--min_length 150 \
--max_length 250 \
--skip_porechop "no" \
--vsearch_perc_id 0.95 \
--local_blastdb "/local/home/shannond/sdickey_minion_pipeline/blastdb/local_blast_db" \
--global_blastdb "/local/home/shannond/sdickey_minion_pipeline/blastdb/global_blast_db/nt_euk" \
--global_filter_taxids "/local/home/shannond/sdickey_minion_pipeline/chordata_no_human.txids" \
--min_alignment_length 90 \
--sample_key "/local/home/shannond/sdickey_minion_pipeline/Tuckerton_Sample_Sheet.xlsx" \
--positive_control_list "/local/home/shannond/sdickey_minion_pipeline/species_lists/positive_control_species.txt" \
--contaminant_list "/local/home/shannond/sdickey_minion_pipeline/species_lists/contaminant_species.txt" \
--ncbi_lineages_loc "/local/home/shannond/sdickey_minion_pipeline/ncbi_lineages_2025-03-10.csv" \
--n_threads 36 \
--os "linux"
