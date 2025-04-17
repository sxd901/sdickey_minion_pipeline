#!/bin/bash

# Main pipeline script to source all modules and execute the steps

# Parse command-line arguments
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
        --ncbi_lineages_loc) ncbi_lineages_loc="$2"; shift ;;
        --n_threads) n_threads="$2"; shift ;;
        --os) os="$2"; shift ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
    shift
done

# Check for missing arguments
if [ -z "$input_fastq" ] || [ -z "$output_dir" ] || [ -z "$quality_score" ] || [ -z "$raw_min_length" ] || [ -z "$raw_max_length" ]; then
    echo "Missing required arguments. Exiting."
    exit 1
fi

# Initialize file name without extension
input_fastq_no_extension=$(basename "$input_fastq" .fastq)

# Record start time
start_time=$(date +%s.%N)

# Log the start of the pipeline
echo "Starting the pipeline at $(date)"

# Source all the modules (Ensure the paths to the scripts are correct)
source ./check_scripts.sh
source ./setup_variables.sh
source ./create_subdirs.sh
source ./log_parameters.sh
source ./helper_functions.sh
source ./process_reads.sh  # New module for processing the reads

# Initialize directories
fasta_dir="${output_dir}/fasta"
read_distrib_dir="${output_dir}/read_length_distribs"
mkdir -p "$fasta_dir" "$read_distrib_dir"

# Call functions to set up environment
check_scripts  # Check if all required scripts exist
setup_variables  # Set up OS-specific variables
create_subdirs  # Create necessary subdirectories
log_parameters  # Log input parameters

# Start processing with logging at each step

# --------------------------------------------------------------------------------------------- #
# Get n reads in input file
get_counts_append "$input_fastq" "$read_counts" "$os"

# --------------------------------------------------------------------------------------------- #
# Quality filter reads
echo "Quality filtering reads in $input_fastq to >=q$quality_score"
qfiltered_fastq="${output_dir}/${input_fastq_no_extension}.q${quality_score}.fastq"
echo "Running NanoFilt on $input_fastq"
run_nanofilt "$input_fastq" "$qfiltered_fastq" "$quality_score" "$os"

# Log the filtered read count
get_counts_append "$qfiltered_fastq" "$read_counts" "$os"

# --------------------------------------------------------------------------------------------- #
# Get initial read count distribution by barcode
echo "Getting initial read count distribution by barcode (of quality filtered reads)"
generate_barcode_distribution "$qfiltered_fastq" "$output_dir" "$input_fastq_no_extension"

# --------------------------------------------------------------------------------------------- #
# Length Filter reads
echo "Length filtering reads in $qfiltered_fastq between $raw_min_length and $raw_max_length"
length_filtered_fastq="${fasta_dir}/${input_fastq_no_extension}.q${quality_score}.l${raw_min_length}.L${raw_max_length}.fastq"
length_filtered_distrib="${read_distrib_dir}/${input_fastq_no_extension}.q${quality_score}.l${raw_min_length}.L${raw_max_length}.distrib.txt"

# Ensure the output directory exists
mkdir -p "$(dirname "$length_filtered_fastq")"

run_length_filter "$qfiltered_fastq" "$length_filtered_fastq" "$raw_min_length" "$raw_max_length"
get_read_len_distribs "$length_filtered_fastq" "$length_filtered_distrib" "$os"
get_counts_append "$length_filtered_fastq" "$read_counts" "$os"

# --------------------------------------------------------------------------------------------- #
# Trim Primers
echo "Detecting and trimming primers using cutadapt"
trimmed_reads_linked="${fasta_dir}/${input_fastq_no_extension}.q${quality_score}.l${raw_min_length}.L${raw_max_length}.${primer_name}_linked.fastq"
trimmed_reads_unlinked="${fasta_dir}/${input_fastq_no_extension}.q${quality_score}.l${raw_min_length}.L${raw_max_length}.${primer_name}_unlinked.fastq"

# Ensure the output directory exists
mkdir -p "$(dirname "$trimmed_reads_linked")"

run_cutadapt "$length_filtered_fastq" "$primer_file" "$trimmed_reads_linked" "$trimmed_reads_unlinked" "$cutadapt_error_rate" "$n_threads"

# Deduplicate linked reads and log counts
deduplicated_linked="${fasta_dir}/${input_fastq_no_extension}.q${quality_score}.l${raw_min_length}.L${raw_max_length}.${primer_name}_linked.dd.fastq"
seqkit rmdup "$trimmed_reads_linked" > "$deduplicated_linked"
get_counts_append "$deduplicated_linked" "$read_counts" "$os"

# --------------------------------------------------------------------------------------------- #
# Generate histogram for linked reads
echo "Generating histogram for linked reads"
linked_distrib="${read_distrib_dir}/${input_fastq_no_extension}.q${quality_score}.l${raw_min_length}.L${raw_max_length}.${primer_name}_linked.dd.distrib.txt"
linked_plot="${read_distrib_dir}/${input_fastq_no_extension}.q${quality_score}.l${raw_min_length}.L${raw_max_length}.${primer_name}_linked.dd.distrib.png"
generate_read_hist_plots "$linked_distrib" "$linked_plot" "$min_length" "$max_length" "$os"

# --------------------------------------------------------------------------------------------- #
# Remove IDs of linked reads from unlinked reads
echo "Removing IDs of linked reads from unlinked reads"
remove_linked_ids "$deduplicated_linked" "$trimmed_reads_unlinked"

# --------------------------------------------------------------------------------------------- #
# Porechop unlinked reads
if [ "$skip_porechop" = "no" ]; then
    echo "Porechopping unlinked reads"
    porechop_unlinked "$trimmed_reads_unlinked" "$fasta_dir" "$input_fastq_no_extension" "$quality_score" "$raw_min_length" "$raw_max_length" "$primer_name" "$os"
fi

# --------------------------------------------------------------------------------------------- #
# Combine linked and unlinked reads
echo "Combining linked and unlinked reads"
combine_linked_unlinked "$deduplicated_linked" "$porechopped_reads" "$fasta_dir" "$input_fastq_no_extension" "$quality_score" "$raw_min_length" "$raw_max_length" "$primer_name" "$os"

# --------------------------------------------------------------------------------------------- #
# Finish timer and calculate total runtime
end_time=$(date +%s.%N)
runtime=$(echo "$end_time - $start_time" | bc)

# Log the end of the pipeline and total runtime
echo "Pipeline finished at $(date)"
echo "Script finished in $runtime seconds"
