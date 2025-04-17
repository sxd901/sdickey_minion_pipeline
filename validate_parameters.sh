#!/bin/bash

# Validate all the parameters
validate_parameters() {
    # Validate input fastq
    if [ -z "$input_fastq" ]; then
        log_error "Input fastq not provided. Exiting."
        exit 1
    fi
    validate_file "$input_fastq"
    validate_extension "$input_fastq" "fastq"

    # Validate output directory
    if [ -z "$output_dir" ]; then
        log_error "Output directory not provided. Exiting."
        exit 1
    fi
    validate_directory "$output_dir"

    # Validate quality score
    if [ -z "$quality_score" ]; then
        log_info "Quality score not provided, assuming 10."
        quality_score=10
    else
        validate_positive_integer "$quality_score"
    fi

    # Validate primer name
    if [ -z "$primer_name" ]; then
        log_error "Primer name not provided. Exiting."
        exit 1
    fi

    # Validate primer file
    if [ -z "$primer_file" ]; then
        log_error "Primer fasta file not provided. Exiting."
        exit 1
    fi
    validate_file "$primer_file"

    # Validate cutadapt error rate
    if [ -z "$cutadapt_error_rate" ]; then
        log_error "Cutadapt error rate not provided. Exiting."
        exit 1
    fi
    validate_float_range "$cutadapt_error_rate"

    # Validate other parameters (similar logic can be applied to others)
    # ... You can add more validation as needed for the remaining parameters ...
}
