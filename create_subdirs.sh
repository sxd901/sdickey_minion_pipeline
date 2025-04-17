#!/bin/bash

# Create required directories
create_subdirs() {
    mkdir -p "${output_dir}/intermediate_fastas"
    mkdir -p "${output_dir}/read_length_distribs"
    mkdir -p "${output_dir}/vsearch"
    mkdir -p "${output_dir}/blast_summary"
}
