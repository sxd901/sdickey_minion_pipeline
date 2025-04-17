#!/bin/bash

# Function to check if a script exists
check_script() {
    local script_path="$1"
    local script_name="$2"
    if [ ! -f "$script_path" ]; then
        echo "$script_name is missing. Exiting."
        logger "ERROR: $script_name is missing. Exiting."
        exit 1
    fi
}

# Check if the required R scripts are present
check_script "/local/home/shannond/sdickey_minion_pipeline/scripts/R/run_cutadapt_CL.R" "run_cutadapt_CL.R"
check_script "/local/home/shannond/sdickey_minion_pipeline/scripts/R/process_read_distribs_CL.R" "process_read_distribs_CL.R"
check_script "/local/home/shannond/sdickey_minion_pipeline/scripts/R/blast_compile_results_CL.R" "blast_compile_results_CL.R"
check_script "/local/home/shannond/sdickey_minion_pipeline/scripts/R/process_positive_negatives_CL.R" "process_positive_negatives_CL.R"
