#!/bin/bash

# Log input parameters
log_parameters() {
    log_file="$output_dir/_input_parameters.txt"
    > "$log_file"  # Clear the file

    for ((i = 0; i < ${#original_args[@]}; i+=2)); do
        key="${original_args[i]}"
        value="${original_args[i+1]}"
        echo "Processing: $key $value"
        printf "Parameter: %-25s Value: %-20s\n" "$key" "$value" >> "$log_file"
    done

    echo "Log file created: $log_file"
}
