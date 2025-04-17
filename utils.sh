#!/bin/bash

# Logging function to log messages
log_info() {
    logger -t "eDNA_pipeline" "[INFO] $1"
}

log_error() {
    logger -t "eDNA_pipeline" "[ERROR] $1"
}

# Validate if a file exists
validate_file() {
    local file_path=$1
    if [ ! -f "$file_path" ]; then
        log_error "File does not exist: $file_path"
        exit 1
    fi
}

# Validate if a directory exists
validate_directory() {
    local dir_path=$1
    if [ ! -d "$dir_path" ]; then
        log_error "Directory does not exist: $dir_path"
        exit 1
    fi
}

# Validate if a parameter is a number
validate_positive_integer() {
    local num=$1
    if ! [[ "$num" =~ ^[0-9]+$ ]]; then
        log_error "Parameter must be a positive integer: $num"
        exit 1
    fi
}

# Validate if a parameter is a float between 0 and 1
validate_float_range() {
    local num=$1
    if (( $(echo "$num < 0" | bc -l) )) || (( $(echo "$num > 1" | bc -l) )); then
        log_error "Parameter must be between 0 and 1: $num"
        exit 1
    fi
}

# Validate if file has the expected extension
validate_extension() {
    local file_path=$1
    local ext=$2
    if [[ "$file_path" != *.$ext ]]; then
        log_error "File $file_path must have .$ext extension"
        exit 1
    fi
}
