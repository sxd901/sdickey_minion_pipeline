#!/bin/bash

# Setup OS-specific variables
setup_variables() {
    if [ "$os" = "linux" ]; then
        seqkit_call="seqkit"
        blastn_call="$HOME/miniconda3/bin/blastn"
        r_env_call="singularity exec /local/home/shannond/sdickey_minion_pipeline/images/r_env.sif Rscript"
    elif [ "$os" = "mac" ]; then
        seqkit_call="seqkit"
        blastn_call="blastn"
        r_env_call="Rscript"
    fi
}
