#!/bin/bash

# ---------------------------------------------------------------------- #

# Linux 
# note blast 2.16.0 doesn't have get_species_taxids.sh need to use 2.15.0

# install edirect (needed for get_species_taxids.sh)
singularity exec /local/home/shannond/sdickey_minion_pipeline/images/blast.sif \
sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)" && \
sh -c "$(wget -q https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -O -)" && \
export PATH=${HOME}/edirect:${PATH}

# get chordate ids
singularity exec /local/home/shannond/sdickey_minion_pipeline/images/blast.sif \
get_species_taxids.sh -t 7711 > chordata.txids

# also download human taxids and remove from chordata
singularity exec /local/home/shannond/sdickey_minion_pipeline/images/blast.sif \
get_species_taxids.sh -t 9604 > hominidae.txids

# rm hominidae from chordata
## first sort files
sort chordata.txids -o chordata.txids
sort hominidae.txids -o hominidae.txids
## then comm to rm human from chordata
comm -23 chordata.txids hominidae.txids > chordata_no_human.txids

# move to blastdb folder
mv chordata_no_human.txids ~/blastdb/

# rm others
rm chordata.txids hominidae.txids


