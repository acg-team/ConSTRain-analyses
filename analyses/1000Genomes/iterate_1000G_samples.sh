#!/usr/bin/env bash
# Iterate over a file containing sample id and FTP links to aligned sequencing data for 1000 Genomes samples
# Meant to figure out how to do this so we can later adapt it for an sbatch script using sbatch-array

set -euo pipefail

# IMPORTANT: 'CSV' should:
#   - be a comme-separated file
#   - contain single a header line
#   - contain the following columns: <sequential index>,triadID,sampleID,role,pgx_id,sex,ENA_FILE_PATH
CSV="../../data/1000Genomes/1000_genomes_triad_pedigrees.csv"


for ID in $(seq 3 5);
do  
    LINE=$(grep -m 1 "^${ID}" ${CSV});    
    TRIAD=$(echo $LINE | cut -d "," -f 2);
    SAMPLE=$(echo $LINE | cut -d "," -f 3);
    SEX=$(echo $LINE | cut -d "," -f 6);
    FTP=$(echo $LINE | cut -d "," -f 7);
    echo ${ID} ${SAMPLE} ${TRIAD} ${SEX} ${FTP};
done;