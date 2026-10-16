#!/usr/bin/env bash

set -euo pipefail

UNIPROTSIGNALOUT="uniprot_signal_peptides.csv"

# Remember where we started
CURRENT_DIR="$(pwd)"

# Parent directory where genename2uniqproteinhaplo.sh expects to run
cd ..

missing_count=0
while IFS=',' read -r gene np_id nm_id aa_length uniprot signal_start signal_end signal_length status || [ -n "$gene" ]; do

    # Skip header
    [ "$gene" = "gene" ] && continue

    # Define paths for your two target files per gene
    file1="results/${gene}/${gene}_protein_haplotypes.fa"
    file2="results/${gene}/${gene}_protein_path_haplotypes.fa"

    # Check if both files exist
    if [ -f "$file1" ] && [ -f "$file2" ]; then
        echo "Both files exist for gene: $gene"
    else
        [ ! -f "$file1" ] && echo "Missing 1kG haplotypes ($file1) for gene: $gene"
        [ ! -f "$file2" ] && echo "Missing patho haplotypes ($file2) for gene: $gene"
        ((missing_count++))
    fi

    #./genename2foldedhaplo.sh "$gene"

done < "$CURRENT_DIR/$UNIPROTSIGNALOUT"

# Exit with error code 1 if any missing files were encountered
if [ "$missing_count" -gt 0 ]; then
    echo "Error: required file(s) missing." >&2
    exit 1
fi

while IFS=',' read -r gene np_id nm_id aa_length uniprot signal_start signal_end signal_length status || [ -n "$gene" ]; do
    [ "$gene" = "gene" ] && continue

    ./genename2foldedhaplo.sh "$gene"

done < "$CURRENT_DIR/$UNIPROTSIGNALOUT"
