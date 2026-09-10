#!/usr/bin/env bash

set -euo pipefail

# File definitions
CLINVAR="../data/clinvar_20260503.vcf.gz"
SMALLGENES="small-genes_mane-select_1-100.csv"
MANEGFF="../data/MANE.GRCh38.v1.3.ensembl_genomic.gff.gz"
OUTPUT="suitable_genes.txt"

while IFS= read -r gene || [ -n "$gene" ]; do
    # Clean whitespace and carriage returns (\r) from Windows/DOS line endings
    gene=$(echo "$gene" | tr -d '\r\n ')

    # Skip empty lines
    [ -z "$gene" ] && continue

    # Process your gene here
    echo "Processing gene: $gene"

done < "suitable_genes.txt"