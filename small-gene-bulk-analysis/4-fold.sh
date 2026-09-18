#!/usr/bin/env bash

set -euo pipefail

# Remember where we started
CURRENT_DIR="$(pwd)"

# Parent directory where genename2uniqproteinhaplo.sh expects to run
cd ..

while IFS=',' read -r gene np_id nm_id aa_length uniprot signal_start signal_end signal_length status || [ -n "$gene" ]; do

    # Skip header
    [ "$gene" = "gene" ] && continue

    ./genename2foldedhaplo.sh "$gene"



done < "$CURRENT_DIR/$UNIPROTSIGNALOUT"
