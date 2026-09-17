#!/usr/bin/env bash

set -euo pipefail

UNIPROTSIGNALOUT="uniprot_signal_peptides.csv"


# Remember where we started
CURRENT_DIR="$(pwd)"

# Parent directory where genename2uniqproteinhaplo.sh expects to run
cd ..

while IFS=',' read -r gene np_id nm_id aa_length uniprot signal_start signal_end signal_length status || [ -n "$gene" ]; do

    # Skip header
    [ "$gene" = "gene" ] && continue

    # Paths to the two resulting FASTA files
    NORMAL_HAPLO="results/$gene/${gene}_protein_haplotypes.fa"
    PATH_HAPLO="results/$gene/${gene}_protein_path_haplotypes.fa"

    # Skip gene if both already exist
    if [ -s "$NORMAL_HAPLO" ] && [ -s "$PATH_HAPLO" ]; then
        echo "=== $gene ==="
        echo "Both haplotype FASTA files already exist."
        echo "Skipping gene."
        continue
    fi

    # Default: no signal peptide
    cutoff=0

    # Use end of signal peptide as the front cutoff
    if [ -n "$signal_end" ]; then
        cutoff="$signal_end"
    fi

    echo "=== $gene ==="
    echo "UniProt: $uniprot"
    echo "Signal peptide end: ${signal_end:-none}"
    echo "Calling: ../genename2uniqproteinhaplo.sh $gene $cutoff"

    ./genename2uniqproteinhaplo.sh "$gene" "$np_id" "$cutoff"
    ./genename2pathohaplo.sh "$gene" "$np_id" "$CURRENT_DIR/genes/$gene/variants.vcf.gz" "$cutoff"
    rm -rf "tmp/$gene" || true # don't stop when file locks have not been fully released

done < "$CURRENT_DIR/$UNIPROTSIGNALOUT"
#done <<< "APOA2,NP_001634.1,NM_001643.2,100,P02652,1,18,18,OK"
