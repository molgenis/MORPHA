#!/bin/bash
set -euo pipefail

UNIPROTSIGNALOUT="uniprot_signal_peptides.csv"

while IFS=',' read -r gene np_id nm_id aa_length uniprot signal_start signal_end signal_length status || [ -n "$gene" ]; do
    [ "$gene" = "gene" ] || [ -z "$gene" ] && continue

    echo "=== Processing gene: $gene ==="

    # Path relative to script location (repo/small-gene-analysis/)
    GENE_DIR="../results/$gene"
    FOLDED_DIR="$GENE_DIR/folded"

    FASTA_FILES=(
        "$GENE_DIR/${gene}_protein_haplotypes.fa"
        "$GENE_DIR/${gene}_protein_path_haplotypes.fa"
    )

    VALID_INPUTS=true
    for fa_file in "${FASTA_FILES[@]}"; do
        if [ ! -f "$fa_file" ]; then
            echo "ERROR: FASTA file missing: $fa_file"
            VALID_INPUTS=false
        fi
    done

    if [ "$VALID_INPUTS" = false ]; then
        echo "Skipping $gene due to missing FASTA files."
        continue
    fi

    # parse FASTA headers to determine EXPECTED .sc files
    EXPECTED_SCORE_FILES=()
    for fa_file in "${FASTA_FILES[@]}"; do
        while read -r line; do
            header="${line#>}"
            header="${header#"${header%%[![:space:]]*}"}"
            header="${header%"${header##*[![:space:]]}"}"
            id_only="${header%%|*}"
            EXPECTED_SCORE_FILES+=("${id_only}_Repair.sc")
        done < <(grep '^>' "$fa_file")
    done

    # collect ACTUAL .sc files in $FOLDED_DIR
    ACTUAL_SCORE_FILES=()
    if [ -d "$FOLDED_DIR" ]; then
        while IFS= read -r -d '' sc_file; do
            ACTUAL_SCORE_FILES+=("$(basename "$sc_file")")
        done < <(find "$FOLDED_DIR" -maxdepth 1 -type f -name "*_Repair.sc" -print0)
    fi

    # validate 1:1 match
    IFS=$'\n' SORTED_EXPECTED=($(sort <<<"${EXPECTED_SCORE_FILES[*]}"))
    IFS=$'\n' SORTED_ACTUAL=($(sort <<<"${ACTUAL_SCORE_FILES[*]}"))
    unset IFS

    VALIDATION_PASSED=true
    if [ "${#SORTED_EXPECTED[@]}" -ne "${#SORTED_ACTUAL[@]}" ]; then
        echo "VALIDATION FAILED for $gene:"
        echo "  Expected ${#SORTED_EXPECTED[@]} scorefiles, but found ${#SORTED_ACTUAL[@]} in $FOLDED_DIR"
        VALIDATION_PASSED=false
    else
        for i in "${!SORTED_EXPECTED[@]}"; do
            if [ "${SORTED_EXPECTED[$i]}" != "${SORTED_ACTUAL[$i]}" ]; then
                VALIDATION_PASSED=false
                break
            fi
        done
    fi

    if [ "$VALIDATION_PASSED" = false ]; then
        echo "Mismatch detail for $gene:"
        echo "  Expected scorefiles from FASTAs: ${SORTED_EXPECTED[*]}"
        echo "  Actual scorefiles found in $FOLDED_DIR: ${SORTED_ACTUAL[*]}"
        echo "Quitting..."
        echo "--------------------------------------------------"
        exit 1
    fi

    echo "Validation PASSED: FASTA headers match existing .sc files (${#SORTED_EXPECTED[@]} files)."

    # extract minimum total_score per scorefile
    # output file target: repo/results/$gene/${gene}_minimum_per_haplotype.csv
    OUTPUT_CSV_NAME="${gene}_minimum_per_haplotype.csv"
    # also, extract all total_scores per scorefile in separate CSV
    ALL_DG_OUTPUT_CSV_NAME="${gene}_all_dg_per_haplotype.csv"

    # Initialize file
    > "$GENE_DIR/$OUTPUT_CSV_NAME"
    > "$GENE_DIR/$ALL_DG_OUTPUT_CSV_NAME"

    # change working dir to folded
    pushd "$FOLDED_DIR" > /dev/null || continue

    for sc_file in "${SORTED_EXPECTED[@]}"; do
        haplotype_id="${sc_file%_Repair.sc}"

        # Get minimum total_score from column 2
        min_score=$(grep '^SCORE:' "$sc_file" | grep -v 'total_score' | awk '{print $2}' | sort -n | head -n 1)

        if [ -n "$min_score" ]; then
            # write to parent directory (results/$gene/)
            echo "${haplotype_id},${min_score}" >> "../$OUTPUT_CSV_NAME"

            # also write all DG scores per haplotype to CSV, i.e. the outcomes of all N trajectories resulting from nstruct=N
			scores=$(grep '^SCORE:' "$sc_file" | grep -v 'total_score' | awk '{print $2}' | paste -sd,)
            echo "${haplotype_id},${scores}" >> "../$ALL_DG_OUTPUT_CSV_NAME"
        else
            echo "ERROR: Could not extract total_score from $sc_file"
            exit 1
        fi
    done

    # change working dir back to current
    popd > /dev/null
    echo "Finished generating CSV for $gene at repo/results/$gene/$OUTPUT_CSV_NAME."
    echo "Finished generating CSV for $gene at repo/results/$gene/$ALL_DG_OUTPUT_CSV_NAME."
    echo "--------------------------------------------------"

done < "$UNIPROTSIGNALOUT"
