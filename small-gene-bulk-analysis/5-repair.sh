#!/bin/bash
set -euo pipefail

# Remember where we started
CURRENT_DIR="$(pwd)"

# Configuration
FOLDX_BIN="/mnt/c/ProteinStuff/wsl/foldx_20270131"
UNIPROTSIGNALOUT="uniprot_signal_peptides.csv"

# step into arent directory
cd ..

while IFS=',' read -r gene np_id nm_id aa_length uniprot signal_start signal_end signal_length status || [ -n "$gene" ]; do
    [ "$gene" = "gene" ] || [ -z "$gene" ] && continue

    echo "=== Processing Gene: $gene ==="
    
    GENE_DIR="results/$gene"
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

	# 1. Parse FASTA headers (Pure Bash - no pipe or xargs limit)
    EXPECTED_PDBS=()
    for fa_file in "${FASTA_FILES[@]}"; do
        while read -r line; do
            # Strip leading '>'
            header="${line#>}"
            
            # Trim leading/trailing whitespace
            header="${header#"${header%%[![:space:]]*}"}"
            header="${header%"${header##*[![:space:]]}"}"
            
            # Strip everything from the first pipe '|' onward
            id_only="${header%%|*}"
            
            EXPECTED_PDBS+=("${id_only}.pdb")
        done < <(grep '^>' "$fa_file")
    done

    # 2. Collect actual PDB files (EXCLUDING *_Repair.pdb)
    ACTUAL_PDBS=()
    if [ -d "$FOLDED_DIR" ]; then
        while IFS= read -r -d '' pdb_file; do
            ACTUAL_PDBS+=("$(basename "$pdb_file")")
        done < <(find "$FOLDED_DIR" -maxdepth 1 -type f -name "*.pdb" ! -name "*_Repair.pdb" -print0)
    fi

    # 3. Validate
    IFS=$'\n' SORTED_EXPECTED=($(sort <<<"${EXPECTED_PDBS[*]}"))
    IFS=$'\n' SORTED_ACTUAL=($(sort <<<"${ACTUAL_PDBS[*]}"))
    unset IFS

    VALIDATION_PASSED=true
    if [ "${#SORTED_EXPECTED[@]}" -ne "${#SORTED_ACTUAL[@]}" ]; then
        echo "VALIDATION FAILED for $gene:"
        echo "  Expected ${#SORTED_EXPECTED[@]} PDBs, but found ${#SORTED_ACTUAL[@]} in $FOLDED_DIR"
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
        echo "  Expected from FASTAs: ${SORTED_EXPECTED[*]}"
        echo "  Found in $FOLDED_DIR: ${SORTED_ACTUAL[*]}"
        echo "Skipping FoldX for $gene."
        echo "--------------------------------------------------"
        continue
    fi

    echo "Validation PASSED: FASTA headers match existing PDB files (${#SORTED_EXPECTED[@]} files)."

    # 4. Execute FoldX Repair (Skip already repaired PDBs)
    pushd "$FOLDED_DIR" > /dev/null || continue
    
    for pdb_name in "${SORTED_EXPECTED[@]}"; do
        repair_output="${pdb_name%.pdb}_Repair.pdb"
        
        # Check if output already exists
        if [ -f "$repair_output" ]; then
            echo "Skipping $pdb_name (repaired file $repair_output already exists)."
            continue
        fi

        echo "Running FoldX Repair on $pdb_name..."
        "$FOLDX_BIN" --command=RepairPDB --pdb="$pdb_name"
    done

    popd > /dev/null
    echo "Finished $gene."
    echo "--------------------------------------------------"

done < "$CURRENT_DIR/$UNIPROTSIGNALOUT"
