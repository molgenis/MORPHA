#!/bin/bash
#SBATCH --gres=gpu:rtx6000:1
#SBATCH --job-name=esmfold2_mefv_asc_foldjobs
#SBATCH --output=esmfold2_mefv_asc_foldjobs.out
#SBATCH --error=esmfold2_mefv_asc_foldjobs.err
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --nodes=1
#SBATCH --export=NONE


[ -z "$1" ] && { echo "Missing argument: gene name"; exit 1; }

GENE=$1

echo "Starting script with gene name = ${GENE}"

RESULTDIR="results/${GENE}"
NORMALHAPLO="${RESULTDIR}/${GENE}_protein_haplotypes.fa"
PATHHAPLO="${RESULTDIR}/${GENE}_protein_path_haplotypes.fa"
MMCIFDIR="${RESULTDIR}/${GENE}/folded"
mkdir -p "${MMCIFDIR}"

# Iterate through both of your target FASTA files
for FASTA_FILE in "${NORMALHAPLO}" "${PATHHAPLO}"; do

    # Safety check: make sure the file actually exists and contains data
    if [ ! -s "${FASTA_FILE}" ]; then
        echo "Warning: ${FASTA_FILE} is missing or empty. Skipping to next file..."
        continue
    fi

    echo "================================================================="
    echo "Processing Input File: $(basename "${FASTA_FILE}")"
    echo "================================================================="

    awk '
    # Split on > if it is at the start of the file (^) OR preceded by a newline (\n)
    BEGIN { RS = "(\n|^)>" }
    NR >= 1 {
        # Clean any Windows carriage returns
        gsub(/\r/, "", $0)

        # Separate the header line from the sequence block
        nl_pos = index($0, "\n")
        if (nl_pos == 0) next
        full_header = substr($0, 1, nl_pos - 1)
        seq = substr($0, nl_pos + 1)

        # 2. Extract clean name (works for both normal and path headers)
        pipe_pos = index(full_header, "|")
        if (pipe_pos > 0) {
            clean_name = substr(full_header, 1, pipe_pos - 1)
        } else {
            clean_name = full_header
        }

        # 3. Strip line breaks and truncate cleanly at stop codon (*)
        gsub(/\n/, "", seq)
        stop_pos = index(seq, "*")
        if (stop_pos > 0) {
            seq = substr(seq, 1, stop_pos - 1)
        }

        # 4. Stream to the while loop
        if (length(clean_name) > 0 && length(seq) > 0) {
            print clean_name " " seq
        }
    }
    ' "${FASTA_FILE}" | while read -r VARIANT_NAME AMINO_ACIDS; do

        [ -z "$VARIANT_NAME" ] && continue

        echo "  Processing variant: ${VARIANT_NAME} (${#AMINO_ACIDS} AA)..."

        MMCIF_OUTPUT="${MMCIFDIR}/${VARIANT_NAME}.mmcif"
        MMCIF_TMP="${MMCIF_OUTPUT}.tmp"

        # Checkpoint: Skip if the finalized MMCIF already exists
        if [ -s "${MMCIF_OUTPUT}" ]; then
            echo "     -> MMCIF exists. Skipping..."
            continue
        fi

        echo "     -> Requesting fold..."

        # Clean up any broken leftover temp files from a previous aborted run
        rm -f "${MMCIF_TMP}"

        # 1. Download directly to the temporary file
        curl -X POST http://localhost:8844/fold \
            -H "Content-Type: application/json" \
            -d "{ \"sequences\": [{\"id\": \"A\", \"sequence\": \"${AMINO_ACIDS}\", \"input_type\": \"protein\"},
                  {\"id\": \"B\", \"sequence\": \"MGRARDAILDALENLTAEELKKFKLKLLSVPLREGYGRIPRGALLSMDALDLTDKLVSFYLETYGAELTANVLRDMGLQEMAGQLQAATHQGSGAAPAGIQAPPQSAAKPGLHFIDQHRAALIARVTNVEWLLDALYGKVLTDEQYQAVRAEPTNPSKMRKLFSFTPAWNWTCKDLLLQALRESQSYLVEDLERS\", \"input_type\": \"protein\"}] }" \
            -o "${MMCIF_TMP}"

          
        # 2. Only promote to a real .MMCIF if the download completed successfully
        if [ -s "${MMCIF_TMP}" ]; then
            mv "${MMCIF_TMP}" "${MMCIF_OUTPUT}"
        else
            echo "     [WARNING] Fold failed or server returned empty data."
            rm -f "${MMCIF_TMP}"
        fi

    done
done

echo "Execution complete. All normal and pathogenic variants folded."


