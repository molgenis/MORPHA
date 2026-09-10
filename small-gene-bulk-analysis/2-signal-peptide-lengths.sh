#!/bin/bash

GENES="genes_with_path_var.txt"
MANE="short_MANE_proteins.csv"

OUTPUT="uniprot_signal_peptides.csv"
ERRORS="uniprot_signal_peptides_errors.tsv"

TMP_DIR=$(mktemp -d)
trap 'rm -rf "$TMP_DIR"' EXIT

# ------------------------------------------------------------
# Check dependencies
# ------------------------------------------------------------

for cmd in curl jq awk sed tr wc; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "ERROR: $cmd not found"
        exit 1
    fi
done

# ------------------------------------------------------------
# Output files
# ------------------------------------------------------------

printf "gene,NP_identifier,NM_identifier,AA_length,UniProt,Signal_start,Signal_end,Signal_length,status\n" \
    > "$OUTPUT"

printf "gene\tNP_identifier\tNM_identifier\tstatus\tdetails\n" \
    > "$ERRORS"

# ------------------------------------------------------------
# Build input table ONLY from genes_with_path_var.txt
# ------------------------------------------------------------

INPUT="$TMP_DIR/input.tsv"

while IFS= read -r raw_gene || [ -n "$raw_gene" ]; do

    gene=$(printf '%s' "$raw_gene" |
        tr -d '\r' |
        sed 's/^ *//;s/ *$//')

    # Skip blank lines
    [ -z "$gene" ] && continue

    # Find this gene in short_MANE_proteins.csv
    row=$(awk -F',' -v g="$gene" '
        NR > 1 && $1 == g {
            print $0
            exit
        }
    ' "$MANE")

    if [ -z "$row" ]; then

        printf "%s\t\t\tNO_MANE_MATCH\tGene not found in %s\n" \
            "$gene" "$MANE" >> "$ERRORS"

        printf "%s,,,,,,,,NO_MANE_MATCH\n" "$gene" >> "$OUTPUT"

        echo "$gene: NO MANE MATCH"
        continue
    fi

    np=$(printf '%s\n' "$row" | awk -F',' '{print $2}')
    nm=$(printf '%s\n' "$row" | awk -F',' '{print $3}')
    aa=$(printf '%s\n' "$row" | awk -F',' '{print $4}')

    printf "%s\t%s\t%s\t%s\n" \
        "$gene" "$np" "$nm" "$aa" >> "$INPUT"

done < "$GENES"

N_GENES=$(wc -l < "$INPUT" | tr -d ' ')

echo
echo "======================================"
echo "UniProt signal peptide analysis"
echo "======================================"
echo "Input gene list: $GENES"
echo "MANE file:       $MANE"
echo "Genes to query:  $N_GENES"
echo

# ------------------------------------------------------------
# Query each NP identifier
# ------------------------------------------------------------

while IFS="$(printf '\t')" read -r gene np nm aa; do

    echo "--------------------------------------"
    echo "$gene"
    echo "  NP: $np"
    echo "  NM: $nm"
    echo "  AA: $aa"

    # --------------------------------------------------------
    # Query UniProt using explicit RefSeq cross-reference
    # --------------------------------------------------------

    SEARCH_JSON="$TMP_DIR/${gene}_search.json"

    #QUERY="(organism_id:9606) AND (database:RefSeq) AND (xref:${np})"
    QUERY="(organism_id:9606) AND (reviewed:true) AND (database:RefSeq) AND (xref:${np})"

    if ! curl -fsS --retry 5 --retry-delay 2 \
        -G "https://rest.uniprot.org/uniprotkb/search" \
        --data-urlencode "query=$QUERY" \
        --data-urlencode "format=json" \
        --data-urlencode "fields=accession,id,gene_names,length" \
        --data-urlencode "size=10" \
        > "$SEARCH_JSON"; then

        printf "%s\t%s\t%s\tUNIPROT_QUERY_ERROR\t%s\n" \
            "$gene" "$np" "$nm" "$np" >> "$ERRORS"

        printf "%s,%s,%s,%s,,,,,UNIPROT_QUERY_ERROR\n" \
            "$gene" "$np" "$nm" "$aa" >> "$OUTPUT"

        echo "  ERROR: UniProt query failed"
        continue
    fi

    # --------------------------------------------------------
    # Number of UniProt matches
    # --------------------------------------------------------

    MATCHES=$(jq -r '.results | length' "$SEARCH_JSON")

    echo "  UniProt matches: $MATCHES"

    # --------------------------------------------------------
    # No match
    # --------------------------------------------------------

    if [ "$MATCHES" -eq 0 ]; then

        printf "%s\t%s\t%s\tNO_UNIPROT_MATCH\t%s\n" \
            "$gene" "$np" "$nm" "$np" >> "$ERRORS"

        printf "%s,%s,%s,%s,,,,,NO_UNIPROT_MATCH\n" \
            "$gene" "$np" "$nm" "$aa" >> "$OUTPUT"

        continue
    fi

    # --------------------------------------------------------
    # Multiple matches
    # --------------------------------------------------------

    if [ "$MATCHES" -gt 1 ]; then

        ACCS=$(jq -r '.results[].primaryAccession' "$SEARCH_JSON" |
            tr '\n' ';' |
            sed 's/;$//')

        printf "%s\t%s\t%s\tMULTIPLE_UNIPROT_MATCHES\t%s\n" \
            "$gene" "$np" "$nm" "$ACCS" >> "$ERRORS"

        printf "%s,%s,%s,%s,,,,,MULTIPLE_UNIPROT_MATCHES\n" \
            "$gene" "$np" "$nm" "$aa" >> "$OUTPUT"

        echo "  WARNING: multiple UniProt matches: $ACCS"
        continue
    fi

    # --------------------------------------------------------
    # Extract UniProt accession
    # --------------------------------------------------------

    uniprot=$(jq -r '.results[0].primaryAccession' "$SEARCH_JSON")

    echo "  UniProt: $uniprot"

    # --------------------------------------------------------
    # Retrieve complete UniProt JSON
    # --------------------------------------------------------

    ENTRY="$TMP_DIR/${gene}_entry.json"

    if ! curl -fsS --retry 5 --retry-delay 2 \
        "https://rest.uniprot.org/uniprotkb/${uniprot}.json" \
        > "$ENTRY"; then

        printf "%s\t%s\t%s\tUNIPROT_RETRIEVAL_ERROR\t%s\n" \
            "$gene" "$np" "$nm" "$uniprot" >> "$ERRORS"

        printf "%s,%s,%s,%s,%s,,,,UNIPROT_RETRIEVAL_ERROR\n" \
            "$gene" "$np" "$nm" "$aa" "$uniprot" >> "$OUTPUT"

        continue
    fi

    # --------------------------------------------------------
    # Extract Signal peptide feature
    # --------------------------------------------------------

    SIGNAL=$(jq -r '
        [
            .features[]?
            | select(.type == "Signal")
            | [
                .location.start.value,
                .location.end.value
              ]
            | @tsv
        ]
        | .[0] // empty
    ' "$ENTRY")

    # --------------------------------------------------------
    # No signal peptide
    # --------------------------------------------------------

    if [ -z "$SIGNAL" ]; then

        printf "%s,%s,%s,%s,%s,,,,NO_SIGNAL_PEPTIDE\n" \
            "$gene" \
            "$np" \
            "$nm" \
            "$aa" \
            "$uniprot" \
            >> "$OUTPUT"

        echo "  Signal peptide: none"
        continue
    fi

    # --------------------------------------------------------
    # Signal peptide coordinates
    # --------------------------------------------------------

    signal_start=$(printf '%s' "$SIGNAL" | awk '{print $1}')
    signal_end=$(printf '%s' "$SIGNAL" | awk '{print $2}')

    signal_length=$((signal_end - signal_start + 1))

    echo "  Signal peptide: ${signal_start}-${signal_end}"
    echo "  Signal length:  ${signal_length} aa"

    # --------------------------------------------------------
    # Write result
    # --------------------------------------------------------

    printf "%s,%s,%s,%s,%s,%s,%s,%s,OK\n" \
        "$gene" \
        "$np" \
        "$nm" \
        "$aa" \
        "$uniprot" \
        "$signal_start" \
        "$signal_end" \
        "$signal_length" \
        >> "$OUTPUT"

done < "$INPUT"

# ------------------------------------------------------------
# Final summary
# ------------------------------------------------------------

echo
echo "======================================"
echo "DONE"
echo "======================================"

echo "Input genes:"
grep -cve '^[[:space:]]*$' "$GENES"

echo "Output rows:"
tail -n +2 "$OUTPUT" | wc -l

echo
echo "Output:"
echo "  $OUTPUT"

echo
echo "Errors:"
echo "  $ERRORS"
