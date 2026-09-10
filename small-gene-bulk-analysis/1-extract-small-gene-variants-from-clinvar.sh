#!/usr/bin/env bash

set -euo pipefail

# File definitions
CLINVAR="../data/clinvar_20260503.vcf.gz"
SMALLGENES="short_MANE_proteins.csv"
MANEGFF="../data/MANE.GRCh38.v1.3.ensembl_genomic.gff.gz"
OUTPUT="genes_with_path_var.txt"

echo "=== Starting Variant Extraction Pipeline ==="

# detect VCF chromosome naming style
sample_vcf_chrom=$( (bcftools view -H "${CLINVAR}" 2>/dev/null || true) | head -n 1 | cut -f1)
echo "Chromosome naming style of ClinVar: $sample_vcf_chrom"

# Clear output file
> "$OUTPUT"

total_genes=$(grep -c '.' "$SMALLGENES" || true)
echo "Loaded $total_genes gene entries from $SMALLGENES"

count=0
matched_genes=0

# either Pathogenic, Likely_pathogenic, Pathogenic/Likely_pathogenic in combination with 1 star or better by ignoring anything with 'no_assertion'
FILTER_EXPR='(CLNSIG ~ "^Pathogenic$" || CLNSIG ~ "^Likely_pathogenic$" || CLNSIG ~ "^Pathogenic/Likely_pathogenic$") && (CLNREVSTAT !~ "no_assertion")'

while IFS=, read -r raw_gene _ raw_tx _ || [ -n "$raw_gene" ]; do
    gene=$(echo "$raw_gene" | tr -d '\r\n ')
    tx_id=$(echo "$raw_tx" | tr -d '\r\n ')

    count=$((count + 1))
    echo -ne "Processing [$count/$total_genes]: $gene ($tx_id)... "

    # Find matching transcript coordinates safely ignoring SIGPIPE from gzip/awk exit
    region_info=$( (gzip -dc "$MANEGFF" 2>/dev/null || true) | awk -v tx="$tx_id" '
    $3 == "transcript" && $0 ~ tx {
        print $1, $4, $5;
        exit
    }')

    if [ -z "$region_info" ]; then
        echo "Region not found in GFF (Skipping)"
        continue
    fi

    gff_chrom=$(echo "$region_info" | awk '{print $1}')
    start=$(echo "$region_info" | awk '{print $2}')
    end=$(echo "$region_info" | awk '{print $3}')

    echo "Region found: $gff_chrom:$start-$end"


    # Adjust chromosome naming to match VCF style
    if [[ "$sample_vcf_chrom" != chr* ]] && [[ "$gff_chrom" == chr* ]]; then
        target_chrom="${gff_chrom#chr}"
    elif [[ "$sample_vcf_chrom" == chr* ]] && [[ "$gff_chrom" != chr* ]]; then
        target_chrom="chr${gff_chrom}"
    else
        target_chrom="$gff_chrom"
    fi

    region="${target_chrom}:${start}-${end}"

    # Extract ONLY variants that match our filter
    vars_output=$(bcftools query -i "$FILTER_EXPR" -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' -r "$region" "${CLINVAR}" 2>/dev/null || true)

    if [ -n "$vars_output" ]; then
        mkdir -p genes/$gene

        echo "##fileformat=VCFv4.1" > "genes/$gene/variants.vcf"
        echo "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO" >> "genes/$gene/variants.vcf"
        # 1	1392	c.1226A>G	T	C	.	.	.
        formatted_vars=""
        var_count=0
        
        while IFS= read -r line; do
            if [ -z "$line" ]; then
                continue
            fi
            var_count=$((var_count + 1))
            printf "%s\t.\t.\t.\n" "$line" >> "genes/$gene/variants.vcf"
        done <<< "$vars_output"

        echo $gene >> $OUTPUT

        bgzip -c "genes/$gene/variants.vcf" > "genes/$gene/variants.vcf.gz"
        tabix -p vcf "genes/$gene/variants.vcf.gz"
        rm "genes/$gene/variants.vcf"

        matched_genes=$((matched_genes + 1))
        echo "Found $var_count Pathogenic variant(s) in $region"
    else
        echo "0 Pathogenic variants in $region"
    fi
done < "$SMALLGENES"
#done <<< "RPS29,NM_001032.5"


echo "--------------------------------------------------"
echo "=== Complete! ==="
echo "Processed: $count genes"
echo "Saved pathogenic variants for: $matched_genes genes -> $OUTPUT"