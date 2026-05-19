#!/bin/bash


[ -z "$1" ] && { echo "Missing argument: gene name"; exit 1; }

GENE=$1

DATADIR=data
TMPDIR=tmp
GENEDIR="${TMPDIR}/${GENE}"

GRCh38REFFA=GRCh38_full_analysis_set_plus_decoy_hla.fa
GRCh38REFFAINDEX=GRCh38_full_analysis_set_plus_decoy_hla.fa.fai
GRCh38MANEGFF=MANE.GRCh38.v1.3.ensembl_genomic.gff
GRCh38MANEGFFGZ="${GRCh38MANEGFF}.gz"

echo "Starting script with gene name = ${GENE}"


# create consensus based on previously extracted samples (separate script?)

exit 0

# Work directly with your original VCF that has clean 'chr16' names
ORIGINAL_VCF="tmp/MEFV/ClinVarPathogenic.local.vcf.gz"
FASTA_REF="tmp/MEFV/samples/HG00620.h1.fa"

# 1. Grab the exact original header to reuse later
original_header=$(grep ">" "$FASTA_REF" | head -n 1)

bcftools query -f '%CHROM:%POS:%ALT\n' "$ORIGINAL_VCF" | while read -r variant; do

    # 1. Isolate the position, chromosome and alternative allele
    chrom=$(echo "$variant" | cut -d':' -f1)
    pos=$(echo "$variant" | cut -d':' -f2)
    alt=$(echo "$variant" | cut -d':' -f3)

    # Create a unique file suffix including the ALT base to prevent overwriting
    out_name="${chrom}_${pos}_${alt}"
    echo "Processing variant at $chrom:$pos with ALT allele $alt..."

    # 2. Make a temporary VCF on the fly with a physical file + index
    tmp_vcf="tmp_single_${out_name}.vcf.gz"
    bcftools view -O z -i "POS==$pos && ALT=='$alt'" "$ORIGINAL_VCF" > "$tmp_vcf"
    bcftools index "$tmp_vcf"

    # 3. Temporarily simplify your FASTA header so it says ">chr16" instead of ">chr16:3242..."
    # This tricks bcftools consensus into matching the VCF and FASTA perfectly.
    sed "s/>.*/>${chrom}/" tmp/MEFV/samples/HG00620.h1.fa > tmp_temp_ref.fa

    # 4. Run consensus
    bcftools consensus \
        -f tmp_temp_ref.fa \
        -c "HG00620.h1.${out_name}.chain" \
        "$tmp_vcf" \
        > "HG00620.h1.${out_name}.fa"

    # 2. FIX THE HEADER: Replace ">chr16" back with the original coordinate header
    # '1s/.*/.../' tells sed to only modify the first line of the file
    sed -i "" "1s/.*/${original_header}/" "HG00620.h1.${out_name}.fa"

    # 5. Clean up temporary files
    rm "$tmp_vcf" "${tmp_vcf}.csi" tmp_temp_ref.fa

done

# liftover MANE coordinates:

# first create a temporary GFF where the chromosome name column is simplified to match the chain file (changes "chr16:3242027-3256633" to "chr16")
awk -F'\t' 'BEGIN{OFS="\t"} {!/^#/ && $1="chr16"; print}' \
    "tmp/MEFV/samples/HG00620.h1.lifted.local.MANE.GRCh38.v1.3.ensembl_genomic.gff" \
    > tmp_normalized_source.gff
# then perform liftover
CrossMap gff \
    "HG00620.h1.${out_name}.chain" \
    "tmp_normalized_source.gff" \
    "HG00620.h1.${out_name}.pathovar.lifted.gff"

# translate to CDS and protein:
gffread "HG00620.h1.${out_name}.pathovar.lifted.gff" -g "HG00620.h1.chr16_14348.fa" -x "HG00620.h1.chr16_14348.cds.fa" -y "HG00620.h1.chr16_14348.protein.fa" -S


# % diff HG00620.h1.chr16_14348.protein.fa tmp/MEFV/samples/HG00620.h1.protein.fa
# 11c11
# < DGPQRFDSCIIVLGSPSFLSGRRYWEVEVGDKTAWILGACKTSISRKGNMTLSPENGYWVVIMVKENEYQ
# ---
# > DGPQRFDSCIIVLGSPSFLSGRRYWEVEVGDKTAWILGACKTSISRKGNMTLSPENGYWVVIMMKENEYQ