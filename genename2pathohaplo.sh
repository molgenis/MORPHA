#!/bin/bash


[ -z "$1" ] && { echo "Missing argument: gene name"; exit 1; }

GENE=$1

DATADIR=data
TMPDIR=tmp
GENEDIR="${TMPDIR}/${GENE}"
SAMPLENAMES="1kGP_sample_names.txt"
CLINVARGENE="${GENEDIR}/ClinVarPathogenic"
CLINVARGENELOCAL="${CLINVARGENE}.local.vcf.gz"
GENESAMPLESHAPLODIR="${GENEDIR}/samples"
GENESAMPLESPATHODIR="${GENEDIR}/samples_wpatho"
mkdir -p "${GENESAMPLESPATHODIR}"

echo "Starting script with gene name = ${GENE}"


# create consensus based on previously extracted samples (separate script?)



bcftools query -f '%CHROM:%POS:%ALT\n' "$CLINVARGENELOCAL" | while read -r variant; do

    # DEBUG
    variant="chr16:1379:TAC"

    # Isolate the position, chromosome and alternative allele
    CHROM=$(echo "$variant" | cut -d':' -f1)
    POS=$(echo "$variant" | cut -d':' -f2)
    ALT=$(echo "$variant" | cut -d':' -f3)

    # Create a unique file suffix including the ALT base to prevent overwriting
    CHR_POS_ALT="${CHROM}_${POS}_${ALT}"
    echo "Processing variant at $CHROM:$POS with ALT allele $ALT..."

    # Make a temporary single-variant VCF on the fly plus index
    SINGLE_VAR_TMP_VCF="${GENEDIR}/tmp_vcf_for_${CHR_POS_ALT}.vcf.gz"
    bcftools view -O z -i "POS==${POS} && ALT=='${ALT}'" "${CLINVARGENELOCAL}" > "${SINGLE_VAR_TMP_VCF}"
    bcftools index "${SINGLE_VAR_TMP_VCF}"
    if [ ! -s "${SINGLE_VAR_TMP_VCF}" ]; then
      echo "File ${SINGLE_VAR_TMP_VCF} missing or empty."
      exit 1
    fi
    echo "Single-variant VCF created at ${SINGLE_VAR_TMP_VCF}."

    # now loop over all samples
    while read -r SAMPLE; do
        [ -z "$SAMPLE" ] && continue # skip empty lines
        for H in 1 2; do
            PREFIX_HAP="${GENESAMPLESHAPLODIR}/${SAMPLE}.h${H}"
            PREFIX_PTH="${GENESAMPLESPATHODIR}/${SAMPLE}.h${H}"

            echo "Applying to sample: ${PREFIX_HAP}"
            SAMPLEFA="${PREFIX_HAP}.fa"
            if [ ! -s "${SAMPLEFA}" ]; then
              echo "File ${SAMPLEFA} missing or empty for ${SAMPLE} h${H}."
              exit 1
            fi

            # Temporarily simplify your FASTA header so it says ">chr16" instead of ">chr16:3242..."
            # to enable bcftools consensus matching of VCF and FASTA
            SAMPLEFAHEADER=$(grep ">" "${SAMPLEFA}" | head -n 1)
            echo "Extracted original FASTA header: ${SAMPLEFAHEADER}"
            sed "s/>.*/>${CHROM}/" "${SAMPLEFA}" > "${PREFIX_PTH}.renameheader.fa"
            if [ ! -s "${PREFIX_PTH}.renameheader.fa" ]; then
              echo "File ${PREFIX_PTH}.renameheader.fa missing or empty for ${SAMPLE} h${H}."
              exit 1
            fi
            echo "Created temporary header rename FASTA at ${PREFIX_PTH}.renameheader.fa"

            # 4. Run consensus
            bcftools consensus \
                -f "${PREFIX_PTH}.renameheader.fa" \
                -c "${PREFIX_PTH}.${CHR_POS_ALT}.chain" \
                "$SINGLE_VAR_TMP_VCF" \
                > "${PREFIX_PTH}.${CHR_POS_ALT}.fa"
            if [ ! -s "${PREFIX_PTH}.${CHR_POS_ALT}.fa" ]; then
              echo "File ${PREFIX_PTH}.${CHR_POS_ALT}.fa missing or empty for ${SAMPLE} h${H}."
              exit 1
            fi
            echo "Created FASTA file for ${SAMPLE} h${H} at ${PREFIX_PTH}.${CHR_POS_ALT}.fa"
            if [ ! -s "${PREFIX_PTH}.${CHR_POS_ALT}.chain" ]; then
              echo "File ${PREFIX_PTH}.${CHR_POS_ALT}.chain missing or empty for ${SAMPLE} h${H}."
              exit 1
            fi
            echo "Created chain file for ${SAMPLE} h${H} at ${PREFIX_PTH}.${CHR_POS_ALT}.chain"


            # 2. FIX THE HEADER: Replace ">chr16" back with the original coordinate header
            # '1s/.*/.../' tells sed to only modify the first line of the file
            sed -i "" "1s/.*/${SAMPLEFAHEADER}/" "${PREFIX_PTH}.${CHR_POS_ALT}.fa"
            echo "Restored header in ${PREFIX_PTH}.${CHR_POS_ALT}.fa"

            # 5. Clean up temporary files
            #rm "$tmp_vcf" "${tmp_vcf}.csi" tmp_temp_ref.fa


            HAPLIFTEDSAMPLEGFF="${PREFIX_HAP}.lifted.local.MANE.GRCh38.v1.3.ensembl_genomic.gff"
            CHROMRENAME_HAPLIFTEDSAMPLEGFF="${PREFIX_PTH}.lifted.local.chromrename.MANE.GRCh38.v1.3.ensembl_genomic.gff"

            if [ ! -s "${HAPLIFTEDSAMPLEGFF}" ]; then
              echo "File ${HAPLIFTEDSAMPLEGFF}or empty for ${SAMPLE} h${H}."
              exit 1
            fi
            # create temporary GFF where the chromosome name column is simplified to match the chain file (changes "chr16:3242027-3256633" to "chr16")
            awk -F'\t' -v chrom="$CHROM" 'BEGIN{OFS="\t"} {!/^#/ && $1=chrom; print}' \
                "${HAPLIFTEDSAMPLEGFF}" \
                > "${CHROMRENAME_HAPLIFTEDSAMPLEGFF}"
            if [ ! -s "${CHROMRENAME_HAPLIFTEDSAMPLEGFF}" ]; then
              echo "File ${CHROMRENAME_HAPLIFTEDSAMPLEGFF} missing or empty for ${SAMPLE} h${H}."
              exit 1
            fi
            echo "Create temporary GFF with chromosome renamed to ${CHROM}"
            # ALSO RENAME THE FA HEADER...
            sed '1s/>.*/>chr16/' "${PREFIX_PTH}.${CHR_POS_ALT}.fa" > "${PREFIX_PTH}.${CHR_POS_ALT}.renamed.fa"

            # then perform liftover to pathogenic variant coordinates
            PATHOLIFTEDSAMPLEGFF="${PREFIX_PTH}.patholifted.local.MANE.GRCh38.v1.3.ensembl_genomic.gff"
            echo "Running CrossMap: CrossMap gff ${PREFIX_PTH}.${CHR_POS_ALT}.chain ${CHROMRENAME_HAPLIFTEDSAMPLEGFF} ${PATHOLIFTEDSAMPLEGFF}"
            CrossMap gff \
              "${PREFIX_PTH}.${CHR_POS_ALT}.chain" \
              "${CHROMRENAME_HAPLIFTEDSAMPLEGFF}" \
              "${PATHOLIFTEDSAMPLEGFF}"

            # translate to CDS and protein:
            echo "Now translating to protein"
            echo "gffread ${PATHOLIFTEDSAMPLEGFF} -g ${PREFIX_PTH}.${CHR_POS_ALT}.renamed.fa -x ${PREFIX_PTH}.${CHR_POS_ALT}.cds.fa -y ${PREFIX_PTH}.${CHR_POS_ALT}.protein.fa -S"
            gffread "${PATHOLIFTEDSAMPLEGFF}" -g "${PREFIX_PTH}.${CHR_POS_ALT}.renamed.fa" -x "${PREFIX_PTH}.${CHR_POS_ALT}.cds.fa" -y "${PREFIX_PTH}.${CHR_POS_ALT}.protein.fa" -S



        done
    done < "${DATADIR}/${SAMPLENAMES}"


done

# liftover MANE coordinates:



# % diff HG00620.h1.chr16_14348.protein.fa tmp/MEFV/samples/HG00620.h1.protein.fa
# 11c11
# < DGPQRFDSCIIVLGSPSFLSGRRYWEVEVGDKTAWILGACKTSISRKGNMTLSPENGYWVVIMVKENEYQ
# ---
# > DGPQRFDSCIIVLGSPSFLSGRRYWEVEVGDKTAWILGACKTSISRKGNMTLSPENGYWVVIMMKENEYQ