#!/bin/bash


[ -z "$1" ] && { echo "Missing argument: gene name"; exit 1; }

GENE=$1

echo "Starting script with gene name = ${GENE}"

DATADIR=data
TMPDIR=tmp
GENEDIR="${TMPDIR}/${GENE}"
SAMPLENAMES="1kGP_sample_names.txt"
CLINVARGENE="${GENEDIR}/ClinVarPathogenic"
CLINVARGENELOCAL="${CLINVARGENE}.local.vcf.gz"
GENESAMPLESHAPLODIR="${GENEDIR}/samples"
GENESAMPLESPATHODIR="${GENEDIR}/samples_patho"
mkdir -p "${GENESAMPLESPATHODIR}"

bcftools query -f '%CHROM:%POS:%ALT\n' "$CLINVARGENELOCAL" | while read -r variant; do

    # DEBUG
    variant="chr16:1379:TAC"

    # Isolate the position, chromosome and alternative allele
    # NOTE: we assume 1 alt per line, check if this is really true!!
    CHROM=$(echo "$variant" | cut -d':' -f1)
    POS=$(echo "$variant" | cut -d':' -f2)
    ALT=$(echo "$variant" | cut -d':' -f3)
    CHR_POS_ALT="${CHROM}_${POS}_${ALT}"

    echo "Processing variant at $CHROM:$POS with ALT allele $ALT..."

    # Make a temporary single-variant VCF on the fly plus index
    SINGLE_VAR_TMP_VCF="${GENEDIR}/tmp_singlevar_${CHR_POS_ALT}.vcf.gz"
    if [ ! -s "${SINGLE_VAR_TMP_VCF}" ]; then
      echo "Creating ${SINGLE_VAR_TMP_VCF}..."
      bcftools view -O z -i "POS==${POS} && ALT=='${ALT}'" "${CLINVARGENELOCAL}" > "${SINGLE_VAR_TMP_VCF}"
      bcftools index "${SINGLE_VAR_TMP_VCF}"
    else
      echo "Skipping creation of ${SINGLE_VAR_TMP_VCF}..."
    fi
    if [ ! -s "${SINGLE_VAR_TMP_VCF}" ]; then
      echo "File ${SINGLE_VAR_TMP_VCF} missing or empty."
      exit 1
    fi

    # now loop over all samples
    while read -r SAMPLE; do
        [ -z "$SAMPLE" ] && continue # skip empty lines
        for H in 1 2; do

            echo "Going to work on sample: ${SAMPLE}.h${H}"
            PREFIX_HAP="${GENESAMPLESHAPLODIR}/${SAMPLE}.h${H}"
            PREFIX_PTH="${GENESAMPLESPATHODIR}/${SAMPLE}.h${H}"
            ORIGINALSAMPLEFA="${PREFIX_HAP}.fa"

            if [ ! -s "${ORIGINALSAMPLEFA}" ]; then
              echo "File ${ORIGINALSAMPLEFA} missing or empty for ${SAMPLE} h${H}."
              exit 1
            fi

            # Temporarily simplify your FASTA header so it says ">chr16" instead of ">chr16:3242..."
            # to enable bcftools consensus matching of VCF and FASTA
            ORIGINALSAMPLEFAHEADER=$(grep ">" "${ORIGINALSAMPLEFA}" | head -n 1)
            echo "Extracted original FASTA header: ${ORIGINALSAMPLEFAHEADER}"
            RENAMEDHEADERSAMPLEFA="${PREFIX_PTH}.hdr.fa"
            if [ ! -s "${RENAMEDHEADERSAMPLEFA}" ]; then
              echo "Creating ${RENAMEDHEADERSAMPLEFA}..."
              #sed "s/>.*/>${CHROM}/" "${ORIGINALSAMPLEFAHEADER}" > "${RENAMEDHEADERSAMPLEFA}"
              sed "1s/>.*/>${CHROM}/" "${ORIGINALSAMPLEFA}" > "${RENAMEDHEADERSAMPLEFA}"
            else
              echo "Skipping creation of ${RENAMEDHEADERSAMPLEFA}..."
            fi
            if [ ! -s "${RENAMEDHEADERSAMPLEFA}" ]; then
              echo "File ${RENAMEDHEADERSAMPLEFA} missing or empty for ${SAMPLE} h${H}."
              exit 1
            fi

            # 4. Run consensus
            SAMPLEPATHCHAIN="${PREFIX_PTH}.${CHR_POS_ALT}.chain"
            SAMPLEPATHFA="${PREFIX_PTH}.${CHR_POS_ALT}.fa"

            if [[ ! -s "${SAMPLEPATHCHAIN}" || ! -s "${SAMPLEPATHFA}" ]]; then
                echo "Creating ${SAMPLEPATHCHAIN} and ${SAMPLEPATHFA} sequence consensus..."
                bcftools consensus \
                    -f "${RENAMEDHEADERSAMPLEFA}" \
                    -c "${SAMPLEPATHCHAIN}" \
                    "${SINGLE_VAR_TMP_VCF}" \
                    > "${SAMPLEPATHFA}"
            else
              echo "${SAMPLEPATHCHAIN} and ${SAMPLEPATHFA} present, skipping consensus..."
            fi
            if [[ ! -s "${SAMPLEPATHCHAIN}" || ! -s "${SAMPLEPATHFA}" ]]; then
               echo "${SAMPLEPATHCHAIN} and/or ${SAMPLEPATHFA} are missing."
               exit 1
            fi

            # 2. FIX THE HEADER: Replace ">chr16" back with the original coordinate header
            # '1s/.*/.../' tells sed to only modify the first line of the file
            sed -i "" "1s/.*/${ORIGINALSAMPLEFAHEADER}/" "${SAMPLEPATHFA}"
            echo "Restored original FA header in ${SAMPLEPATHFA}"

            # 5. Clean up temporary files
            #rm "$tmp_vcf" "${tmp_vcf}.csi" tmp_temp_ref.fa

            HAPLIFTEDSAMPLEGFF="${PREFIX_HAP}.lifted.local.MANE.GRCh38.v1.3.ensembl_genomic.gff"
            CHROMRENAME_HAPLIFTEDSAMPLEGFF="${PREFIX_PTH}.lifted.local.chromrename.MANE.GRCh38.v1.3.ensembl_genomic.gff"

            if [ ! -s "${HAPLIFTEDSAMPLEGFF}" ]; then
              echo "File ${HAPLIFTEDSAMPLEGFF} missing or empty for ${SAMPLE} h${H}."
              exit 1
            fi
            # create temporary GFF where the chromosome name column is simplified to match the chain file (changes "chr16:3242027-3256633" to "chr16")
            if [ ! -s "${CHROMRENAME_HAPLIFTEDSAMPLEGFF}" ]; then
              echo "Creating ${CHROMRENAME_HAPLIFTEDSAMPLEGFF}..."
              awk -F'\t' -v chrom="$CHROM" 'BEGIN{OFS="\t"} {!/^#/ && $1=chrom; print}' \
                  "${HAPLIFTEDSAMPLEGFF}" \
                  > "${CHROMRENAME_HAPLIFTEDSAMPLEGFF}"
            else
              echo "Skipping creation of ${CHROMRENAME_HAPLIFTEDSAMPLEGFF}..."
            fi
            if [ ! -s "${CHROMRENAME_HAPLIFTEDSAMPLEGFF}" ]; then
              echo "File ${CHROMRENAME_HAPLIFTEDSAMPLEGFF} missing or empty for ${SAMPLE} h${H}."
              exit 1
            fi


            # also rename from the sample patho FA file...
            CHROMRENAME_SAMPLEPATHFA="${SAMPLEPATHFA}.chromrename.fa"
            if [ ! -s "${CHROMRENAME_SAMPLEPATHFA}" ]; then
               echo "Creating renamed chrom sample patho FA file ${CHROMRENAME_SAMPLEPATHFA}..."
               sed "1s/>.*/>${CHROM}/" "${SAMPLEPATHFA}" > "${CHROMRENAME_SAMPLEPATHFA}"
            else
              echo "Skipping creation of ${CHROMRENAME_SAMPLEPATHFA}..."
            fi
            if [ ! -s "${CHROMRENAME_SAMPLEPATHFA}" ]; then
              echo "File ${CHROMRENAME_SAMPLEPATHFA} missing or empty for ${SAMPLE} h${H}."
              exit 1
            fi

            # then perform liftover to pathogenic variant coordinates
            PATHOLIFTEDSAMPLEGFF="${PREFIX_PTH}.patholifted.local.MANE.GRCh38.v1.3.ensembl_genomic.gff"
            #echo "Running CrossMap: CrossMap gff ${SAMPLEPATHCHAIN} ${CHROMRENAME_HAPLIFTEDSAMPLEGFF} ${PATHOLIFTEDSAMPLEGFF}"
            CrossMap gff \
              "${SAMPLEPATHCHAIN}" \
              "${CHROMRENAME_HAPLIFTEDSAMPLEGFF}" \
              "${PATHOLIFTEDSAMPLEGFF}"

            # translate to CDS and protein:
            echo "Now translating to protein"
            SAMPLEPATH_CDS="${PREFIX_PTH}.${CHR_POS_ALT}.cds.fa"
            SAMPLEPATH_PROT="${PREFIX_PTH}.${CHR_POS_ALT}.protein.fa"
            #echo "gffread ${PATHOLIFTEDSAMPLEGFF} -g ${PREFIX_PTH}.${CHR_POS_ALT}.renamed.fa -x ${PREFIX_PTH}.${CHR_POS_ALT}.cds.fa -y ${PREFIX_PTH}.${CHR_POS_ALT}.protein.fa -S"
            gffread "${PATHOLIFTEDSAMPLEGFF}" -g "${CHROMRENAME_SAMPLEPATHFA}" -x "${SAMPLEPATH_CDS}" -y "${SAMPLEPATH_PROT}" -S


        done
    done < "${DATADIR}/${SAMPLENAMES}"


done

# liftover MANE coordinates:



# % diff HG00620.h1.chr16_14348.protein.fa tmp/MEFV/samples/HG00620.h1.protein.fa
# 11c11
# < DGPQRFDSCIIVLGSPSFLSGRRYWEVEVGDKTAWILGACKTSISRKGNMTLSPENGYWVVIMVKENEYQ
# ---
# > DGPQRFDSCIIVLGSPSFLSGRRYWEVEVGDKTAWILGACKTSISRKGNMTLSPENGYWVVIMMKENEYQ