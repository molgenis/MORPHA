#!/bin/bash

# Check for required arguments
if [ -z "$1" ] || [ -z "$2" ]; then
    echo "Usage: $0 <gene name> <pathovariants local vcf.gz>"
    exit 1
fi

GENE=$1 #e.g. "MEFV"
PATHOVARIANTSLOCAL=$2 #e.g. "results/MEFV/infevers/MEFV-infevers-one-inframedel.vcf.gz"

echo "Starting script with gene name = ${GENE}"
echo "Using pathovariants local file = ${PATHOVARIANTSLOCAL}"

DATADIR=data
TMPDIR=tmp
GENEDIR="${TMPDIR}/${GENE}"
SAMPLENAMES="1kGP_sample_names.txt"
GENESAMPLESHAPLODIR="${GENEDIR}/samples"
GENESAMPLESPATHODIR="${GENEDIR}/samples_patho"
mkdir -p "${GENESAMPLESPATHODIR}"

localize_chain() {
    local GLOBAL_CHAIN=$1
    local LOCAL_CHAIN=$2

    awk '
    BEGIN {FS=OFS=" "}

    /^chain/ {
        # Reference (Target) coordinates
        r_start = $6
        $6 = 0
        $7 = $7 - r_start    # Shift target end to local 0-based
        $4 = $7              # Update target total size to local length

        # Sample (Query) coordinates
        s_start = $11
        $11 = 0
        $12 = $12 - s_start  # Shift query end to local 0-based
        $9 = $12             # Update query total size to local length

        print $0
        next
    }

    # Keep the alignment data blocks exactly as they are
    { print $0 }
    ' "$GLOBAL_CHAIN" > "$LOCAL_CHAIN"
}



bcftools query -f '%CHROM:%POS:%ID:%ALT\n' "${PATHOVARIANTSLOCAL}" | while read -r variant; do

    # DEBUG
    #variant="chr16:1377:c.2081_2083del:T" # --> VCF variant "chr16	1377	c.2081_2083del	TTCA	T" causes exon 10 drop in CrossMap GFF liftover, so good test

    # Isolate the position, chromosome and alternative allele
    # NOTE: we assume 1 alt per line, check if this is really true!!
    CHROM=$(echo "$variant" | cut -d':' -f1)
    POS=$(echo "$variant" | cut -d':' -f2)
    ID=$(echo "$variant" | cut -d':' -f3)
    ALT=$(echo "$variant" | cut -d':' -f4)
    CHR_POS_ALT="${CHROM}_${POS}_${ALT}"

    echo "Processing variant at $CHROM:$POS with ALT allele $ALT..."

    # Make a temporary single-variant VCF on the fly plus index
    SINGLE_VAR_TMP_VCF="${GENEDIR}/tmp_singlevar_${CHR_POS_ALT}.vcf.gz"
    if [ ! -s "${SINGLE_VAR_TMP_VCF}" ]; then
      echo "Creating ${SINGLE_VAR_TMP_VCF}..."
      bcftools view -O z -i "CHROM=='${CHROM}' && POS==${POS} && ID=='${ID}' && ALT=='${ALT}'" "${PATHOVARIANTSLOCAL}" > "${SINGLE_VAR_TMP_VCF}"
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

          # debug
          # SAMPLE="HG00663"
          # H="2"

            echo "Going to work on sample: ${SAMPLE}.h${H}"
            PREFIX_HAP="${GENESAMPLESHAPLODIR}/${SAMPLE}.h${H}"
            PREFIX_PTH="${GENESAMPLESPATHODIR}/${SAMPLE}.h${H}"

            # protein file output location, can already check and skip if present!
            SAMPLEPATH_PROT="${PREFIX_PTH}.${CHR_POS_ALT}.${ID}.protein.fa"
            if [ -s "${SAMPLEPATH_PROT}" ]; then
              echo "File ${SAMPLEPATH_PROT} already present for ${SAMPLE} h${H}, skipping to next..."
              continue
            fi

            ORIGINALSAMPLEFA="${PREFIX_HAP}.fa"

            if [ ! -s "${ORIGINALSAMPLEFA}" ]; then
              echo "File ${ORIGINALSAMPLEFA} missing or empty for ${SAMPLE} h${H}."
              exit 1
            fi

            # Temporarily simplify FASTA header so it says ">chr16" instead of ">chr16:3242..."
            # to enable bcftools consensus matching of VCF and FASTA
            ORIGINALSAMPLEFAHEADER=$(grep ">" "${ORIGINALSAMPLEFA}" | head -n 1)
            echo "Extracted original FASTA header: ${ORIGINALSAMPLEFAHEADER}"
            RENAMEDHEADERSAMPLEFA="${PREFIX_PTH}.hdr.fa"
            if [ ! -s "${RENAMEDHEADERSAMPLEFA}" ]; then
              echo "Creating ${RENAMEDHEADERSAMPLEFA}..."
              sed "1s/>.*/>${CHROM}/" "${ORIGINALSAMPLEFA}" > "${RENAMEDHEADERSAMPLEFA}"
            else
              echo "Skipping creation of ${RENAMEDHEADERSAMPLEFA}..."
            fi
            if [ ! -s "${RENAMEDHEADERSAMPLEFA}" ]; then
              echo "File ${RENAMEDHEADERSAMPLEFA} missing or empty for ${SAMPLE} h${H}."
              exit 1
            fi

            # coordinates from FASTA header
            COORDS="${ORIGINALSAMPLEFAHEADER#*:}"
            START="${COORDS%-*}"
            END="${COORDS#*-}"
            REGION_LEN=$(( END - START + 1 ))
            SAMPLE_START=0

            # create  local chain file
            LOCAL_CHAIN="${PREFIX_PTH}.local.chain"
            localize_chain "${PREFIX_HAP}.chain" "${LOCAL_CHAIN}"
            # "${START}" "${SAMPLE_START}" "${REGION_LEN}"

            # liftover the tmp pathogenic variant VCF to match sample coordinates
            LIFTOVER_SINGLE_VAR_TMP_VCF="${PREFIX_PTH}.${CHR_POS_ALT}.liftover.vcf"
            if [ ! -s "${LIFTOVER_SINGLE_VAR_TMP_VCF}" ]; then
              echo "Creating ${LIFTOVER_SINGLE_VAR_TMP_VCF}..."
              CrossMap vcf "${LOCAL_CHAIN}" "${SINGLE_VAR_TMP_VCF}" "${RENAMEDHEADERSAMPLEFA}" "${LIFTOVER_SINGLE_VAR_TMP_VCF}"
            else
              echo "Skipping creation of ${LIFTOVER_SINGLE_VAR_TMP_VCF}..."
            fi
            if [ ! -s "${LIFTOVER_SINGLE_VAR_TMP_VCF}" ]; then
              echo "File ${LIFTOVER_SINGLE_VAR_TMP_VCF} missing or empty for ${SAMPLE} h${H}."
              exit 1
            fi

            # if ALT is equal to haplotype base (e.g. non-reference, variant already present in sample)
            # we want to flag as no-op
            NOOP=1
            if awk -F'\t' '
                BEGIN {
                    found = 0
                    valid = 0
                }

                !/^#/ {
                    found = 1
                    ref = $4
                    alt = $5

                    if (ref != "" && ref != "." &&
                        alt != "" && alt != "." &&
                        ref != alt) {
                        valid = 1
                    }

                    exit
                }

                END {
                    exit !(found && valid)
                }
            ' "${LIFTOVER_SINGLE_VAR_TMP_VCF}"; then
                NOOP=0
            fi
            if [ "${NOOP}" -eq 1 ]; then
                echo "REF/ALT invalid, missing, or identical; treating as no-op."
            else
                echo "REF and ALT are present and different; applying variant."
            fi

            # create GZ and index regardless for consistency
            if [ ! -s "${LIFTOVER_SINGLE_VAR_TMP_VCF}.gz" ]; then
              echo "Creating ${LIFTOVER_SINGLE_VAR_TMP_VCF}.gz..."
              bgzip -c "${LIFTOVER_SINGLE_VAR_TMP_VCF}" > "${LIFTOVER_SINGLE_VAR_TMP_VCF}.gz"
              bcftools index "${LIFTOVER_SINGLE_VAR_TMP_VCF}.gz"
            else
              echo "Skipping creation of ${LIFTOVER_SINGLE_VAR_TMP_VCF}.gz..."
            fi
            if [ ! -s "${LIFTOVER_SINGLE_VAR_TMP_VCF}.gz" ]; then
              echo "File ${LIFTOVER_SINGLE_VAR_TMP_VCF}.gz missing or empty for ${SAMPLE} h${H}."
              exit 1
            fi


            # run consensus to create FA file with pathogenic variant
            SAMPLEPATHCHAIN="${PREFIX_PTH}.${CHR_POS_ALT}.chain"
            SAMPLEPATHFA="${PREFIX_PTH}.${CHR_POS_ALT}.fa"

            # if noop, allele is already present, so sequence is unchanged
            if [ "${NOOP}" -eq 1 ]; then
                echo "No allele change - copying existing haplotype FASTA"
                cp "${RENAMEDHEADERSAMPLEFA}" "${SAMPLEPATHFA}"
            else
              if [[ ! -s "${SAMPLEPATHCHAIN}" || ! -s "${SAMPLEPATHFA}" ]]; then
                  echo "Creating ${SAMPLEPATHCHAIN} and ${SAMPLEPATHFA} sequence consensus..."
                  bcftools consensus \
                      -f "${RENAMEDHEADERSAMPLEFA}" \
                      -c "${SAMPLEPATHCHAIN}" \
                      "${LIFTOVER_SINGLE_VAR_TMP_VCF}.gz" \
                      > "${SAMPLEPATHFA}"
              else
                echo "${SAMPLEPATHCHAIN} and ${SAMPLEPATHFA} present, skipping consensus..."
              fi
              if [[ ! -s "${SAMPLEPATHCHAIN}" || ! -s "${SAMPLEPATHFA}" ]]; then
                 echo "${SAMPLEPATHCHAIN} and/or ${SAMPLEPATHFA} are missing."
                 exit 1
              fi
            fi

            # now also need to lift GFF files, first create tmp GFF and then do liftover
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

            # now perform liftover to pathogenic variant coordinates - cannot use CrossMap here, it does not handle indels correctly!
            PATHOLIFTEDSAMPLEGFF="${PREFIX_PTH}.patholifted.local.MANE.GRCh38.v1.3.ensembl_genomic.gff"

            # if noop, no liftover magic is needed, just reuse existing one
            if [ "${NOOP}" -eq 1 ]; then
                echo "No allele change - copying existing GFF"
                cp "${CHROMRENAME_HAPLIFTEDSAMPLEGFF}" "${PATHOLIFTEDSAMPLEGFF}"
            else
                ./project_gff_by_chain_preserve_features.py \
                  "${SAMPLEPATHCHAIN}" \
                  "${CHROMRENAME_HAPLIFTEDSAMPLEGFF}" \
                  "${PATHOLIFTEDSAMPLEGFF}"
            fi


            # translate to CDS and protein:
            echo "Now translating to protein"
            SAMPLEPATH_CDS="${PREFIX_PTH}.${CHR_POS_ALT}.cds.fa"
            #echo "gffread ${PATHOLIFTEDSAMPLEGFF} -g ${PREFIX_PTH}.${CHR_POS_ALT}.renamed.fa -x ${PREFIX_PTH}.${CHR_POS_ALT}.cds.fa -y ${PREFIX_PTH}.${CHR_POS_ALT}.protein.fa -S"
            gffread "${PATHOLIFTEDSAMPLEGFF}" -g "${SAMPLEPATHFA}" -x "${SAMPLEPATH_CDS}" -y "${SAMPLEPATH_PROT}" -S


            # validations
            if [ ! -s "${SAMPLEPATH_PROT}" ]; then
                echo "ERROR: gffread did not produce a protein for ${SAMPLE}.h${H}." >&2
                exit 1
            fi
            PROTEIN_SEQUENCE=$(
                awk '!/^>/ {
                    gsub(/[[:space:]]/, "")
                    printf "%s", toupper($0)
                }' "${SAMPLEPATH_PROT}"
            )
            if [[ "${PROTEIN_SEQUENCE}" == *X* ]]; then
                echo "ERROR: Protein contains X for ${SAMPLE}.h${H}, variant ${ID}." >&2
                echo "Keeping intermediate files for debugging." >&2
                exit 1
            fi

            # need to clean up since we process so many samples x haplotypes x mutations, only keep protein file
            if [ -n "${SAMPLEPATH_CDS}" ] && [ -f "${SAMPLEPATH_CDS}" ]; then
                rm "${SAMPLEPATH_CDS}"
            fi
            if [ -n "${SAMPLEPATHCHAIN}" ] && [ -f "${SAMPLEPATHCHAIN}" ]; then
                rm "${SAMPLEPATHCHAIN}"
            fi
            if [ -n "${SAMPLEPATHFA}" ] && [ -f "${SAMPLEPATHFA}" ]; then
                rm "${SAMPLEPATHFA}"
                rm "${SAMPLEPATHFA}.fai"
            fi
            if [ -n "${RENAMEDHEADERSAMPLEFA}" ] && [ -f "${RENAMEDHEADERSAMPLEFA}" ]; then
                rm "${RENAMEDHEADERSAMPLEFA}"
                rm "${RENAMEDHEADERSAMPLEFA}.fai"
            fi
            if [ -n "${CHROMRENAME_HAPLIFTEDSAMPLEGFF}" ] && [ -f "${CHROMRENAME_HAPLIFTEDSAMPLEGFF}" ]; then
                rm "${CHROMRENAME_HAPLIFTEDSAMPLEGFF}"
            fi
            if [ -n "${PATHOLIFTEDSAMPLEGFF}" ] && [ -f "${PATHOLIFTEDSAMPLEGFF}" ]; then
                rm "${PATHOLIFTEDSAMPLEGFF}"
                rm "${PATHOLIFTEDSAMPLEGFF}".*.unmap
            fi
            if [ -n "${LIFTOVER_SINGLE_VAR_TMP_VCF}" ] && [ -f "${LIFTOVER_SINGLE_VAR_TMP_VCF}" ]; then
                rm "${LIFTOVER_SINGLE_VAR_TMP_VCF}"
            fi
            if [ -n "${LIFTOVER_SINGLE_VAR_TMP_VCF}.gz" ] && [ -f "${LIFTOVER_SINGLE_VAR_TMP_VCF}.gz" ]; then
                rm "${LIFTOVER_SINGLE_VAR_TMP_VCF}.gz"
            fi
            if [ -n "${LIFTOVER_SINGLE_VAR_TMP_VCF}.gz.csi" ] && [ -f "${LIFTOVER_SINGLE_VAR_TMP_VCF}.gz.csi" ]; then
                rm "${LIFTOVER_SINGLE_VAR_TMP_VCF}.gz.csi"
            fi
            if [ -n "${LIFTOVER_SINGLE_VAR_TMP_VCF}.unmap" ] && [ -f "${LIFTOVER_SINGLE_VAR_TMP_VCF}.unmap" ]; then
                rm "${LIFTOVER_SINGLE_VAR_TMP_VCF}.unmap"
            fi
            if [ -n "${LOCAL_CHAIN}" ] && [ -f "${LOCAL_CHAIN}" ]; then
                rm "${LOCAL_CHAIN}"
            fi


        done
    done < "${DATADIR}/${SAMPLENAMES}"

    if [ -n "${SINGLE_VAR_TMP_VCF}" ] && [ -f "${SINGLE_VAR_TMP_VCF}" ]; then
        rm "${SINGLE_VAR_TMP_VCF}"
        rm "${SINGLE_VAR_TMP_VCF}.csi"
    fi

done

# collect all protein sequences in one TSV - truncate sequences with stop gain
PROTPERSAMPLEANDPATHO="${GENEDIR}/${GENE}_protein_sequences_per_sample_and_patho.tsv"
if [ -s "${PROTPERSAMPLEANDPATHO}" ]; then
    echo "File ${PROTPERSAMPLEANDPATHO} exists and not empty, skipping collection..."
else
  : > "${PROTPERSAMPLEANDPATHO}"
  for f in "${GENESAMPLESPATHODIR}/"*.protein.fa; do
      id=$(basename "$f" .protein.fa)

      # Use awk to build the full sequence, then truncate at the first asterisk
      seq=$(awk '
          /^>/ {next}
          {
              gsub(/[[:space:]]/, "", $0)
              full_seq = full_seq toupper($0)
          }
          END {
              stop_pos = index(full_seq, "*")
              if (stop_pos > 0) {
                  # Cut off everything from the * onwards
                  print substr(full_seq, 1, stop_pos - 1)
              } else {
                  print full_seq
              }
          }
      ' "$f")

      if [ -z "$seq" ]; then
          # Note: If a stop codon is at the very first position, the sequence becomes empty
          echo "Warning: Protein sequence for $id is empty (possible immediate stop codon)" >&2
          continue
      fi

      printf "%s\t%s\n" "$id" "$seq" >> "${PROTPERSAMPLEANDPATHO}"
  done
  echo "Collection and truncation complete."
fi



# group identical protein sequences together and store in results
RESULTDIR="results/${GENE}"
mkdir -p "${RESULTDIR}"
UNIQUEPROTGROUPED="${RESULTDIR}/${GENE}_protein_path_haplotypes.fa"

if [ -s "${UNIQUEPROTGROUPED}" ]; then
    echo "File ${UNIQUEPROTGROUPED} exists and not empty, skipping grouping of sequences..."
else
  echo "Grouping protein sequence into ${UNIQUEPROTGROUPED}"
  awk -F'\t' -v gene="$GENE" '
  {
      id = $1
      seq = $2

      if (!(seq in seen)) {
          seen[seq] = ++n
          order[n] = seq
      }

      if (samples[seq] == "") {
          samples[seq] = id
      } else {
          samples[seq] = samples[seq] "," id
      }

      count[seq]++
  }

  END {
      for (i = 1; i <= n; i++) {
          seq = order[i]

          printf ">%s_protein_path_variant_%d|n=%d|samples=%s\n", gene, i, count[seq], samples[seq]

          for (j = 1; j <= length(seq); j += 60) {
              print substr(seq, j, 60)
          }
      }
  }
  ' "${PROTPERSAMPLEANDPATHO}" > "${UNIQUEPROTGROUPED}"
fi

echo "...all done"

# --> inspect file and number of haplotypes:
# tail -50 MEFV_protein_path_haplotypes.fa
# --> check based on sequences:
# awk '/^>/ {printf("\n%s\n",$0);next} {printf("%s",$0)} END {printf("\n")}' MEFV_protein_path_haplotypes.fa | grep -v "^>" | grep -v "^$" | sort | uniq | wc -l
# --> check based on headers:
# grep -c "^>" MEFV_protein_path_haplotypes.fa
