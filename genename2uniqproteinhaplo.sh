#!/bin/bash
set -euo pipefail

# tools required:
#% samtools --version
#samtools 1.23.1
#% bcftools --version
#bcftools 1.23.1
#% gffread --version
#0.12.9
#% wget --version
#GNU Wget 1.25.0 built on darwin24.1.0.
#% zgrep --version
#zgrep (BSD grep, GNU compatible) 2.6.0-FreeBSD

if [ -z "$1" ] || [ -z "$2" ]; then
    echo "Error: Missing required arguments."
    echo "Usage: $0 <GENE_NAME> <REFSEQ_PROTEIN_ID> [SIGNAL_CUTOFF]"
    exit 1
fi

GENE=$1
PROTEIN_ID=$2

# Parse optional second argument for removing signal peptides from the front of
# reconstructed amino acid sequences before aggregation (SIGNAL_CUTOFF, default to 0)
SIGNAL_CUTOFF="${3:-0}"

DATADIR=data
TMPDIR=tmp
mkdir -p "${DATADIR}"
mkdir -p "${TMPDIR}"

GRCh38REFFA="GRCh38_full_analysis_set_plus_decoy_hla.fa"
GRCh38REFFAINDEX="GRCh38_full_analysis_set_plus_decoy_hla.fa.fai"
GRCh38MANEGFF="MANE.GRCh38.v1.5.refseq_genomic.gff"
GRCh38MANEGFFGZ="${GRCh38MANEGFF}.gz"

echo "Starting script with gene name = ${GENE}"
echo "RefSeq protein ID = ${PROTEIN_ID}"
echo "Signal peptide cutoff = ${SIGNAL_CUTOFF}"

if [ ! -s "${DATADIR}/${GRCh38REFFA}" ]; then
    echo "File ${DATADIR}/${GRCh38REFFA} not found, downloading..."
    wget -P "${DATADIR}" "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/${GRCh38REFFA}"
else
    echo "File ${DATADIR}/${GRCh38REFFA} exists, skipping download..."
fi
if [[ ! -s "${DATADIR}/${GRCh38REFFA}" ]]; then
    echo "Error: File ${DATADIR}/${GRCh38REFFA} not found or empty."
    exit 1
fi

if [ ! -s "${DATADIR}/${GRCh38REFFAINDEX}" ]; then
    echo "File ${DATADIR}/${GRCh38REFFAINDEX} not found, downloading..."
    wget -P "${DATADIR}" "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/${GRCh38REFFAINDEX}"
else
    echo "File ${DATADIR}/${GRCh38REFFAINDEX} exists, skipping download..."
fi
if [[ ! -s "${DATADIR}/${GRCh38REFFAINDEX}" ]]; then
    echo "Error: File ${DATADIR}/${GRCh38REFFAINDEX} not found or empty."
    exit 1
fi

if [ ! -s "${DATADIR}/${GRCh38MANEGFFGZ}" ]; then
    echo "File ${DATADIR}/${GRCh38MANEGFFGZ} not found, downloading..."
    wget -P "${DATADIR}" "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5/${GRCh38MANEGFFGZ}"
else
    echo "File ${DATADIR}/${GRCh38MANEGFFGZ} exists, skipping download..."
fi
if [[ ! -s "${DATADIR}/${GRCh38MANEGFFGZ}" ]]; then
    echo "Error: File ${DATADIR}/${GRCh38MANEGFFGZ} not found or empty."
    exit 1
fi

filter_sv_for_gene() {
    set -o pipefail
    local chrom="$1"
    local start="$2"
    local end="$3"

    bcftools view \
        -r "${chrom}:${start}-${end}" \
        "${DATADIR}/${THGENSAMPLES}" |
    awk -v start="$start" \
        -v end="$end" '
    BEGIN { FS=OFS="\t" }

    /^#/ {
        print
        next
    }

    {
        alt = $5

        # Ordinary variants
        if (alt !~ /^</) {
            print
            next
        }

        # Symbolic alleles
        svtype = ""
        svend = ""

        n = split($8, info, ";")

        for (i = 1; i <= n; i++) {
            if (info[i] ~ /^SVTYPE=/)
                svtype = substr(info[i], 8)

            if (info[i] ~ /^END=/)
                svend = substr(info[i], 5) + 0
        }

        # a DEL, supported by bcftools consensus
        if (svtype == "DEL") {
            print
            next
        }

        # a DUP encompassing the entire gene
        # ignore because the localized FASTA is not intended to represent the extra copy
        if (svtype == "DUP" &&  $2 <= start && svend >= end) {
            next
        }

        # a DUP completely outside region
        if (svtype == "DUP" && (svend < start || $2 > end)) {
            next
        }

        # insertions, throw error
        if (svtype == "INS" || svtype ~ /^INS:/) {
            print "ERROR: unsupported SV: " svtype " " $3 " " $1 ":" $2 "-" svend > "/dev/stderr"
            print "OFFENDING VCF RECORD:\n" $0 > "/dev/stderr"
            exit 2
        }

        # Everything else: also fail
        print "ERROR: unsupported SV: " svtype " " $3 " " $1 ":" $2 "-" svend > "/dev/stderr"
        print "OFFENDING VCF RECORD:\n" $0 > "/dev/stderr"
        exit 2
    }'
}

# create gene-specific MANE GFF
GENEDIR="${TMPDIR}/${GENE}"
mkdir -p "${GENEDIR}"
if [ ! -s "${GENEDIR}/${GRCh38MANEGFF}" ]; then
    echo "File ${GENEDIR}/${GRCh38MANEGFF} not found, extracting for Gene: ${GENE}, Protein ID: ${PROTEIN_ID}..."
    # Pass 1: find the parent ID
    # Pass 2: filter & print
    # 'exon' technically not needed but good to keep compliance
    zgrep -E "gene=${GENE}(;|$)" "${DATADIR}/${GRCh38MANEGFFGZ}" | awk -v prot="${PROTEIN_ID}" '
    BEGIN { FS = OFS = "\t" }
    {
        lines[NR] = $0
        if ($9 ~ ("protein_id=" prot "(;|$)")) {
            split($9, attrs, ";")
            for (i in attrs) {
                if (attrs[i] ~ /^Parent=/) {
                    target_tx[substr(attrs[i], 8)] = 1
                }
            }
        }
    }
    END {
        for (i = 1; i <= NR; i++) {
            $0 = lines[i]
            if ($3 == "gene") {
                print $0
                continue
            }
            split($9, attrs, ";")
            id = ""; parent = ""
            for (j in attrs) {
                if (attrs[j] ~ /^ID=/) id = substr(attrs[j], 4)
                if (attrs[j] ~ /^Parent=/) parent = substr(attrs[j], 8)
            }
            if ((id in target_tx) || (parent in target_tx)) {
                print $0
            }
        }
    }' > "${GENEDIR}/${GRCh38MANEGFF}"
else
    echo "File ${GENEDIR}/${GRCh38MANEGFF} exists, skipping extraction..."
fi
if [ ! -s "${GENEDIR}/${GRCh38MANEGFF}" ]; then
    echo "Error: File '${GENEDIR}/${GRCh38MANEGFF}' is empty or does not exist. Bad gene name?"
    exit 1
fi


# validate it: check for required GFF feature types in column 3
MISSING_TYPES=$(awk -F'\t' '
BEGIN { gene=mRNA=CDS=exon=1 }
$3 == "gene" { gene=0 }
$3 == "mRNA" { mRNA=0 }
$3 == "CDS"  { CDS=0 }
$3 == "exon" { exon=0 }
END {
    if (gene) print "gene"
    if (mRNA) print "mRNA"
    if (CDS)  print "CDS"
    if (exon) print "exon"
}' "${GENEDIR}/${GRCh38MANEGFF}")

if [ -n "${MISSING_TYPES}" ]; then
    echo "Error: Validation failed for '${GENEDIR}/${GRCh38MANEGFF}'." >&2
    echo "Missing required feature type(s): $(echo ${MISSING_TYPES} | tr '\n' ' ')" >&2
    # Optional: remove invalid file so it gets retried next time
    rm -f "${GENEDIR}/${GRCh38MANEGFF}"
    exit 1
fi
echo "Validation successful for ${GENE} (${PROTEIN_ID})."

# get region
DATA=$(awk -F'\t' '$3=="mRNA" {print $1, $4, $5; exit}' "${GENEDIR}/${GRCh38MANEGFF}")
if [ -z "$DATA" ]; then
    echo "Error: Could not find a 'mRNA' feature in ${GENEDIR}/${GRCh38MANEGFF}."
    exit 1
fi
read CHR START END <<< "$DATA"
if [[ -z "$CHR" || -z "$START" || -z "$END" ]]; then
    echo "Error: Coordinate extraction failed (one or more fields are empty)."
    exit 1
fi
REGION="${CHR}:${START}-${END}"
echo "Extracted region for this gene is $REGION"

# now we know which 1000G VCF file we need, download if not present
THGENSAMPLES="1kGP_high_coverage_Illumina.${CHR}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
THGENSAMPLESINDEX="${THGENSAMPLES}.tbi"
BASEURL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV"
if [ ! -s "${DATADIR}/${THGENSAMPLES}" ]; then
    echo "File ${DATADIR}/${THGENSAMPLES} not found, downloading..."
    wget -P "${DATADIR}" "${BASEURL}/${THGENSAMPLES}"
else
    echo "File ${DATADIR}/${THGENSAMPLES} exists, skipping download..."
fi
if [[ ! -s "${DATADIR}/${THGENSAMPLES}" ]]; then
    echo "Error: File ${DATADIR}/${THGENSAMPLES} not found or empty."
    exit 1
fi
if [ ! -s "${DATADIR}/${THGENSAMPLESINDEX}" ]; then
    echo "File ${DATADIR}/${THGENSAMPLESINDEX} not found, downloading..."
    wget -P "${DATADIR}" "${BASEURL}/${THGENSAMPLESINDEX}"
else
    echo "File ${DATADIR}/${THGENSAMPLESINDEX} exists, skipping download..."
fi
if [[ ! -s "${DATADIR}/${THGENSAMPLESINDEX}" ]]; then
    echo "Error: File ${DATADIR}/${THGENSAMPLESINDEX} not found or empty."
    exit 1
fi

# get sample names
SAMPLENAMES="1kGP_sample_names.txt"
if [ ! -s "${DATADIR}/${SAMPLENAMES}" ]; then
    echo "File ${DATADIR}/${SAMPLENAMES} not found, extracting from 1KG sample VCF..."
    bcftools query -l "${DATADIR}/$THGENSAMPLES" > "${DATADIR}/${SAMPLENAMES}"
else
    echo "File ${DATADIR}/${SAMPLENAMES} exists, skipping extraction..."
fi
if [[ ! -s "${DATADIR}/${SAMPLENAMES}" ]]; then
    echo "Error: File ${DATADIR}/${SAMPLENAMES} not found or empty."
    exit 1
fi

# regional genomic reference file for this gene
GENEREFFA="${GENEDIR}/${GRCh38REFFA}"
if [ ! -s "${GENEREFFA}" ]; then
    echo "File ${GENEREFFA} not found, extracting from reference genome sequence..."
    samtools faidx "${DATADIR}/${GRCh38REFFA}" "$REGION" > "${GENEREFFA}"
    samtools faidx "${GENEREFFA}" # index it
else
    echo "File ${GENEREFFA} exists, skipping extraction..."
fi
if [[ ! -s "${GENEREFFA}" ]]; then
    echo "Error: File ${GENEREFFA} not found or empty."
    exit 1
fi

GENESAMPLESDIR="${GENEDIR}/samples"
mkdir -p "$GENESAMPLESDIR"

# loop over samples and haplotype 1/2
# perhaps change    -e 'ALT~"<INS>"'   to    -e 'ALT~"<.*>"' 
while read -r SAMPLE; do
    [ -z "$SAMPLE" ] && continue # skip empty lines
    for H in 1 2; do

        PREFIX="${GENESAMPLESDIR}/${SAMPLE}.h${H}"

        if [[ -s "${PREFIX}.fa" && -s "${PREFIX}.chain" ]]; then
            echo "File ${PREFIX}.fa and ${PREFIX}.chain exists and not empty, skipping ${SAMPLE} h${H}..."
            continue
        fi

        echo "Filtering VCF for sample=${SAMPLE} haplotype=${H} in region $CHR:$START-$END"

        FILTERED_VCF="${GENESAMPLESDIR}/${SAMPLE}.h${H}.filtered.vcf.gz"

        filter_sv_for_gene \
            "$CHR" \
            "$START" \
            "$END" |
            bgzip -c > "${FILTERED_VCF}"

        bcftools index -f "${FILTERED_VCF}"

        echo "Creating consensus and chain file for sample=${SAMPLE}, haplotype=${H}..."

        bcftools consensus \
            -f "${GENEREFFA}" \
            -s "$SAMPLE" \
            -H "$H" \
            --regions-overlap 1 \
            --chain "${PREFIX}.chain" \
            "${FILTERED_VCF}" \
            > "${PREFIX}.fa"

      if [ ! -s "${PREFIX}.fa" ]; then
          echo "File ${PREFIX}.fa empty for ${SAMPLE} h${H}."
          exit 1
      fi
      if [ ! -s "${PREFIX}.chain" ]; then
          echo "File ${PREFIX}.chain empty for ${SAMPLE} h${H}."
          exit 1
      fi
    done
done < "${DATADIR}/${SAMPLENAMES}"

# create indices for sample FAs
while read -r SAMPLE; do
    [ -z "$SAMPLE" ] && continue # skip empty lines
    for H in 1 2; do
        PREFIX="${GENESAMPLESDIR}/${SAMPLE}.h${H}"
        if [ -s "${PREFIX}.fa.fai" ]; then
            echo "File ${PREFIX}.fa.fai exists and not empty, skipping indexing of ${SAMPLE} h${H}..."
            continue
        fi
        echo "Creating index for sample=${SAMPLE}, haplotype=${H}..."
        samtools faidx "${PREFIX}.fa"
    done
done < "${DATADIR}/${SAMPLENAMES}"

# create lifted GFF per sample
while read -r SAMPLE; do
    [ -z "$SAMPLE" ] && continue # skip empty lines
    for H in 1 2; do
        PREFIX="${GENESAMPLESDIR}/${SAMPLE}.h${H}"
        LIFTEDMANE="${GENESAMPLESDIR}/${SAMPLE}.h${H}.lifted.${GRCh38MANEGFF}"
        if [ -s "${LIFTEDMANE}" ]; then
            echo "File ${LIFTEDMANE} exists and not empty, skipping ${SAMPLE} h${H}..."
            continue
        fi
        echo "Creating lifted GFF for sample=${SAMPLE}, haplotype=${H}..."
        ./project_gff_by_chain_preserve_features.py \
            "${PREFIX}.chain" \
            "${GENEDIR}/${GRCh38MANEGFF}" \
            "${LIFTEDMANE}"

        if [[ ! -s "${LIFTEDMANE}" ]]; then
            echo "Error: File ${LIFTEDMANE} not found or empty."
            exit 1
        fi
    done
done < "${DATADIR}/${SAMPLENAMES}"


# fix coordinates
while read -r SAMPLE; do
    [ -z "$SAMPLE" ] && continue # skip empty lines
    for H in 1 2; do
        PREFIX="${GENESAMPLESDIR}/${SAMPLE}.h${H}"
        LIFTEDMANE="${GENESAMPLESDIR}/${SAMPLE}.h${H}.lifted.${GRCh38MANEGFF}"
        LIFTEDLOCALMANE="${GENESAMPLESDIR}/${SAMPLE}.h${H}.lifted.local.${GRCh38MANEGFF}"
        if [ -s "${LIFTEDLOCALMANE}" ]; then
            echo "File ${LIFTEDLOCALMANE} exists and not empty, skipping ${SAMPLE} h${H}..."
            continue
        fi
        NSEQ=$(wc -l < "${PREFIX}.fa.fai")
        if [ "$NSEQ" -ne 1 ]; then
            echo "Error: expected exactly 1 FASTA record in ${PREFIX}.fa.fai, found ${NSEQ}." >&2
            exit 1
        fi
        echo "Creating local lifted GFF with coordinate fix for sample=${SAMPLE}, haplotype=${H}"
        SEQNAME=$(cut -f1 "${PREFIX}.fa.fai")
        awk -v seq="$SEQNAME" -v start="$START" '
        BEGIN {
            FS = OFS = "\t"
        }

        /^#/ {
            print
            next
        }

        {
            $1 = seq
            $4 = $4 - start + 1
            $5 = $5 - start + 1
            print
        }
        ' "${LIFTEDMANE}" > "${LIFTEDLOCALMANE}"
        if [[ ! -s "${LIFTEDLOCALMANE}" ]]; then
            echo "Error: File ${LIFTEDLOCALMANE} not found or empty."
            exit 1
        fi
    done
done < "${DATADIR}/${SAMPLENAMES}"



# create CDS and translate to protein using lifted local GFF and sample FA
while read -r SAMPLE; do
    [ -z "$SAMPLE" ] && continue # skip empty lines
    for H in 1 2; do
        PREFIX="${GENESAMPLESDIR}/${SAMPLE}.h${H}"
        LIFTEDLOCALMANE="${GENESAMPLESDIR}/${SAMPLE}.h${H}.lifted.local.${GRCh38MANEGFF}"

        # Check if output files already exist
        if [ -s "${PREFIX}.cds.fa" ] && [ -s "${PREFIX}.protein.fa" ]; then
            echo "CDS and Protein files for ${SAMPLE} h${H} already exist and not empty. Skipping..."
            continue
        fi

        echo "Creating CDS and protein translation for sample=${SAMPLE}, haplotype=${H}..."
        gffread "${LIFTEDLOCALMANE}" -g "${PREFIX}.fa" -x "${PREFIX}.cds.fa" -y "${PREFIX}.protein.fa" -S

        if [ ! -s "${PREFIX}.protein.fa" ]; then
            echo "Error: Translation failed for ${SAMPLE} h${H}"
            exit 1
        fi

    done
done < "${DATADIR}/${SAMPLENAMES}"



# collect all protein sequences in one TSV - truncate sequences with stop gain
PROTPERSAMPLE="${GENEDIR}/${GENE}_protein_sequences_per_sample.tsv"
if [ -s "${PROTPERSAMPLE}" ]; then
    echo "File ${PROTPERSAMPLE} exists and not empty, skipping collection of sequences..."
else
  : > "${PROTPERSAMPLE}"
  for f in "${GENESAMPLESDIR}/"*.protein.fa; do
      id=$(basename "$f" .protein.fa)

      # Use awk to build the full sequence, truncate at stop, and cut signal peptide
      seq=$(awk -v cutoff="$SIGNAL_CUTOFF" '
          /^>/ {next}
          {
              gsub(/[[:space:]]/, "", $0)
              full_seq = full_seq toupper($0)
          }
          END {
              stop_pos = index(full_seq, "*")
              if (stop_pos > 0) {
                  # Cut off everything from the * onwards
                  full_seq = substr(full_seq, 1, stop_pos - 1)
              }

              # Trim signal peptide from the front (start at position cutoff + 1)
              if (cutoff > 0 && length(full_seq) > cutoff) {
                  full_seq = substr(full_seq, cutoff + 1)
              } else if (cutoff > 0 && length(full_seq) <= cutoff) {
                  full_seq = ""
              }

              print full_seq
          }
      ' "$f")

      if [ -z "$seq" ]; then
          # If stop codon or signal cutoff renders sequence empty
          echo "ERROR: Protein sequence for $id is empty (possible immediate stop codon or cutoff exceeds length)" >&2
          exit 1
      fi

      printf "%s\t%s\n" "$id" "$seq" >> "${PROTPERSAMPLE}"
  done
  echo "Collection and truncation complete."
fi

# group identical protein sequences together and store in results
RESULTDIR="results/${GENE}"
mkdir -p "${RESULTDIR}"
UNIQUEPROTGROUPED="${RESULTDIR}/${GENE}_protein_haplotypes.fa"

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

          printf ">%s_protein_variant_%d|n=%d|samples=%s\n", gene, i, count[seq], samples[seq]

          for (j = 1; j <= length(seq); j += 60) {
              print substr(seq, j, 60)
          }
      }
  }
  ' "${PROTPERSAMPLE}" > "${UNIQUEPROTGROUPED}"
fi

echo "...all done"
