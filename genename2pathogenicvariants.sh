#!/bin/bash

[ -z "$1" ] && { echo "Missing argument: gene name"; exit 1; }

GENE=$1

DATADIR=data
TMPDIR=tmp
GENEDIR="${TMPDIR}/${GENE}"

CLINVARVCFGZ="clinvar_20260503.vcf.gz"
GENEMANEGFF="MANE.GRCh38.v1.5.refseq_genomic.gff"

if [ ! -s "${DATADIR}/${CLINVARVCFGZ}" ]; then
    echo "File ${DATADIR}/${CLINVARVCFGZ} not found, downloading..."
    wget -O "${DATADIR}/${CLINVARVCFGZ}" "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2026/${CLINVARVCFGZ}"
else
    echo "File ${DATADIR}/${CLINVARVCFGZ} exists, skipping download..."
fi
if [[ ! -s "${DATADIR}/${CLINVARVCFGZ}" ]]; then
    echo "Error: File ${DATADIR}/${CLINVARVCFGZ} not found or empty."
    exit 1
fi

# get region
DATA=$(awk -F'\t' '$3=="mRNA" {print $1, $4, $5; exit}' "${GENEDIR}/${GENEMANEGFF}")
if [ -z "$DATA" ]; then
    echo "Error: Could not find a 'mRNA' feature in ${GENEDIR}/${GENEMANEGFF}."
    exit 1
fi
read CHR START END <<< "$DATA"
if [[ -z "$CHR" || -z "$START" || -z "$END" ]]; then
    echo "Error: Coordinate extraction failed (one or more fields are empty)."
    exit 1
fi
REGION="${CHR}:${START}-${END}"
echo "Extracted region for this gene is $REGION"

CLINVARGENE="${GENEDIR}/ClinVarPathogenic"
CLINVARGENEVCF="${CLINVARGENE}.vcf"
CLINVARGENEVCFGZ="${CLINVARGENEVCF}.gz"

if [ ! -s "${CLINVARGENEVCF}" ]; then
    echo "File ${CLINVARGENEVCF} not found, extracting from full ClinVar file using coordinates..."
    # Ensure the master ClinVar file is indexed for coordinate slicing
    if [ ! -f "${DATADIR}/${CLINVARVCFGZ}.tbi" ] && [ ! -f "${DATADIR}/${CLINVARVCFGZ}.csi" ]; then
        echo "Master ClinVar index missing. Generating tabix index..."
        tabix -p vcf "${DATADIR}/${CLINVARVCFGZ}"
    fi
    # --- THE CHROMOSOME TRAP FIX ---
    # MANE GFF files use UCSC format ("chr5"), but NCBI ClinVar often uses plain integers ("5").
    # This snippet checks the master VCF format and cleans up your REGION string if needed.
    if bcftools view -h "${DATADIR}/${CLINVARVCFGZ}" | grep -q "ID=chr"; then
        CLINVAR_REGION="${REGION}"
    else
        CLINVAR_REGION="${REGION#chr}" # Strips the "chr" prefix if ClinVar expects plain numbers
    fi
    # Extract the target coordinates and isolate Pathogenic variants natively
    # Using 'CLNSIG ~ "Pathogenic"' captures "Pathogenic", "Pathogenic/Likely_pathogenic", etc.
    bcftools view \
        -r "${CLINVAR_REGION}" \
        -i 'CLNSIG=="Pathogenic"' \
        -O v \
        "${DATADIR}/${CLINVARVCFGZ}" > "${CLINVARGENEVCF}"
else
    echo "File ${CLINVARGENEVCF} exists, skipping extraction..."
fi
if [[ ! -s "${CLINVARGENEVCF}" ]]; then
    echo "Error: File ${CLINVARGENEVCF} not found or empty after coordinate extraction."
    exit 1
fi

# shift coordinates to local frame and rename chromosome to hg38 nomenclature
if [ ! -s  "${CLINVARGENE}.local.vcf" ]; then
  echo "Creating ${CLINVARGENE}.local.vcf..."
  awk -v seq="${REGION}" -v start="${START}" '
  BEGIN{FS=OFS="\t"}

  # keep headers
  /^#/ {print; next}

  {
      $1 = "chr"$1
      $2 = $2 - start + 1
      print
  }
  ' "${CLINVARGENEVCF}" > "${CLINVARGENE}.local.vcf"
else
    echo "File ${CLINVARGENE}.local.vcf exists, skipping creation..."
fi
if [[ ! -s "${CLINVARGENE}.local.vcf" ]]; then
    echo "Error: File ${CLINVARGENE}.local.vcf not found or empty."
    exit 1
fi

if [ ! -s  "${CLINVARGENE}.local.vcf.gz" ]; then
  echo "Creating ${CLINVARGENE}.local.vcf.gz..."
  bgzip -c "${CLINVARGENE}.local.vcf" > "${CLINVARGENE}.local.vcf.gz"
  tabix -p vcf "${CLINVARGENE}.local.vcf.gz"
else
    echo "File ${CLINVARGENE}.local.vcf.gz exists, skipping creation..."
fi
if [[ ! -s "${CLINVARGENE}.local.vcf.gz" ]]; then
    echo "Error: File ${CLINVARGENE}.local.vcf.gz not found or empty."
    exit 1
fi

echo "...all done"
