#!/bin/bash

[ -z "$1" ] && { echo "Missing argument: gene name"; exit 1; }

GENE=$1

DATADIR=data
TMPDIR=tmp
GENEDIR="${TMPDIR}/${GENE}"
mkdir -p "${DATADIR}"
mkdir -p "${TMPDIR}"
mkdir -p "${GENEDIR}"

CLINVARVCFGZ=clinvar_20260503.vcf.gz

if [ ! -s "${DATADIR}/${CLINVARVCFGZ}" ]; then
    echo "File ${DATADIR}/${CLINVARVCFGZ} not found, downloading..."
    wget -P "${DATADIR}" "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2026/${GRCh38REFFA}"
else
    echo "File ${DATADIR}/${CLINVARVCFGZ} exists, skipping download..."
fi
if [[ ! -s "${DATADIR}/${CLINVARVCFGZ}" ]]; then
    echo "Error: File ${DATADIR}/${CLINVARVCFGZ} not found or empty."
    exit 1
fi

CLINVARGENE="${GENEDIR}/ClinVarPathogenic"
CLINVARGENEVCF="${CLINVARGENE}.vcf"
CLINVARGENEVCFGZ="${CLINVARGENEVCF}.gz"

if [ ! -s "${CLINVARGENEVCF}" ]; then
    echo "File ${CLINVARGENEVCF} not found, extracting from full ClinVar file..."
    printf "##fileformat=VCFv4.2\n" > "${CLINVARGENEVCF}"
    printf "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" >> "${CLINVARGENEVCF}"
    zgrep -i "${GENE}" "${DATADIR}/${CLINVARVCFGZ}" | grep -i "CLNSIG=Pathogenic;" >> "${CLINVARGENEVCF}"
else
    echo "File ${CLINVARGENEVCF} exists, skipping extraction..."
fi
if [[ ! -s "${CLINVARGENEVCF}" ]]; then
    echo "Error: File ${CLINVARGENEVCF} not found or empty."
    exit 1
fi

# get region and replace in VCF
GRCh38MANEGFF=MANE.GRCh38.v1.3.ensembl_genomic.gff
DATA=$(awk -F'\t' '$3=="transcript" {print $1, $4, $5; exit}' "${GENEDIR}/${GRCh38MANEGFF}")
if [ -z "$DATA" ]; then
    echo "Error: Could not find a 'transcript' feature in ${GENEDIR}/${GRCh38MANEGFF}."
    exit 1
fi
read CHR START END <<< "$DATA"
if [[ -z "$CHR" || -z "$START" || -z "$END" ]]; then
    echo "Error: Coordinate extraction failed (one or more fields are empty)."
    exit 1
fi
REGION="${CHR}:${START}-${END}"
echo "Extracted region for this gene is $REGION"

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

exit 0

# create consensus based on previously extracted samples (separate script?)

bcftools consensus \
    -f tmp/MEFV/samples/HG00620.h1.fa \
    tmp/MEFV/ClinVarPathogenic.vcf.gz \
    > CVtest_MEFV.HG00620.h1.fa

gffread lifted.local.gff \
    -g MEFV.HG00620.h1.clinvar449657.fa \
    -x tmp.cds.fa \
    -y tmp.protein.fa \
    -S
