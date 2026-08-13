#!/usr/bin/env bash

# generate_foldx_jobs.sh
#
# Creates and optionally submits one Slurm job per .mmcif file.
#
# Usage:
#   ./generate_foldx_jobs.sh /path/to/mmcif_dir
#
# Requirements:
#   - mmcif_to_foldx.py
#   - FoldX executable
#   - Slurm cluster

set -euo pipefail

MMCIF_DIR="${1:-.}"

# Path to your Python script
PIPELINE="/groups/umcg-gcc/tmp04/users/umcg-tniemeijer/ESM/MORPHA/complex_analysis_scripts/analyse_complex.py"

# Path to FoldX
FOLDX="/groups/umcg-gcc/tmp04/users/umcg-tniemeijer/ESM/foldx_results/foldx_20261231"

OUTPUT_DIR="/groups/umcg-gcc/tmp04/users/umcg-tniemeijer/ESM/foldx_results/results_ASC"

# Directory to store generated sbatch scripts
JOB_DIR="${MMCIF_DIR}/slurm_jobs"
mkdir -p "${JOB_DIR}"

for MMCIF in "${MMCIF_DIR}"/*.mmcif; do

    [ -e "$MMCIF" ] || continue

    BASENAME=$(basename "$MMCIF" .mmcif)

    JOB_SCRIPT="${JOB_DIR}/${BASENAME}.sbatch"

    cat > "$JOB_SCRIPT" << EOF
#!/bin/bash
#SBATCH --job-name=fx_${BASENAME}
#SBATCH --output=${MMCIF_DIR}/${BASENAME}.out
#SBATCH --error=${MMCIF_DIR}/${BASENAME}.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G

set -euo pipefail

module purge

# Load your environment here
module load Python/3.12.3-GCCcore-13.3.0 

source /groups/umcg-gcc/tmp04/users/umcg-tniemeijer/ESM/MORPHA/.venv/bin/activate

WORKDIR=/groups/umcg-gcc/tmp04/users/umcg-tniemeijer/ESM/foldx_results/work

python "${PIPELINE}" \
"${MMCIF}" \
--foldx-path "${FOLDX}" \
--output-dir "${OUTPUT_DIR}" \
--work-dir "\$WORKDIR"

# Copy any FoldX output files generated in the work directory
find "${MMCIF_DIR}" -type f -name "*.fxout" -exec cp -f {} "${OUTPUT_DIR}/" \;

EOF

    echo "Created: $JOB_SCRIPT"

    # Uncomment to automatically submit:
    # sbatch "$JOB_SCRIPT"

done

echo "Done."