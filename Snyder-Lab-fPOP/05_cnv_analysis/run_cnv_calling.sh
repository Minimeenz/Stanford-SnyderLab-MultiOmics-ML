#!/bin/bash
#SBATCH --job-name=fpop_cnv
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=../logs/cnv_%j.log
#
# Snyder Lab fPOP — Copy Number Variant (CNV) Analysis
#
# Detects copy number variants (deletions/duplications) using GATK's
# germline CNV calling pipeline (gCNV), an important complementary
# variant class to SNVs/indels for epilepsy gene discovery, since many
# epilepsy-associated genes have documented pathogenic CNVs.
#
# NOTE: No sample identifiers or patient data are included.

set -euo pipefail

REFERENCE="${REFERENCE:-../references/GRCh38.fa}"
INPUT_DIR="${INPUT_DIR:-../01_mapping/aligned_bams}"
INTERVALS="${INTERVALS:-../references/exome_targets.interval_list}"
OUTPUT_DIR="${OUTPUT_DIR:-./cnv_calls}"

mkdir -p "${OUTPUT_DIR}"

echo "[fPOP:cnv] Collecting read counts per sample..."
for bam_file in "${INPUT_DIR}"/*.final.bam; do
  sample_id=$(basename "${bam_file}" .final.bam)

  gatk CollectReadCounts \
    -I "${bam_file}" \
    -L "${INTERVALS}" \
    --interval-merging-rule OVERLAPPING_ONLY \
    -O "${OUTPUT_DIR}/${sample_id}.counts.hdf5"
done

echo "[fPOP:cnv] Building panel of normals from unaffected family members..."
gatk CreateReadCountPanelOfNormals \
  -I "${OUTPUT_DIR}"/*family*.counts.hdf5 \
  -O "${OUTPUT_DIR}/family_pon.hdf5"

echo "[fPOP:cnv] Denoising proband read counts against panel of normals..."
gatk DenoiseReadCounts \
  -I "${OUTPUT_DIR}"/*proband*.counts.hdf5 \
  --count-panel-of-normals "${OUTPUT_DIR}/family_pon.hdf5" \
  --standardized-copy-ratios "${OUTPUT_DIR}/proband.standardizedCR.tsv" \
  --denoised-copy-ratios "${OUTPUT_DIR}/proband.denoisedCR.tsv"

echo "[fPOP:cnv] Segmenting copy ratios..."
gatk ModelSegments \
  --denoised-copy-ratios "${OUTPUT_DIR}/proband.denoisedCR.tsv" \
  --output "${OUTPUT_DIR}" \
  --output-prefix proband

echo "[fPOP:cnv] Calling copy-neutral / amplified / deleted segments..."
gatk CallCopyRatioSegments \
  --input "${OUTPUT_DIR}/proband.cr.seg" \
  --output "${OUTPUT_DIR}/proband.called.seg"

echo "[fPOP:cnv] Complete. Calls at ${OUTPUT_DIR}/proband.called.seg"
