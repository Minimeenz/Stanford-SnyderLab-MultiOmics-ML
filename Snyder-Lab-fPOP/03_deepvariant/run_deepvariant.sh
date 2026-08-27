#!/bin/bash
#SBATCH --job-name=fpop_deepvariant
#SBATCH --partition=batch
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --gres=gpu:1
#SBATCH --time=12:00:00
#SBATCH --output=../logs/deepvariant_%j.log
#
# Snyder Lab fPOP - DeepVariant Deep-Learning Variant Calling
#
# Runs Google's DeepVariant as an orthogonal, deep-learning-based variant
# caller alongside the GATK-based pipeline (02_variant_calling), providing
# a second independent variant call set for concordance filtering - variants
# supported by both callers are prioritized as high-confidence candidates.
#
# NOTE: No sample identifiers or patient data are included.

set -euo pipefail

REFERENCE="${REFERENCE:-../references/GRCh38.fa}"
INPUT_DIR="${INPUT_DIR:-../01_mapping/aligned_bams}"
OUTPUT_DIR="${OUTPUT_DIR:-./deepvariant_calls}"
MODEL_TYPE="${MODEL_TYPE:-WGS}"  # WGS, WES, or PACBIO depending on input data
THREADS="${SLURM_CPUS_PER_TASK:-16}"

mkdir -p "${OUTPUT_DIR}"

for bam_file in "${INPUT_DIR}"/*.final.bam; do
  sample_id=$(basename "${bam_file}" .final.bam)
  echo "[fPOP:deepvariant] Calling variants for ${sample_id} (model: ${MODEL_TYPE})..."

  singularity run --nv docker://google/deepvariant:latest \
    /opt/deepvariant/bin/run_deepvariant \
    --model_type="${MODEL_TYPE}" \
    --ref="${REFERENCE}" \
    --reads="${bam_file}" \
    --output_vcf="${OUTPUT_DIR}/${sample_id}.deepvariant.vcf.gz" \
    --output_gvcf="${OUTPUT_DIR}/${sample_id}.deepvariant.g.vcf.gz" \
    --num_shards="${THREADS}"

  echo "[fPOP:deepvariant] Completed ${sample_id}"
done

echo "[fPOP:deepvariant] Computing concordance with GATK call set..."
python3 concordance_check.py \
  --deepvariant-dir "${OUTPUT_DIR}" \
  --gatk-dir ../02_variant_calling/gvcfs \
  --output "${OUTPUT_DIR}/caller_concordance_summary.csv"

echo "[fPOP:deepvariant] Complete."
