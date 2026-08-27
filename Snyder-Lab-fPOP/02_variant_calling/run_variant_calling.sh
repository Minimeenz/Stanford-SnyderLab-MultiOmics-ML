#!/bin/bash
#SBATCH --job-name=fpop_variant_calling
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=../logs/variant_calling_%j.log
#
# fPOP Epilepsy Family Study — Stage 02: Variant Calling
# Stanford School of Medicine, Department of Genetics (Snyder Lab)
#
# Calls per-sample variants from recalibrated BAMs using GATK
# HaplotypeCaller in GVCF mode, ahead of joint genotyping across the
# family in Stage 03. Runs on Stanford's SCG compute cluster.
#
# NOTE: No sample identifiers or patient data are included.

set -euo pipefail

REFERENCE="${REFERENCE:-../references/GRCh38.fa}"
INPUT_DIR="${INPUT_DIR:-../01_mapping/aligned_bams}"
OUTPUT_DIR="${OUTPUT_DIR:-./gvcfs}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

mkdir -p "${OUTPUT_DIR}"

for bam_file in "${INPUT_DIR}"/*.final.bam; do
  sample_id=$(basename "${bam_file}" .final.bam)
  echo "[fPOP:variant_calling] Calling variants for ${sample_id}..."

  gatk HaplotypeCaller \
    -R "${REFERENCE}" \
    -I "${bam_file}" \
    -O "${OUTPUT_DIR}/${sample_id}.g.vcf.gz" \
    -ERC GVCF \
    --native-pair-hmm-threads "${THREADS}"

  echo "[fPOP:variant_calling] Completed ${sample_id}"
done

echo "[fPOP:variant_calling] All per-sample GVCFs generated."
echo "[fPOP:variant_calling] Proceed to 03_family_joint_phasing for joint genotyping."
