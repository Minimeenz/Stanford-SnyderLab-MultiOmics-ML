#!/bin/bash
#SBATCH --job-name=fpop_joint_phasing
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=../logs/joint_phasing_%j.log
#
# fPOP Epilepsy Family Study — Stage 03: Family Joint Genotyping & Phasing
# Stanford School of Medicine, Department of Genetics (Snyder Lab)
#
# Jointly genotypes the proband and unaffected family members, then
# performs pedigree-aware phasing to determine variant inheritance
# patterns (de novo, inherited, compound heterozygous) — critical for
# distinguishing candidate causal variants in the proband from background
# family variation.
#
# NOTE: No sample identifiers, family IDs, or patient data are included.
# Pedigree file (PED format) is expected to contain de-identified IDs only.

set -euo pipefail

REFERENCE="${REFERENCE:-../references/GRCh38.fa}"
GVCF_DIR="${GVCF_DIR:-../02_variant_calling/gvcfs}"
PEDIGREE_FILE="${PEDIGREE_FILE:-../metadata/family_pedigree.ped}"
OUTPUT_DIR="${OUTPUT_DIR:-./joint_calls}"

mkdir -p "${OUTPUT_DIR}"

echo "[fPOP:joint_phasing] Combining GVCFs across family members..."
gvcf_args=()
for gvcf in "${GVCF_DIR}"/*.g.vcf.gz; do
  gvcf_args+=(-V "${gvcf}")
done

gatk CombineGVCFs \
  -R "${REFERENCE}" \
  "${gvcf_args[@]}" \
  -O "${OUTPUT_DIR}/family_combined.g.vcf.gz"

echo "[fPOP:joint_phasing] Joint genotyping across family..."
gatk GenotypeGVCFs \
  -R "${REFERENCE}" \
  -V "${OUTPUT_DIR}/family_combined.g.vcf.gz" \
  -O "${OUTPUT_DIR}/family_joint.vcf.gz"

echo "[fPOP:joint_phasing] Running pedigree-aware phasing..."
# Determines inheritance pattern per variant (de novo / inherited / compound het)
# using family structure defined in the pedigree file.
gatk CalculateGenotypePosteriors \
  -V "${OUTPUT_DIR}/family_joint.vcf.gz" \
  -ped "${PEDIGREE_FILE}" \
  -O "${OUTPUT_DIR}/family_phased.vcf.gz"

echo "[fPOP:joint_phasing] Flagging de novo candidates in proband..."
gatk VariantFiltration \
  -V "${OUTPUT_DIR}/family_phased.vcf.gz" \
  --genotype-filter-expression "isHet == 1 && GQ < 20" \
  --genotype-filter-name "lowGQdenovo" \
  -O "${OUTPUT_DIR}/family_phased_flagged.vcf.gz"

echo "[fPOP:joint_phasing] Complete. Phased, pedigree-aware calls ready for annotation."
