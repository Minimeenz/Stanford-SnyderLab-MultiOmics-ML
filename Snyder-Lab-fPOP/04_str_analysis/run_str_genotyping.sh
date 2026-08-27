#!/bin/bash
#SBATCH --job-name=fpop_str_analysis
#SBATCH --partition=batch
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=../logs/str_analysis_%j.log
#
# Snyder Lab fPOP — Short Tandem Repeat (STR) / Repeat Expansion Analysis
#
# Genotypes known disease-associated repeat expansion loci using
# ExpansionHunter, relevant given known STR-associated epilepsy and
# neurological disease loci (e.g., in genes such as ATXN2, ATXN1, FMR1,
# and other repeat-associated neurological disease genes).
#
# NOTE: No sample identifiers or patient data are included.

set -euo pipefail

REFERENCE="${REFERENCE:-../references/GRCh38.fa}"
INPUT_DIR="${INPUT_DIR:-../01_mapping/aligned_bams}"
OUTPUT_DIR="${OUTPUT_DIR:-./str_calls}"
VARIANT_CATALOG="${VARIANT_CATALOG:-../references/expansion_hunter_variant_catalog_hg38.json}"

mkdir -p "${OUTPUT_DIR}"

for bam_file in "${INPUT_DIR}"/*.final.bam; do
  sample_id=$(basename "${bam_file}" .final.bam)
  echo "[fPOP:str] Genotyping repeat expansion loci for ${sample_id}..."

  ExpansionHunter \
    --reads "${bam_file}" \
    --reference "${REFERENCE}" \
    --variant-catalog "${VARIANT_CATALOG}" \
    --output-prefix "${OUTPUT_DIR}/${sample_id}"

  echo "[fPOP:str] Completed ${sample_id}"
done

echo "[fPOP:str] Aggregating repeat sizes across family and flagging pathogenic-range expansions..."
python3 flag_pathogenic_expansions.py \
  --str-dir "${OUTPUT_DIR}" \
  --groups ../metadata/proband_vs_family_groups.csv \
  --output "${OUTPUT_DIR}/str_summary.csv"

echo "[fPOP:str] Complete."
