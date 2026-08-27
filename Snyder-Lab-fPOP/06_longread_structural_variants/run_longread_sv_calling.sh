#!/bin/bash
#SBATCH --job-name=fpop_longread_sv
#SBATCH --partition=batch
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=16:00:00
#SBATCH --output=../logs/longread_sv_%j.log
#
# Snyder Lab fPOP — Long-Read Structural Variant Calling
#
# Aligns long-read sequencing data (PacBio HiFi / ONT) using minimap2 and
# calls structural variants (large insertions, deletions, inversions,
# translocations) with Sniffles2. Long-read SV calling complements the
# short-read CNV pipeline (05) with improved sensitivity for complex and
# repeat-region structural variation relevant to ASD/epilepsy etiology.
#
# NOTE: No sample identifiers or patient data are included.

set -euo pipefail

REFERENCE="${REFERENCE:-../references/GRCh38.fa}"
LONGREAD_DIR="${LONGREAD_DIR:-./raw_longreads}"
OUTPUT_DIR="${OUTPUT_DIR:-./sv_calls}"
PRESET="${PRESET:-map-hifi}"  # map-hifi for PacBio HiFi, map-ont for Oxford Nanopore
THREADS="${SLURM_CPUS_PER_TASK:-16}"

mkdir -p "${OUTPUT_DIR}/bams"

for longread_fastq in "${LONGREAD_DIR}"/*.fastq.gz; do
  sample_id=$(basename "${longread_fastq}" .fastq.gz)
  echo "[fPOP:longread_sv] Aligning long reads for ${sample_id} (preset: ${PRESET})..."

  minimap2 -ax "${PRESET}" -t "${THREADS}" \
    "${REFERENCE}" "${longread_fastq}" \
    | samtools sort -@ "${THREADS}" -o "${OUTPUT_DIR}/bams/${sample_id}.sorted.bam" -

  samtools index "${OUTPUT_DIR}/bams/${sample_id}.sorted.bam"

  echo "[fPOP:longread_sv] Calling structural variants for ${sample_id}..."
  sniffles \
    --input "${OUTPUT_DIR}/bams/${sample_id}.sorted.bam" \
    --vcf "${OUTPUT_DIR}/${sample_id}.sniffles.vcf.gz" \
    --reference "${REFERENCE}" \
    --threads "${THREADS}"

  echo "[fPOP:longread_sv] Completed ${sample_id}"
done

echo "[fPOP:longread_sv] Merging family SV calls for joint analysis..."
sniffles \
  --input "${OUTPUT_DIR}/bams"/*.sorted.bam \
  --vcf "${OUTPUT_DIR}/family_merged.sniffles.vcf.gz" \
  --reference "${REFERENCE}" \
  --threads "${THREADS}"

echo "[fPOP:longread_sv] Complete."
