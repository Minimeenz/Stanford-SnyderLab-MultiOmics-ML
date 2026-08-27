#!/bin/bash
#SBATCH --job-name=fpop_mapping
#SBATCH --partition=batch
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=../logs/mapping_%j.log
#
# fPOP Epilepsy Family Study - Stage 01: Read Mapping
# Stanford School of Medicine, Department of Genetics (Snyder Lab)
#
# Aligns raw sequencing reads (proband + family members) to the reference
# genome using BWA-MEM2, followed by sorting, duplicate marking, and BQSR.
# Runs on Stanford's SCG compute cluster via SLURM.
#
# NOTE: No sample identifiers, family IDs, or patient data are included.
# Sample list is read from a de-identified manifest at runtime.

set -euo pipefail

REFERENCE="${REFERENCE:-../references/GRCh38.fa}"
FASTQ_DIR="${FASTQ_DIR:-./raw_fastq}"
OUTPUT_DIR="${OUTPUT_DIR:-./aligned_bams}"
MANIFEST="${MANIFEST:-../metadata/sample_manifest.tsv}"  # de-identified sample IDs only
THREADS="${SLURM_CPUS_PER_TASK:-16}"

mkdir -p "${OUTPUT_DIR}"

while IFS=$'\t' read -r sample_id fastq_r1 fastq_r2; do
  echo "[fPOP:mapping] Aligning ${sample_id}..."

  bwa-mem2 mem -t "${THREADS}" -R "@RG\tID:${sample_id}\tSM:${sample_id}\tPL:ILLUMINA" \
    "${REFERENCE}" "${FASTQ_DIR}/${fastq_r1}" "${FASTQ_DIR}/${fastq_r2}" \
    | samtools sort -@ "${THREADS}" -o "${OUTPUT_DIR}/${sample_id}.sorted.bam" -

  samtools index "${OUTPUT_DIR}/${sample_id}.sorted.bam"

  # Mark duplicates
  gatk MarkDuplicates \
    -I "${OUTPUT_DIR}/${sample_id}.sorted.bam" \
    -O "${OUTPUT_DIR}/${sample_id}.dedup.bam" \
    -M "${OUTPUT_DIR}/${sample_id}.dup_metrics.txt"

  # Base quality score recalibration
  gatk BaseRecalibrator \
    -I "${OUTPUT_DIR}/${sample_id}.dedup.bam" \
    -R "${REFERENCE}" \
    --known-sites ../references/known_sites.vcf.gz \
    -O "${OUTPUT_DIR}/${sample_id}.recal.table"

  gatk ApplyBQSR \
    -I "${OUTPUT_DIR}/${sample_id}.dedup.bam" \
    -R "${REFERENCE}" \
    --bqsr-recal-file "${OUTPUT_DIR}/${sample_id}.recal.table" \
    -O "${OUTPUT_DIR}/${sample_id}.final.bam"

  echo "[fPOP:mapping] Completed ${sample_id}"
done < "${MANIFEST}"

echo "[fPOP:mapping] All samples mapped."
