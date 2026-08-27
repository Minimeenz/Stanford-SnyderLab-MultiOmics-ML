#!/bin/bash
#SBATCH --job-name=fpop_atacseq
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=../logs/atacseq_%j.log
#
# fPOP Epilepsy Family Study — ATAC-seq (Chromatin Accessibility)
# Stanford School of Medicine, Department of Genetics (Snyder Lab)
#
# Processes ATAC-seq data across the proband and family members to assess
# whether chromatin accessibility differs at regulatory regions near
# candidate genes identified in the genomic (04) and methylation (06)
# analyses — testing whether variants/DMRs correspond to altered
# regulatory activity.
#
# NOTE: No sample identifiers or patient data are included.

set -euo pipefail

REFERENCE="${REFERENCE:-../references/GRCh38.fa}"
FASTQ_DIR="${FASTQ_DIR:-./raw_fastq}"
OUTPUT_DIR="${OUTPUT_DIR:-./results}"
MANIFEST="${MANIFEST:-../metadata/atacseq_sample_manifest.tsv}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

mkdir -p "${OUTPUT_DIR}/peaks" "${OUTPUT_DIR}/bams"

while IFS=$'\t' read -r sample_id fastq_r1 fastq_r2; do
  echo "[fPOP:atacseq] Processing ${sample_id}..."

  # Alignment with Bowtie2 (standard for ATAC-seq)
  bowtie2 -p "${THREADS}" -X 2000 --very-sensitive \
    -x "${REFERENCE%.fa}" \
    -1 "${FASTQ_DIR}/${fastq_r1}" -2 "${FASTQ_DIR}/${fastq_r2}" \
    | samtools sort -@ "${THREADS}" -o "${OUTPUT_DIR}/bams/${sample_id}.sorted.bam" -

  samtools index "${OUTPUT_DIR}/bams/${sample_id}.sorted.bam"

  # Remove mitochondrial reads and duplicates (standard ATAC-seq QC)
  samtools view -@ "${THREADS}" -b -F 1024 \
    "${OUTPUT_DIR}/bams/${sample_id}.sorted.bam" chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 \
    chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY \
    > "${OUTPUT_DIR}/bams/${sample_id}.filtered.bam"

  # Peak calling with MACS2
  macs2 callpeak \
    -t "${OUTPUT_DIR}/bams/${sample_id}.filtered.bam" \
    -f BAMPE -g hs \
    -n "${sample_id}" \
    --outdir "${OUTPUT_DIR}/peaks" \
    -q 0.01 --nomodel --shift -100 --extsize 200

  echo "[fPOP:atacseq] Completed ${sample_id}"
done < "${MANIFEST}"

echo "[fPOP:atacseq] Generating consensus peak set and differential accessibility (proband vs. family)..."
python3 differential_accessibility.py \
  --peaks-dir "${OUTPUT_DIR}/peaks" \
  --groups ../metadata/proband_vs_family_groups.csv \
  --output "${OUTPUT_DIR}/differential_accessibility.csv"

echo "[fPOP:atacseq] Complete."
