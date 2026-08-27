#!/bin/bash
#
# fPOP Epilepsy Family Study - GTF Chromosome Naming Conversion Utility
# Stanford School of Medicine, Department of Genetics (Snyder Lab)
#
# Converts chromosome naming convention in a GTF annotation file between
# Ensembl-style ("1", "2", ..., "X") and UCSC-style ("chr1", "chr2", ...,
# "chrX"), needed to keep annotation files consistent across tools that
# expect different reference genome naming conventions (common friction
# point when combining GATK/UCSC-based tools with Ensembl-based ones).
#
# Usage: ./convert_gtf_chr_naming.sh <input.gtf> <output.gtf> <to_ucsc|to_ensembl>

set -euo pipefail

INPUT_GTF="${1:?Usage: $0 <input.gtf> <output.gtf> <to_ucsc|to_ensembl>}"
OUTPUT_GTF="${2:?Usage: $0 <input.gtf> <output.gtf> <to_ucsc|to_ensembl>}"
DIRECTION="${3:?Usage: $0 <input.gtf> <output.gtf> <to_ucsc|to_ensembl>}"

echo "[fPOP] Converting ${INPUT_GTF} (${DIRECTION})..."

if [[ "${DIRECTION}" == "to_ucsc" ]]; then
  # Prepend "chr" to chromosome names that don't already have it
  awk 'BEGIN{OFS="\t"} /^#/{print; next} {
    if ($1 !~ /^chr/) { $1 = "chr" $1 }
    print
  }' "${INPUT_GTF}" > "${OUTPUT_GTF}"

elif [[ "${DIRECTION}" == "to_ensembl" ]]; then
  # Strip leading "chr" from chromosome names
  awk 'BEGIN{OFS="\t"} /^#/{print; next} {
    sub(/^chr/, "", $1)
    print
  }' "${INPUT_GTF}" > "${OUTPUT_GTF}"

else
  echo "[fPOP] ERROR: direction must be 'to_ucsc' or 'to_ensembl'" >&2
  exit 1
fi

echo "[fPOP] Converted GTF written to ${OUTPUT_GTF}"
