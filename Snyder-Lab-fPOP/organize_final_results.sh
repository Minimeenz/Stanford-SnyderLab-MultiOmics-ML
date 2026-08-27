#!/bin/bash
#
# fPOP Epilepsy Family Study - Final Results Organization
# Stanford School of Medicine, Department of Genetics (Snyder Lab)
#
# Aggregates key outputs from each omic pipeline stage (genomic
# candidates, differentially abundant proteins, DMRs, differential
# cytokines, differential accessibility, viral screening, multiome DE)
# into a single organized final_results directory for cross-omic
# integration and review.
#
# NOTE: No sample identifiers or patient data are included.

set -euo pipefail

FINAL_RESULTS_DIR="${FINAL_RESULTS_DIR:-../final_results}"

mkdir -p "${FINAL_RESULTS_DIR}"/{genomic,proteomic,epigenomic,immune,chromatin,virome,multiome}

echo "[fPOP] Collecting genomic candidates..."
cp ../04_annotation/annotated/candidate_disease_genes.csv \
   "${FINAL_RESULTS_DIR}/genomic/candidate_disease_genes.csv"

echo "[fPOP] Collecting proteomic results..."
cp ../05_olink_proteomics/results/olink_differential_abundance.csv \
   "${FINAL_RESULTS_DIR}/proteomic/olink_differential_abundance.csv"

echo "[fPOP] Collecting methylation results..."
cp ../06_methylation_analysis/results/differentially_methylated_regions.csv \
   "${FINAL_RESULTS_DIR}/epigenomic/differentially_methylated_regions.csv"

echo "[fPOP] Collecting cytokine results..."
cp ../07_cytokine_analysis/results/cytokine_differential_analysis.csv \
   "${FINAL_RESULTS_DIR}/immune/cytokine_differential_analysis.csv"

echo "[fPOP] Collecting ATAC-seq results..."
cp ../fPOP_Epi_ATACSeq/results/differential_accessibility.csv \
   "${FINAL_RESULTS_DIR}/chromatin/differential_accessibility.csv"

echo "[fPOP] Collecting viral screening results..."
cp ../CSF96_Viral_enrichment/results/viral_abundance_summary.csv \
   "${FINAL_RESULTS_DIR}/virome/viral_abundance_summary.csv"

echo "[fPOP] Collecting multiome results..."
cp ../multiome_analysis/results/cluster_differential_expression.csv \
   "${FINAL_RESULTS_DIR}/multiome/cluster_differential_expression.csv"

echo "[fPOP] Final results organized at ${FINAL_RESULTS_DIR}"
echo "[fPOP] Ready for cross-omic integration and candidate gene prioritization."
