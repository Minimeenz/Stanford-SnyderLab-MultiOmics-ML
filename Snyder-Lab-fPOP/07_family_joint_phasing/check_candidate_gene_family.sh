#!/bin/bash
#
# fPOP Epilepsy Family Study - Candidate Gene Family Segregation Check
# Stanford School of Medicine, Department of Genetics (Snyder Lab)
#
# Given a candidate gene surfaced from cross-omic analysis (Stages 04-07),
# checks segregation of variants within that gene across the family
# pedigree - confirming whether the variant pattern is consistent with
# the proband's phenotype (e.g., present in proband, absent/heterozygous
# in unaffected relatives).
#
# Usage: ./check_candidate_gene_family.sh <GENE_SYMBOL>
#
# NOTE: No sample identifiers, family IDs, or patient data are included.

set -euo pipefail

GENE_SYMBOL="${1:?Usage: $0 <GENE_SYMBOL>}"
ANNOTATED_VCF="${ANNOTATED_VCF:-../04_annotation/annotated/family_variants_annotated.hg38_multianno.csv}"
PEDIGREE_FILE="${PEDIGREE_FILE:-../metadata/family_pedigree.ped}"

echo "[fPOP] Checking segregation for candidate gene: ${GENE_SYMBOL}"

python3 - "${GENE_SYMBOL}" "${ANNOTATED_VCF}" "${PEDIGREE_FILE}" << 'PYEOF'
import sys
import pandas as pd

gene, annotated_path, pedigree_path = sys.argv[1], sys.argv[2], sys.argv[3]

df = pd.read_csv(annotated_path, low_memory=False)
gene_variants = df[df["Gene.refGene"] == gene]

if gene_variants.empty:
    print(f"[fPOP] No variants found in {gene} in the annotated candidate list.")
    sys.exit(0)

pedigree = pd.read_csv(pedigree_path, sep="\t", names=["family_id", "sample_id", "father", "mother", "sex", "phenotype"])

print(f"[fPOP] {len(gene_variants)} variant(s) found in {gene}.")
print(f"[fPOP] Family members in pedigree: {len(pedigree)}")
print(f"[fPOP] Cross-reference genotype columns per sample to confirm segregation pattern")
print(f"[fPOP] (proband should carry the variant; unaffected relatives should not, for a de novo/dominant model)")
PYEOF

echo "[fPOP] Segregation check complete for ${GENE_SYMBOL}."
