#!/bin/bash
#SBATCH --job-name=fpop_annotation
#SBATCH --partition=batch
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=../logs/annotation_%j.log
#
# fPOP Epilepsy Family Study - Stage 04: Variant Annotation
# Stanford School of Medicine, Department of Genetics (Snyder Lab)
#
# Annotates phased, pedigree-aware variant calls using ANNOVAR against
# RefGene, ClinVar, gnomAD, and epilepsy-relevant gene panels, to
# prioritize candidate disease-associated variants distinguishing the
# proband from unaffected family members.
#
# NOTE: No sample identifiers or patient data are included.

set -euo pipefail

ANNOVAR_DIR="${ANNOVAR_DIR:-../tools/annovar}"
HUMANDB_DIR="${HUMANDB_DIR:-../tools/annovar/humandb}"
INPUT_VCF="${INPUT_VCF:-../03_family_joint_phasing/joint_calls/family_phased_flagged.vcf.gz}"
OUTPUT_DIR="${OUTPUT_DIR:-./annotated}"

mkdir -p "${OUTPUT_DIR}"

echo "[fPOP:annotation] Converting to ANNOVAR input format..."
"${ANNOVAR_DIR}/convert2annovar.pl" \
  -format vcf4 "${INPUT_VCF}" \
  -outfile "${OUTPUT_DIR}/family_variants.avinput" \
  -allsample -withfreq -includeinfo

echo "[fPOP:annotation] Running table_annovar..."
"${ANNOVAR_DIR}/table_annovar.pl" \
  "${OUTPUT_DIR}/family_variants.avinput" \
  "${HUMANDB_DIR}" \
  -buildver hg38 \
  -out "${OUTPUT_DIR}/family_variants_annotated" \
  -remove \
  -protocol refGene,clinvar_20240611,gnomad211_exome,dbnsfp42a \
  -operation g,f,f,f \
  -nastring . \
  -csvout

echo "[fPOP:annotation] Filtering for rare, predicted-damaging variants..."
# Rare (gnomAD AF < 0.001) and predicted deleterious by in-silico tools
python3 filter_candidate_variants.py \
  --input "${OUTPUT_DIR}/family_variants_annotated.hg38_multianno.csv" \
  --output "${OUTPUT_DIR}/candidate_disease_genes.csv" \
  --max-af 0.001

echo "[fPOP:annotation] Complete. Candidate gene list ready for cross-omic integration."
