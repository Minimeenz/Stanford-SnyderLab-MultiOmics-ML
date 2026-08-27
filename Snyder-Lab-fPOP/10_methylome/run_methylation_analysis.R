#!/usr/bin/env Rscript
#
# fPOP Epilepsy Family Study - Stage 06: DNA Methylation Analysis
# Stanford School of Medicine, Department of Genetics (Snyder Lab)
#
# Processes whole-genome bisulfite sequencing (or array-based) methylation
# data across the proband and family members, identifies differentially
# methylated regions (DMRs), and cross-references DMRs against candidate
# genes from genomic (Stage 04) and proteomic (Stage 05) analyses.
#
# NOTE: No sample identifiers or patient data are included.

suppressPackageStartupMessages({
  library(methylKit)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
methylation_dir <- ifelse(length(args) >= 1, args[1], "./methylation_calls")
sample_manifest  <- ifelse(length(args) >= 2, args[2], "../metadata/proband_vs_family_groups.csv")
output_dir       <- ifelse(length(args) >= 3, args[3], "./results")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("[fPOP:methylation] Loading sample manifest...\n")
manifest <- read.csv(sample_manifest)
sample_ids <- as.list(manifest$sample_id)
treatment  <- ifelse(manifest$group == "proband", 1, 0)

cat("[fPOP:methylation] Reading per-CpG methylation calls...\n")
cpg_files <- as.list(file.path(methylation_dir, paste0(manifest$sample_id, ".CpG_report.txt.gz")))

meth_obj <- methRead(
  cpg_files,
  sample.id = sample_ids,
  assembly  = "hg38",
  treatment = treatment,
  context   = "CpG",
  mincov    = 10
)

cat("[fPOP:methylation] Filtering and normalizing coverage...\n")
filtered <- filterByCoverage(meth_obj, lo.count = 10, lo.perc = NULL, hi.perc = 99.9)
normalized <- normalizeCoverage(filtered)

cat("[fPOP:methylation] Merging samples at common CpG sites...\n")
merged <- unite(normalized, destrand = FALSE)

cat("[fPOP:methylation] Calling differentially methylated CpGs (proband vs. family)...\n")
diff_meth <- calculateDiffMeth(merged, mc.cores = 4)
dmrs <- getMethylDiff(diff_meth, difference = 25, qvalue = 0.05)

cat("[fPOP:methylation] Annotating DMRs to nearest genes...\n")
dmr_df <- getData(dmrs) %>%
  arrange(qvalue)

write.csv(dmr_df, file.path(output_dir, "differentially_methylated_regions.csv"), row.names = FALSE)

cat(sprintf(
  "[fPOP:methylation] Complete. %d significant DMRs identified (>25%% diff, q < 0.05).\n",
  nrow(dmr_df)
))
