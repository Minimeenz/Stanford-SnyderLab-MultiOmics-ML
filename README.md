# Snyder Lab fPOP - Familial Personal Omics Profiling: Epilepsy Family Study

Computational pipeline for **fPOP** (familial Personal Omics Profiling), directed by Dr. Michael Snyder at Stanford School of Medicine's Department of Genetics (Snyder Lab). This study extends the Snyder Lab's iPOP (integrated Personal Omics Profiling) approach to a **family-based design**, applied here to an epilepsy family case: a **proband** (affected individual) compared against **unaffected family members** across a comprehensive multi-omic panel, to identify candidate disease-associated molecular signatures.

> **Privacy note:** This repository contains pipeline code and directory structure only. No sample identifiers, family IDs, patient data, or PHI is in this repo.

## Study Design

Proband vs. family comparison across the following omics layers:

| # | Omics Layer | Folder | Method |
|---|---|---|---|
| 1 | **Read mapping** | `01_mapping` | BWA-MEM2 alignment, dedup, BQSR |
| 2 | **Variant calling** | `02_variant_calling` | GATK HaplotypeCaller (GVCF mode) |
| 3 | **Deep-learning variant calling** | `03_deepvariant` | Google DeepVariant, cross-caller concordance vs. GATK |
| 4 | **Short tandem repeats (STR)** | `04_str_analysis` | ExpansionHunter, pathogenic-range repeat expansion flagging |
| 5 | **Copy number variants (CNV)** | `05_cnv_analysis` | GATK gCNV (read-depth based, family panel of normals) |
| 6 | **Long-read structural variants** | `06_longread_structural_variants` | minimap2 + Sniffles2 (PacBio HiFi / ONT) |
| 7 | **Family joint phasing** | `07_family_joint_phasing` | Pedigree-aware joint genotyping, de novo/inherited variant classification, candidate-gene family segregation check |
| 8 | **Annotation** | `08_annotation` | ANNOVAR (RefGene, ClinVar, gnomAD, dbNSFP), rare/damaging variant filtering |
| 9 | **Proteomics** | `09_proteomics` | Olink NPX differential abundance |
| 10 | **Methylome** | `10_methylome` | Whole-genome bisulfite sequencing, DMR calling (methylKit) |
| 11 | **Cytokines / immune profiling** | `11_cytokines` | Multiplex cytokine panel, differential analyte abundance |
| 12 | **Chromatin accessibility** | `12_atac_seq` | ATAC-seq, MACS2 peak calling, differential accessibility |
| 13 | **Single-cell multiome** | `13_multiome_analysis` | Joint scRNA + scATAC, Scanpy clustering, cell-type-resolved differential expression |
| 14 | **Structural biology** | `14_structural_biology` | ESM-based variant effect prediction, AlphaFold structural context for candidate missense variants |
| 15 | **Spatial transcriptomics** | `15_spatial_transcriptomics` | 10x Visium, spatial autocorrelation (Moran's I) of candidate gene expression |
| 16 | **Metabolomics** | `16_metabolomics` | LC-MS/MS peak table analysis, differential metabolite abundance |
| 17 | **Cross-omic integration** | `17_cross_omic_integration` | Rank-normalized multi-omic feature fusion, gradient-boosted candidate gene ranking, per-gene interpretability reports |

## Pipeline Structure

```
Snyder-Lab-fPOP/
├── 01_mapping/
├── 02_variant_calling/
├── 03_deepvariant/
├── 04_str_analysis/
├── 05_cnv_analysis/
├── 06_longread_structural_variants/
├── 07_family_joint_phasing/            # includes candidate-gene segregation check
├── 08_annotation/                       # includes GTF chr-naming conversion utility
├── 09_proteomics/
├── 10_methylome/
├── 11_cytokines/
├── 12_atac_seq/
├── 13_multiome_analysis/
├── 14_structural_biology/
├── 15_spatial_transcriptomics/
├── 16_metabolomics/
├── 17_cross_omic_integration/    # unifies all omics layers into one ranked candidate list
├── organize_final_results.sh            # aggregates outputs across all stages
├── metadata/                            # de-identified sample manifests, pedigree files (not tracked)
├── references/                          # reference genome, known sites (not tracked)
└── logs/                                # SLURM job logs (not tracked)
```

## Infrastructure

All compute-intensive stages run on **Stanford's SCG compute cluster** via SLURM job scheduling.

## Tools & Libraries

- **Alignment & short-variant calling:** BWA-MEM2, GATK (HaplotypeCaller, CombineGVCFs, GenotypeGVCFs, CalculateGenotypePosteriors), DeepVariant, samtools
- **Repeat expansions:** ExpansionHunter
- **CNV:** GATK gCNV (CollectReadCounts, DenoiseReadCounts, ModelSegments, CallCopyRatioSegments)
- **Long-read SV:** minimap2, Sniffles2
- **Annotation:** ANNOVAR (RefGene, ClinVar, gnomAD, dbNSFP)
- **Proteomics:** Olink NPX analysis (pandas, SciPy, statsmodels)
- **Epigenomics:** methylKit (R/Bioconductor)
- **Chromatin accessibility:** MACS2
- **Single-cell multiome:** Scanpy, AnnData
- **Structural biology:** ESM-2 (variant effect prediction), AlphaFold Protein Structure Database
- **Spatial transcriptomics:** Scanpy, Squidpy (spatial autocorrelation)
- **Metabolomics:** pandas, SciPy, scikit-learn (PCA)
- **Statistics:** SciPy, statsmodels (Mann-Whitney U, Benjamini-Hochberg FDR correction)

## Status

Active - ongoing multi-omic analysis and cross-omic integration for candidate gene prioritization.
