# Computational Oncology Pipeline for Tumor Subtype Discovery

A modular transcriptomic analysis framework for unsupervised tumor subtype discovery, functional enrichment analysis, and survival modeling using TCGA RNA-seq data.

---

## Overview

Cancer is a highly heterogeneous disease in which patients with similar histopathological diagnoses may exhibit distinct molecular profiles and clinical outcomes.

This project implements an end-to-end computational oncology pipeline that integrates transcriptomic analysis with unsupervised machine learning to identify biologically and clinically meaningful tumor subtypes.

The framework is designed as a reproducible and modular research pipeline for exploratory analysis in computational oncology and bioinformatics.

---

## Objectives

- Discover latent tumor subtypes from RNA-seq expression data
- Evaluate whether transcriptomic structure correlates with patient survival
- Identify pathway-level differences between subtypes
- Compute aging-related transcriptional signatures
- Provide a reproducible pipeline for multi-stage bioinformatics analysis

---

## Dataset

### TCGA Breast Cancer (BRCA)

The pipeline uses publicly available data from:

- TCGA RNA-seq gene expression profiles (BRCA cohort)
- Matched clinical survival metadata

---

## Pipeline Overview

```text
RNA-seq expression + clinical data
        ->
Preprocessing (normalization + filtering)
        ->
Feature engineering (survival gene selection)
        ->
PCA dimensionality reduction
        ->
Clustering (model selection + optimal k)
        ->
UMAP visualization embedding
        ->
Subtype assignment
        ->
GSEA functional enrichment analysis
        ->
Pathway aggregation (NES matrix)
        ->
Aging signature scoring
        ->
Survival analysis (Kaplan-Meier)
```

This project was developed iteratively using AI-assisted coding tools for refactoring and documentation, with full human supervision and validation.
