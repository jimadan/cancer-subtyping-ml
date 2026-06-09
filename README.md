# Computational Oncology Pipeline for Tumor Subtype Discovery

A modular transcriptomic analysis framework for unsupervised tumor subtype discovery,
functional enrichment analysis, aging-program scoring, and survival modeling using
TCGA RNA-seq data.

## Overview

Cancer is a heterogeneous disease in which patients with similar histopathological
diagnoses may show distinct molecular profiles and clinical outcomes. This project
implements an end-to-end computational oncology pipeline that integrates RNA-seq
preprocessing, survival-informed feature selection, unsupervised clustering, pathway
enrichment, biological process annotation, aging signature analysis, and
Kaplan-Meier survival evaluation.

The codebase is organized as reusable pipeline stages so each analysis step can be
run, inspected, or extended independently.

## Objectives

- Discover latent tumor subtypes from RNA-seq expression data.
- Select survival-associated genes before clustering.
- Evaluate whether transcriptomic subtypes correlate with patient survival.
- Identify pathway-level differences between subtypes using GSEA.
- Annotate enriched pathways into broad biological process groups.
- Compute aging-related transcriptional program scores per cluster.
- Generate reproducible figures and CSV outputs for downstream interpretation.

## Dataset

The pipeline expects TCGA Breast Cancer (BRCA) inputs under `data/raw/`:

- `HiSeqV2.txt`: RNA-seq gene expression matrix.
  - Dataset page: [TCGA BRCA HiSeqV2 on UCSC Xena](https://xenabrowser.net/datapages/?dataset=TCGA.BRCA.sampleMap%2FHiSeqV2&host=https%3A%2F%2Ftcga.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443)
  - Direct download: [TCGA.BRCA.sampleMap/HiSeqV2.gz](https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.BRCA.sampleMap%2FHiSeqV2.gz)
- `BRCA_survival.txt`: matched clinical survival metadata.
  - Dataset page: [TCGA BRCA survival metadata on UCSC Xena](https://xenabrowser.net/datapages/?dataset=survival%2FBRCA_survival.txt&host=https%3A%2F%2Ftcga.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443)
  - Direct download: [survival/BRCA_survival.txt](https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/survival%2FBRCA_survival.txt)

These paths are configured in `config.py`.

## Pipeline

```text
RNA-seq expression + clinical survival data
        ->
Preprocessing and sample alignment
        ->
Variance-based gene filtering
        ->
Univariate Cox feature selection
        ->
Expression scaling
        ->
PCA dimensionality reduction
        ->
GMM clustering in PCA space with k=2..5 silhouette selection
        ->
UMAP visualization embedding
        ->
GSEA against Hallmark, KEGG, Reactome, and GO Biological Process
        ->
Biological process annotation and NES matrices
        ->
Aging program scoring
        ->
Kaplan-Meier survival analysis
```

Run the complete workflow with:

```bash
python scripts/run_full_pipeline.py
```

## Biological Process Annotation

GSEA results are annotated with shared pathway groups from
`analysis/pathway_analysis.py`. The current biological process categories are:

- `senescence`: P53, apoptosis, senescence, telomere, and DNA repair pathways.
- `inflammation`: TNF, interferon, immune, cytokine, JAK/IL6/NF-related pathways.
- `cell_cycle`: E2F, G2M, mitotic, checkpoint, MYC, and cell-cycle pathways.
- `metabolism`: MTOR, ROS, oxidative, mitochondrial, hypoxia, and autophagy pathways.
- `stemness`: WNT, EMT, NOTCH, Hedgehog, and TGF-related pathways.
- `other`: enriched pathways that do not match the broad keyword groups.

The GSEA stage currently searches:

- `MSigDB_Hallmark_2020`
- `KEGG_2021_Human`
- `Reactome_2022`
- `GO_Biological_Process_2023`

## Aging Analysis

Aging scores are computed from enriched pathways by matching relaxed biological
program keywords and calculating weighted NES scores. The current aging programs are:

- `inflammaging`
- `mitochondrial_stress`
- `proteostasis_loss`
- `cell_cycle_dysregulation`
- `identity_drift`
- `dna_damage_response`

For each cluster, the pipeline stores both the score and the number of contributing
pathways for each aging program.

## Main Outputs

The full pipeline writes results under `results/`:

- `results/gsea/all_gsea_results.csv`: combined significant GSEA results with
  `pathway_group` annotations.
- `results/aging/aging_scores.csv`: cluster-level aging program scores and pathway
  support counts.
- `results/figures/umap.png`: UMAP visualization of selected tumor subtypes.
- `results/figures/pathway_heatmap.png`: pathway NES heatmap by cluster.
- `results/figures/pathway_support_heatmap.png`: pathway support counts by cluster.
- `results/figures/aging_scores.png`: bar plot of aging program scores with support.
- `results/figures/aging_heatmaps.png`: aging score and support heatmaps.

Additional survival outputs are produced by `pipelines/survival_pipeline.py`.

## Project Structure

```text
analysis/      Pathway and aging analysis helpers
pipelines/     End-to-end pipeline stages
plots/         GSEA and aging visualization functions
scripts/       Full pipeline entry point
src/           Data, feature, clustering, biology, and evaluation modules
notebooks/     Exploratory and methods notebooks
config.py      Central project paths
```

## Notes

This project was developed iteratively using AI-assisted coding tools for refactoring
and documentation, with full human supervision and validation.
