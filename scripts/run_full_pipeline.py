from pathlib import Path
import sys

import matplotlib.pyplot as plt
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

# Pipelines
from pipelines.preprocessing_pipeline import run_preprocessing
from pipelines.feature_pipeline import run_feature_pipeline
from pipelines.clustering_pipeline import run_clustering_pipeline
from pipelines.gsea_pipeline import run_gsea_pipeline
from pipelines.survival_pipeline import run_survival_pipeline

# Analysis
from analysis.pathway_analysis import build_nes_matrix
from analysis.aging_analysis import compute_aging_scores

# Plots
from plots.plot_gsea import plot_pathway_heatmap, plot_pathway_support_heatmap
from plots.plot_aging import plot_aging_scores, plot_aging_heatmaps


def main():

    Path("results/figures").mkdir(parents=True, exist_ok=True)
    Path("results/gsea").mkdir(parents=True, exist_ok=True)
    Path("results/aging").mkdir(parents=True, exist_ok=True)

    bundle = {}

    # -------------------------
    # 1. PREPROCESSING
    # -------------------------
    expr, clinical = run_preprocessing()

    expr_filtered, expr_survival, X_pca = run_feature_pipeline(expr, clinical)

    # -------------------------
    # 2. CLUSTERING
    # -------------------------
    X_umap, labels, best_k, best_model = run_clustering_pipeline(X_pca)

    print(f"[INFO] Selected k = {best_k}")

    clusters = pd.DataFrame({
        "sample": expr_filtered.index,
        "cluster": labels
    }).set_index("sample")

    bundle["clusters"] = clusters
    bundle["clinical"] = clinical

    # -------------------------
    # 3. GSEA
    # -------------------------
    gsea_df = run_gsea_pipeline(expr_survival, clusters)

    bundle["gsea_raw"] = gsea_df

    # -------------------------
    # 4. PATHWAY ANALYSIS
    # -------------------------
    bundle["nes_matrix"] = build_nes_matrix(gsea_df)

    # -------------------------
    # 5. AGING ANALYSIS
    # -------------------------
    aging_scores = compute_aging_scores(gsea_df)
    bundle["aging_scores"] = aging_scores

    # -------------------------
    # 6. PLOTTING
    # -------------------------

    # GSEA heatmap
    plot_pathway_heatmap(gsea_df)
    plot_pathway_support_heatmap(gsea_df)

    # Aging plots
    plot_aging_scores(aging_scores)
    plot_aging_heatmaps(aging_scores)

    # UMAP plot
    plt.figure(figsize=(8,6))
    plt.scatter(X_umap[:,0], X_umap[:,1], c=labels, cmap="tab10", s=15)
    plt.title(f"UMAP tumor subtypes (k={best_k})")
    plt.savefig("results/figures/umap.png")
    plt.close()

    # -------------------------
    # 7. SURVIVAL
    # -------------------------
    run_survival_pipeline(clinical, clusters, X_umap, labels)

    # -------------------------
    # 8. SAVE RESULTS
    # -------------------------
    gsea_df.to_csv("results/gsea/all_gsea_results.csv", index=False)

    aging_scores.to_csv("results/aging/aging_scores.csv", index=False)

    print("[INFO] Pipeline completed successfully")

    return bundle


if __name__ == "__main__":
    main()
