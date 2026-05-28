import numpy as np
import pandas as pd
from scipy.stats import ttest_ind


def differential_expression(expr_df: pd.DataFrame, clusters: pd.DataFrame, cluster_id: int):
    """
    expr_df: samples × genes
    clusters: DataFrame with column 'cluster'
    """

    group1 = expr_df.loc[clusters["cluster"] == cluster_id]
    group2 = expr_df.loc[clusters["cluster"] != cluster_id]

    results = []

    for gene in expr_df.columns:

        # statistical test
        stat, pval = ttest_ind(
            group1[gene],
            group2[gene],
            equal_var=False
        )

        # log fold change
        logFC = np.log2(
            group1[gene].mean() + 1e-9
        ) - np.log2(
            group2[gene].mean() + 1e-9
        )

        results.append([gene, logFC, pval])

    res_df = pd.DataFrame(results, columns=["gene", "logFC", "pvalue"])

    # ranking score (optional)
    res_df["score"] = -np.log10(res_df["pvalue"] + 1e-12)

    return res_df.sort_values("pvalue")


def get_top_markers(de_results: pd.DataFrame, top_n: int = 100):
    """
    Select top marker genes for pathway analysis
    """

    return (
        de_results
        .sort_values("pvalue")
        .head(top_n)["gene"]
        .tolist()
    )


def run_biological_analysis(expr_df, clusters, n_clusters=3):
    """
    Full wrapper: DE + marker extraction per cluster
    """

    all_results = {}
    all_markers = {}

    unique_clusters = np.sort(clusters["cluster"].unique())

    for k in unique_clusters:

        print(f"[BIO] Running DE for cluster {k}")

        de = differential_expression(expr_df, clusters, k)

        # skip empty / tiny clusters safely
        if de is None or len(de) == 0:
            print(f"[BIO] Skipping cluster {k} (no results)")
            continue

        markers = get_top_markers(de, top_n=100)

        all_results[k] = de
        all_markers[k] = markers

        de.to_csv(f"results/de_cluster_{k}.csv", index=False)


    return all_results, all_markers