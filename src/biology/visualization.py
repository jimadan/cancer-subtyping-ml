import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path


def plot_marker_heatmap(expr_df, clusters, markers, output_path="results/figures/marker_heatmap.png"):
    """
    expr_df: samples × genes
    clusters: DataFrame with 'cluster'
    markers: list of genes
    """

    # subset = expr_df[markers].copy()

    # # attach cluster labels
    # subset["cluster"] = clusters["cluster"].values

    # # sort by cluster
    # subset = subset.sort_values("cluster")

    # data = subset.drop(columns=["cluster"])

    # plt.figure(figsize=(10, 8))

    # sns.heatmap(
    #     data.T,          # genes × samples (better for biology)
    #     cmap="vlag",
    #     center=0,
    #     cbar=True
    # )

    # plt.title("Cluster-specific marker gene expression")
    # plt.xlabel("Samples")
    # plt.ylabel("Marker genes")

    # plt.tight_layout()
    # plt.savefig(output_path, dpi=300)
    # plt.close()
    

    # 1. subset genes
    subset = expr_df[markers].copy()

    # 2. attach clusters
    subset["cluster"] = clusters["cluster"].values

    # 3. sort by cluster
    subset = subset.sort_values("cluster")

    # 4. split data
    data = subset.drop(columns=["cluster"])

    # 5. cluster colors
    palette = sns.color_palette("tab10", n_colors=subset["cluster"].nunique())
    color_map = {
        cluster: palette[i]
        for i, cluster in enumerate(sorted(subset["cluster"].unique()))
    }
    cluster_colors = subset["cluster"].map(color_map)

    # 6. plot heatmap WITH cluster bar
    g = sns.clustermap(
        data.T,
        col_colors=cluster_colors,
        col_cluster=False,
        row_cluster=False,
        cmap="vlag",
        center=0,
        figsize=(10, 8)
    )

    g.fig.suptitle("Cluster-specific marker gene expression", y=1.02)

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
