from src.biology.differential_expression import run_biological_analysis
from src.biology.visualization import plot_marker_heatmap


def run_biology_pipeline(expr_df, clusters):

    print("\n[STEP] Biological interpretation...")

    results, markers = run_biological_analysis(
        expr_df,
        clusters,
        n_clusters=3
    )

    plot_marker_heatmap(expr_df, clusters, markers[0])

    return results, markers