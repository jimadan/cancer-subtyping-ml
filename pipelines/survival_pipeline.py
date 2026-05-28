import matplotlib.pyplot as plt
from src.evaluation.survival import SurvivalAnalysis
from lifelines.statistics import multivariate_logrank_test


def run_survival_pipeline(clinical, clusters, X_umap, labels):

    print("\n[STEP] Survival analysis...")

    sa = SurvivalAnalysis()
    df = sa.prepare_data(clinical, clusters)

    plt.figure()
    sa.plot_km(df)
    plt.savefig("results/figures/km_plot.png")

    stats = sa.logrank_test_between_clusters(df)
    print("[INFO] Log-rank test:", stats)

    cox_results = sa.cox_model(df)
    print(cox_results)

    results = multivariate_logrank_test(
        event_durations=df["time"],
        groups=df["cluster"],
        event_observed=df["event"]
    )

    print(results.summary)

    # UMAP plot
    # plt.figure(figsize=(8,6))
    # plt.scatter(X_umap[:,0], X_umap[:,1], c=labels, cmap="tab10", s=15)
    # plt.savefig("results/figures/umap_clusters.png")