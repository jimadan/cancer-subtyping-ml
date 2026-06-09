import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from pathlib import Path


def plot_aging_scores(aging_df, output="results/figures/aging_scores.png"):

    if aging_df is None or aging_df.empty:
        print("[WARNING] Skipping aging score plot: no aging scores")
        return None

    df = aging_df.copy()

    # ----------------------------
    # split score vs n columns
    # ----------------------------
    score_cols = [
        c
        for c in df.columns
        if c != "cluster"
        and not c.endswith("_n")
        and pd.api.types.is_numeric_dtype(df[c])
    ]

    if not score_cols:
        print("[WARNING] Skipping aging score plot: no numeric aging scores")
        return None

    score_df = df[["cluster", *score_cols]]
    n_df = df[[c for c in df.columns if c.endswith("_n")]]

    # rename n columns to match processes
    n_df.columns = [c.replace("_n", "") for c in n_df.columns]

    n_df = n_df.copy()
    n_df["cluster"] = df["cluster"].to_numpy()
    n_df = n_df.set_index("cluster")

    # ----------------------------
    # reshape scores for barplot
    # ----------------------------
    melted = score_df.melt(
        id_vars="cluster",
        var_name="aging_process",
        value_name="score"
    )

    # ----------------------------
    # figure layout
    # ----------------------------
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # ----------------------------
    # LEFT: aging scores
    # ----------------------------
    sns.barplot(
        data=melted,
        x="aging_process",
        y="score",
        hue="cluster",
        ax=axes[0]
    )

    axes[0].set_title("Aging-related pathway activity")
    axes[0].set_xlabel("Aging process")
    axes[0].set_ylabel("Score")
    axes[0].tick_params(axis='x', rotation=30)

    # ----------------------------
    # RIGHT: pathway support (n)
    # ----------------------------
    sns.heatmap(
        n_df,
        cmap="Greys",
        annot=True,
        fmt="d",
        ax=axes[1]
    )

    axes[1].set_title("Pathway support (n)")
    axes[1].set_xlabel("Aging process")
    axes[1].set_ylabel("Cluster")

    # ----------------------------
    # save
    # ----------------------------
    plt.tight_layout()

    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output, dpi=300)
    plt.close()

    return output

def plot_aging_heatmaps(aging_df, output="results/figures/aging_heatmaps.png"):

    if aging_df is None or aging_df.empty:
        print("[WARNING] Skipping aging heatmap: no aging scores")
        return None

    # -----------------------------
    # Split score vs counts
    # -----------------------------
    score_df = aging_df.set_index("cluster")
    count_df = aging_df.set_index("cluster")

    score_cols = [
        c
        for c in score_df.columns
        if not c.endswith("_n")
        and pd.api.types.is_numeric_dtype(score_df[c])
    ]
    count_cols = [c for c in score_df.columns if c.endswith("_n")]

    if not score_cols or not count_cols:
        print("[WARNING] Skipping aging heatmap: missing score or count columns")
        return None

    score_df = score_df[score_cols]
    count_df = count_df[count_cols]

    # Rename count columns (remove "_n")
    count_df.columns = [c.replace("_n", "") for c in count_df.columns]

    # -----------------------------
    # Plot
    # -----------------------------
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    sns.heatmap(
        score_df,
        cmap="coolwarm",
        center=0,
        annot=True,
        ax=axes[0]
    )
    axes[0].set_title("Aging program activity")

    sns.heatmap(
        count_df,
        cmap="Greys",
        annot=True,
        ax=axes[1]
    )
    axes[1].set_title("Number of contributing pathways")

    plt.tight_layout()

    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output, dpi=300)
    plt.close()

    return output
