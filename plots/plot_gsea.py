import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path


def plot_pathway_heatmap(df):

    if df is None or df.empty:
        print("[WARNING] Skipping pathway heatmap: no GSEA results")
        return None

    df = df.copy()

    # ensure correct types
    df["NES"] = pd.to_numeric(df["NES"], errors="coerce")

    df = df.dropna(subset=["NES"])

    pivot = df.pivot_table(
        index="Term",
        columns="cluster",
        values="NES",
        aggfunc="mean"
    )

    if pivot.empty:
        print("[WARNING] Skipping pathway heatmap: empty NES matrix")
        return None

    # safety
    pivot = pivot.astype(float)

    plt.figure(figsize=(10, 12))

    sns.heatmap(
        pivot,
        cmap="coolwarm",
        center=0
    )

    plt.title("Pathway activity (NES)")
    plt.tight_layout()
    output = Path("results/figures/pathway_heatmap.png")
    output.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output, dpi=300)
    plt.close()
    return output
