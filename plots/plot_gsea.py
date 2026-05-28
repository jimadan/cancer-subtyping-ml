import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def plot_pathway_heatmap(df):

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
    plt.savefig("results/figures/pathway_heatmap.png", dpi=300)
    #plt.close()