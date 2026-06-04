import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path


def plot_pathway_heatmap(df):

    if df is None or df.empty:
        print("[WARNING] Skipping pathway heatmap: no GSEA results")
        return None

    df = df.copy()

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

    pivot = pivot.astype(float)

    import matplotlib.pyplot as plt
    import seaborn as sns
    from pathlib import Path

    fig, ax = plt.subplots(
        figsize=(16, 14),
        constrained_layout=True
    )

    sns.heatmap(
        pivot,
        cmap="coolwarm",
        center=0,
        ax=ax,
        cbar_kws={"label": "NES"}
    )

    ax.set_title("Pathway activity (NES)", pad=20)
    ax.tick_params(axis='y', labelsize=8)
    ax.tick_params(axis='x', labelrotation=0)

    output = Path("results/figures/pathway_heatmap.png")
    output.parent.mkdir(parents=True, exist_ok=True)

    plt.savefig(output, dpi=300, bbox_inches="tight")
    plt.close()

    return output

def plot_pathway_support_heatmap(df):

    if df is None or df.empty:
        print("[WARNING] Skipping pathway support heatmap: no GSEA results")
        return None

    df = df.copy()

    n_matrix = df.groupby(["Term", "cluster"]).size().unstack(fill_value=0)

    if n_matrix.empty:
        print("[WARNING] Skipping pathway support heatmap: empty matrix")
        return None

    n_matrix = n_matrix.sort_index()

    import matplotlib.pyplot as plt
    import seaborn as sns
    from pathlib import Path

    fig, ax = plt.subplots(
        figsize=(16, 14),
        constrained_layout=True
    )

    sns.heatmap(
        n_matrix,
        cmap="Greys",
        annot=True,
        fmt="d",
        ax=ax,
        cbar_kws={"label": "Pathway count"}
    )

    ax.set_title("Pathway support (number of contributing pathways)", pad=20)
    ax.tick_params(axis='y', labelsize=8)
    ax.set_xlabel("Cluster")
    ax.set_ylabel("Pathway (Term)")

    output = Path("results/figures/pathway_support_heatmap.png")
    output.parent.mkdir(parents=True, exist_ok=True)

    plt.savefig(output, dpi=300, bbox_inches="tight")
    plt.close()

    return output