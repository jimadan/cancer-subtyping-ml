import matplotlib.pyplot as plt
import seaborn as sns


def plot_aging_scores(aging_df, output="results/figures/aging_scores.png"):

    melted = aging_df.melt(
        id_vars="cluster",
        var_name="aging_process",
        value_name="score"
    )

    plt.figure(figsize=(8,5))

    sns.barplot(
        data=melted,
        x="aging_process",
        y="score",
        hue="cluster"
    )

    plt.title("Aging-related pathway activity by cluster")
    plt.xticks(rotation=30)
    plt.tight_layout()

    plt.savefig(output, dpi=300)
    #plt.close()

def plot_aging_heatmap(aging_df, output="results/figures/aging_heatmap.png"):

    matrix = aging_df.set_index("cluster")

    plt.figure(figsize=(6,4))

    sns.heatmap(
        matrix,
        cmap="coolwarm",
        center=0,
        annot=True
    )

    plt.title("Aging program activity")

    plt.tight_layout()
    plt.savefig(output, dpi=300)
    #plt.close()