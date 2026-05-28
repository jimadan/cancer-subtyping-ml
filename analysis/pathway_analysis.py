import pandas as pd
import numpy as np


# ============================================================
# LOAD GSEA RESULTS
# ============================================================

def load_gsea_results(
    path="results/gsea/all_gsea_results.csv"
):

    df = pd.read_csv(path)

    if df is None or df.empty:

        raise ValueError(
            "GSEA results file is empty"
        )

    return df


# ============================================================
# FILTER SIGNIFICANT PATHWAYS
# ============================================================

def get_significant_pathways(
    df,
    fdr_threshold=0.25
):

    if df is None or df.empty:

        return pd.DataFrame()

    return df[
        df["FDR q-val"] < fdr_threshold
    ]


# ============================================================
# BUILD NES MATRIX
# ============================================================

def build_nes_matrix(
    df,
    fdr_threshold=None
):
    """
    Returns:
        pathways x clusters NES matrix
    """

    if df is None or df.empty:

        print("[WARNING] Empty dataframe")

        return pd.DataFrame()

    # optional filtering
    if fdr_threshold is not None:

        df = df[
            df["FDR q-val"] < fdr_threshold
        ]

    if df.empty:

        print("[WARNING] No pathways after filtering")

        return pd.DataFrame()

    matrix = df.pivot_table(

        index="Term",

        columns="cluster",

        values="NES",

        aggfunc="mean"

    )

    matrix = matrix.apply(pd.to_numeric, errors="coerce")
    matrix = matrix.fillna(0)

    return matrix


# ============================================================
# GET TOP PATHWAYS
# ============================================================

def get_top_pathways(
    df,
    cluster,
    n=10,
    ascending=False
):
    """
    Returns top enriched pathways
    """

    if df is None or df.empty:

        return pd.DataFrame()

    sub = df[
        df["cluster"] == cluster
    ]

    if sub.empty:

        return pd.DataFrame()

    top = (
        sub
        .sort_values(
            "NES",
            ascending=ascending
        )
        .head(n)
    )

    return top


# ============================================================
# GET TOP AGING PATHWAYS
# ============================================================

AGING_KEYWORDS = [

    "AGING",
    "SENESC",
    "P53",
    "DNA REPAIR",
    "INFLAM",
    "TNF",
    "INTERFERON",
    "ROS",
    "OXID",
    "HYPOX",
    "AUTOPH",
    "MITO",
    "MTOR",
    "IMMUNE",
    "CYTOKINE",
    "IL6",
    "JAK",
    "SASP",
    "FOXO",
    "TELO",
    "EMT"
]


def get_aging_pathways(
    df,
    cluster=None,
    fdr_threshold=0.35
):

    if df is None or df.empty:

        return pd.DataFrame()

    sub = df.copy()

    if cluster is not None:

        sub = sub[
            sub["cluster"] == cluster
        ]

    sub = sub[
        sub["FDR q-val"] < fdr_threshold
    ]

    sub["term_upper"] = (
        sub["Term"]
        .astype(str)
        .str.upper()
    )

    mask = sub["term_upper"].apply(

        lambda x: any(
            k in x
            for k in AGING_KEYWORDS
        )
    )

    return sub[mask]


# ============================================================
# PATHWAY GROUPS
# ============================================================

PATHWAY_GROUPS = {

    "senescence": [
        "P53",
        "APOPT",
        "SENESC",
        "TELO",
        "DNA REPAIR"
    ],

    "inflammation": [
        "TNF",
        "INTERFERON",
        "INFLAM",
        "IMMUNE",
        "CYTOKINE",
        "JAK",
        "IL6",
        "NF"
    ],

    "cell_cycle": [
        "CELL CYCLE",
        "E2F",
        "G2M",
        "MITOTIC",
        "CHECKPOINT",
        "MYC"
    ],

    "metabolism": [
        "MTOR",
        "ROS",
        "OXID",
        "MITO",
        "HYPOX",
        "AUTOPH"
    ],

    "stemness": [
        "WNT",
        "EMT",
        "NOTCH",
        "HEDGEHOG",
        "TGF"
    ]
}


# ============================================================
# ASSIGN PATHWAY GROUP
# ============================================================

def assign_pathway_group(term):

    term_upper = str(term).upper()

    for group, keywords in PATHWAY_GROUPS.items():

        if any(
            k in term_upper
            for k in keywords
        ):

            return group

    return "other"


# ============================================================
# ADD PATHWAY GROUPS
# ============================================================

def annotate_pathway_groups(df):

    if df is None or df.empty:

        return pd.DataFrame()

    df = df.copy()

    df["pathway_group"] = (
        df["Term"]
        .apply(assign_pathway_group)
    )

    return df


# ============================================================
# BUILD GROUP-LEVEL MATRIX
# ============================================================

def build_group_matrix(
    df,
    fdr_threshold=0.25
):
    """
    Returns:
        pathway_group x cluster NES matrix
    """

    if df is None or df.empty:

        return pd.DataFrame()

    df = df[
        df["FDR q-val"] < fdr_threshold
    ]

    df = annotate_pathway_groups(df)

    summary = (

        df

        .groupby(
            ["pathway_group", "cluster"]
        )["NES"]

        .mean()

        .reset_index()
    )

    matrix = summary.pivot_table(

        index="pathway_group",

        columns="cluster",

        values="NES",

        aggfunc="mean"

    )

    matrix = matrix.apply(pd.to_numeric, errors="coerce")
    matrix = matrix.fillna(0)

    return matrix


# ============================================================
# CLUSTER SUMMARY
# ============================================================

def summarize_cluster(df, cluster):

    sub = df[
        df["cluster"] == cluster
    ]

    if sub.empty:

        return {}

    summary = {

        "n_pathways": len(sub),

        "mean_NES": sub["NES"].mean(),

        "top_positive": (
            sub
            .sort_values(
                "NES",
                ascending=False
            )
            .iloc[0]["Term"]
        ),

        "top_negative": (
            sub
            .sort_values(
                "NES",
                ascending=True
            )
            .iloc[0]["Term"]
        )
    }

    return summary