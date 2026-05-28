import pandas as pd


# ============================================================
# BROAD AGING PROGRAMS (NOT STRICT KEYWORDS)
# ============================================================

AGING_PROGRAMS = {

    "inflammaging": [
        "INFLAMMATORY RESPONSE",
        "TNFA SIGNALING",
        "INTERFERON",
        "IL6",
        "JAK",
        "NF",
        "CYTOKINE"
    ],

    "mitochondrial_stress": [
        "OXIDATIVE PHOSPHORYLATION",
        "REACTIVE OXYGEN",
        "MITOCHONDRIA"
    ],

    "proteostasis_loss": [
        "UNFOLDED PROTEIN",
        "PROTEASOME",
        "MTORC1",
        "AUTOPHAGY"
    ],

    "cell_cycle_dysregulation": [
        "G2M CHECKPOINT",
        "E2F TARGETS",
        "CELL CYCLE"
    ],

    "identity_drift": [
        "EPITHELIAL MESENCHYMAL",
        "WNT",
        "TGF",
        "NOTCH"
    ],

    "dna_damage_response": [
        "DNA REPAIR",
        "P53",
        "APOPTOSIS"
    ]
}


# ============================================================
# LOAD RESULTS
# ============================================================

def load_results(path="results/gsea/all_gsea_results.csv"):
    return pd.read_csv(path)


# ============================================================
# SOFT AGING SCORING (NO OVER-FILTERING)
# ============================================================

def compute_aging_scores(df):

    scores = []

    for cl in df["cluster"].unique():

        sub = df[df["cluster"] == cl].copy()

        row = {"cluster": cl}

        sub["term_upper"] = sub["Term"].str.upper()

        for program, keywords in AGING_PROGRAMS.items():

            matched = sub[
                sub["term_upper"].apply(
                    lambda x: any(k in x for k in keywords)
                )
            ]

            # IMPORTANT: very relaxed
            matched = matched.loc[sub["FDR q-val"] < 0.50].copy()

            row[program] = matched["NES"].mean() if len(matched) > 0 else 0

        scores.append(row)

    return pd.DataFrame(scores)


# ============================================================
# SIMPLE INTERPRETATION
# ============================================================

def interpret_aging(row):

    if row["inflammaging"] > 0 and row["cell_cycle_dysregulation"] > 0:
        return "Inflammatory proliferative aging"

    if row["mitochondrial_stress"] > 0:
        return "Mitochondrial stress aging"

    if row["proteostasis_loss"] > 0:
        return "Proteostasis aging"

    if row["identity_drift"] > 0:
        return "Identity drift aging"

    return "Mixed / early aging state"


# ============================================================
# PIPELINE
# ============================================================

def run_aging_analysis():

    df = load_results()

    scores = compute_aging_scores(df)

    scores["state"] = scores.apply(interpret_aging, axis=1)

    scores.to_csv("results/gsea/aging_analysis.csv", index=False)

    return scores