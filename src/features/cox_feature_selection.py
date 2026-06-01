import pandas as pd
from lifelines import CoxPHFitter


def cox_feature_selection(expr: pd.DataFrame, clinical: pd.DataFrame, top_n=2000): #500

    results = []

    for gene in expr.columns:

        df = pd.DataFrame({
            "time": clinical["time"],
            "event": clinical["event"],
            "expr": expr[gene]
        }).dropna()

        try:
            cph = CoxPHFitter()
            cph.fit(df, duration_col="time", event_col="event")

            p = cph.summary.loc["expr", "p"]

            results.append((gene, p))

        except Exception:
            continue

    ranked = pd.DataFrame(results, columns=["gene", "p"])

    if ranked.empty:
        raise ValueError("Cox feature selection did not fit any gene models")

    ranked = ranked.sort_values("p")

    selected_genes = ranked.head(top_n)["gene"].tolist()

    return expr[selected_genes], ranked
