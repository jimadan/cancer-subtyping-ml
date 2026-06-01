import pandas as pd
from lifelines import CoxPHFitter


class SurvivalFeatureSelector:

    def __init__(self, p_threshold=0.05):

        self.p_threshold = p_threshold

    def fit_transform(self, expr, clinical):

        selected_genes = []

        for gene in expr.columns:

            df = pd.DataFrame({
                "gene": expr[gene],
                "time": clinical["time"],
                "event": clinical["event"]
            })

            try:
                cph = CoxPHFitter()
                cph.fit(df,
                        duration_col="time",
                        event_col="event")

                p = cph.summary.loc["gene", "p"]

                if p < self.p_threshold:
                    selected_genes.append(gene)

            except Exception:
                continue

        print(f"[INFO] Selected genes: {len(selected_genes)}")

        return expr[selected_genes]
