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
                "OS.time": clinical["OS.time"],
                "OS": clinical["OS"]
            })

            try:
                cph = CoxPHFitter()
                cph.fit(df,
                        duration_col="OS.time",
                        event_col="OS")

                p = cph.summary.loc["gene", "p"]

                if p < self.p_threshold:
                    selected_genes.append(gene)

            except:
                continue

        print(f"[INFO] Selected genes: {len(selected_genes)}")

        return expr[selected_genes]