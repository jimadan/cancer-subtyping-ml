import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test

class SurvivalAnalysis:
    """
    Survival analysis for tumor subtypes.

    Inputs:
    - clinical data (OS time + event)
    - cluster labels

    Outputs:
    - Kaplan-Meier curves
    - log-rank test
    - Cox regression
    """

    def __init__(self, time_col="time", event_col="event"):
        self.time_col = time_col
        self.event_col = event_col

    def prepare_data(self, clinical: pd.DataFrame, clusters: pd.DataFrame):

        # 1. merge safely
        df = clinical.join(clusters, how="inner")

        # 2. drop fully missing rows
        df = df.dropna()

        # 3. -----------------------------
        # STANDARDIZE TIME COLUMN
        # -----------------------------
        possible_time_cols = [
            "time",
            "OS.time",
            "overall_survival_time",
            "days_to_death",
            "days_to_last_followup"
        ]

        if self.time_col not in df.columns:
            found = False
            for col in possible_time_cols:
                if col in df.columns:
                    df[self.time_col] = df[col]
                    found = True
                    break

            if not found:
                raise ValueError(
                    f"No survival time column found. "
                    f"Expected one of: {possible_time_cols}"
                )

        # 4. -----------------------------
        # STANDARDIZE EVENT COLUMN
        # -----------------------------
        possible_event_cols = [
            "event",
            "OS",
            "vital_status",
            "death_event"
        ]

        if self.event_col not in df.columns:
            found = False
            for col in possible_event_cols:
                if col in df.columns:
                    df[self.event_col] = df[col]
                    found = True
                    break

            if not found:
                raise ValueError(
                    f"No event column found. "
                    f"Expected one of: {possible_event_cols}"
                )

        # 5. convert event to binary (VERY IMPORTANT)
        df[self.event_col] = (
            df[self.event_col]
            .replace({
                "Dead": 1,
                "Alive": 0,
                "deceased": 1,
                "living": 0
            })
            .astype(int)
        )

        # 6. ensure numeric time
        df[self.time_col] = pd.to_numeric(df[self.time_col], errors="coerce")

        # 7. remove invalid rows AFTER conversion
        df = df.dropna(subset=[self.time_col, self.event_col])

        # 8. safety check
        if df[self.time_col].min() < 0:
            raise ValueError("Survival time contains negative values")

        # drop fully missing rows
        df = df.dropna()

        return df[[self.time_col, self.event_col, "cluster"]]

        # AttributeError: 'SurvivalAnalysis' object has no attribute 'kmf'"
    def plot_km(self, df: pd.DataFrame, ax=None):

        # create axis if not provided
        if ax is None:
            fig, ax = plt.subplots(figsize=(7, 5))

        clusters = df["cluster"].unique()

        kmf_dict = {}  # store fitted models (important for reproducibility)

        for c in sorted(clusters):

            group = df[df["cluster"] == c]

            kmf = KaplanMeierFitter()

            kmf.fit(
                durations=group[self.time_col],
                event_observed=group[self.event_col],
                label=f"Cluster {c}"
            )

            kmf.plot_survival_function(ax=ax, ci_show=True)

            kmf_dict[c] = kmf

        ax.set_title("Kaplan-Meier Survival by Cluster")
        ax.set_xlabel("Time")
        ax.set_ylabel("Survival probability")

        ax.legend()

        return ax, kmf_dict

    def logrank_test_between_clusters(self, df: pd.DataFrame):

        clusters = df["cluster"].unique()

        if len(clusters) != 2:
            return "Log-rank test implemented for 2 clusters only"

        c1 = df[df["cluster"] == clusters[0]]
        c2 = df[df["cluster"] == clusters[1]]

        result = logrank_test(
            c1[self.time_col], c2[self.time_col],
            event_observed_A=c1[self.event_col],
            event_observed_B=c2[self.event_col]
        )

        return {
            "test_statistic": result.test_statistic,
            "p_value": result.p_value
        }

    def cox_model(self, df: pd.DataFrame):

        df_encoded = df.copy()

        # Cluster as categorical variable
        df_encoded = pd.get_dummies(df_encoded, columns=["cluster"], drop_first=True)

        cph = CoxPHFitter()
        cph.fit(df_encoded, duration_col=self.time_col, event_col=self.event_col)

        return cph.summary