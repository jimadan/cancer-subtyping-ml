import pandas as pd
from lifelines import CoxPHFitter
from joblib import Parallel, delayed


def _fit_gene(gene, values, time, event):

    df = pd.DataFrame({
        "time": time,
        "event": event,
        "expr": values
    }).dropna()

    if len(df) < 10:
        return None

    try:
        cph = CoxPHFitter()
        cph.fit(
            df,
            duration_col="time",
            event_col="event"
        )

        p = cph.summary.loc["expr", "p"]

        return gene, p

    except Exception:
        return None
    
def cox_feature_selection(
    expr,
    clinical,
    top_n=2000,
    n_jobs=-1
):

    time = clinical["time"].values
    event = clinical["event"].values

    results = Parallel(
        n_jobs=n_jobs,
        backend="loky",
        verbose=10
    )(
        delayed(_fit_gene)(
            gene,
            expr[gene].values,
            time,
            event
        )
        for gene in expr.columns
    )

    results = [r for r in results if r is not None]

    ranked = pd.DataFrame(
        results,
        columns=["gene", "p"]
    )

    ranked = ranked.sort_values("p")

    selected_genes = ranked.head(top_n)["gene"]

    return expr[selected_genes], ranked