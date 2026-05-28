import pandas as pd


def align_datasets(expr: pd.DataFrame,
                   clinical: pd.DataFrame):

    common = expr.index.intersection(clinical.index)

    expr = expr.loc[common].sort_index()
    clinical = clinical.loc[common].sort_index()

    assert expr.index.equals(clinical.index), \
        "Alignment error: sample mismatch"

    return expr, clinical