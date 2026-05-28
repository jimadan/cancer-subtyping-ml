import pandas as pd

def load_expression(path: str) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")

def load_clinical(path: str) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")