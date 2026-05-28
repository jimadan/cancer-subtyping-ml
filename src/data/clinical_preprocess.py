import pandas as pd


class ClinicalPreprocessor:
    """
    Clinical metadata preprocessing for survival analysis.
    """

    def __init__(self,
                 required_cols=None):

        if required_cols is None:
            required_cols = ["time", "event"]

        self.required_cols = required_cols

    # ---------------------------------
    def transform(self, clinical_raw: pd.DataFrame) -> pd.DataFrame:

        clinical = self._format(clinical_raw)
        clinical = self._select_columns(clinical)
        clinical = self._drop_missing(clinical)

        return clinical

    # ---------------------------------
    def _format(self, clinical: pd.DataFrame) -> pd.DataFrame:

        clinical = clinical.rename(columns={
            "sample": "sample_id",
            "OS.time": "time",
            "OS": "event"
        })

        clinical = clinical.set_index("sample_id")

        return clinical

    # ---------------------------------
    def _select_columns(self, clinical: pd.DataFrame) -> pd.DataFrame:

        return clinical[self.required_cols]

    # ---------------------------------
    def _drop_missing(self, clinical: pd.DataFrame) -> pd.DataFrame:

        return clinical.dropna()