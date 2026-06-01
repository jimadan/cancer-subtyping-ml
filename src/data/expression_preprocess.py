import pandas as pd
import numpy as np


class ExpressionPreprocessor:
    """
    Preprocessing pipeline for bulk RNA-seq gene expression data.

    Designed for cancer subtyping workflows:
    - sample alignment
    - missing value handling
    - normalization
    - gene filtering
    """

    def __init__(self,
                 log_transform=True,
                 min_variance=0.0):
        self.log_transform = log_transform
        self.min_variance = min_variance
    
    # -------------------------
    # Main pipeline
    # -------------------------
    def fit_transform(self, expr_raw: pd.DataFrame) -> pd.DataFrame:

        expr = self._format_matrix(expr_raw)
        expr = self._clean(expr)
        expr = self._impute(expr)

        if self.log_transform:
            expr = self._log_transform(expr)

        expr = self._filter_low_variance(expr)

        return expr

    # -------------------------
    # Step 1: format
    # -------------------------
    def _format_matrix(self, expr_raw: pd.DataFrame) -> pd.DataFrame:
        """
        Convert raw file into:
        rows = samples
        columns = genes
        """

        expr = expr_raw.set_index(expr_raw.columns[0])
        expr = expr.T

        expr = expr.apply(pd.to_numeric, errors="coerce")

        return expr

    # -------------------------
    # Step 2: clean
    # -------------------------
    def _clean(self, expr: pd.DataFrame) -> pd.DataFrame:

        expr = expr.dropna(axis=0, how="all")
        expr = expr.dropna(axis=1, how="all")

        return expr

    # -------------------------
    # Step 3: impute missing values
    # -------------------------
    def _impute(self, expr: pd.DataFrame) -> pd.DataFrame:

        # gene-wise mean imputation (biologically reasonable)
        return expr.fillna(expr.mean())

    # -------------------------
    # Step 4: log transform
    # -------------------------
    def _log_transform(self, expr: pd.DataFrame) -> pd.DataFrame:

        # RNA-seq standard transformation
        return np.log2(expr + 1)

    # -------------------------
    # Step 5: MAD filtering
    # -------------------------
    def _filter_low_variance(self, expr: pd.DataFrame) -> pd.DataFrame:

        #variances = expr.var(axis=0)
        #expr = expr.loc[:, variances > self.min_variance]
        
        mad = (expr - expr.median()).abs().median(axis=0)
        expr = expr.loc[:, mad > self.min_variance]

        return expr