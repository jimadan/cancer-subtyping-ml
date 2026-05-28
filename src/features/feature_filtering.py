import pandas as pd
import numpy as np

class FeatureFilter:
    """
    Feature selection for gene expression data.

    Goal:
    Reduce dimensionality while keeping biologically informative genes
    for cancer subtyping (clustering).
    """

    def __init__(self,
                 min_variance=0.5,
                 top_n_genes=None,
                 method="variance"):
        """
        Parameters
        ----------
        min_variance : float
            Minimum variance threshold for gene filtering.

        top_n_genes : int or None
            If set, keeps only top N most variable genes.

        method : str
            'variance' or 'mad'
        """

        self.min_variance = min_variance
        self.top_n_genes = top_n_genes
        self.method = method

    # -------------------------
    # Main function
    # -------------------------
    def fit_transform(self, expr: pd.DataFrame) -> pd.DataFrame:
        """
        Apply feature filtering to gene expression matrix.
        """

        if self.method == "variance":
            expr = self._variance_filter(expr)

        elif self.method == "mad":
            expr = self._mad_filter(expr)

        else:
            raise ValueError("Method must be 'variance' or 'mad'")

        if self.top_n_genes is not None:
            expr = self._top_n(expr)

        return expr

    def _variance_filter(self, expr: pd.DataFrame) -> pd.DataFrame:
        """
        Remove low-variance genes.
        """

        variances = expr.var(axis=0)

        filtered = expr.loc[:, variances > self.min_variance]

        print(f"[INFO] Variance filtering: {expr.shape[1]} → {filtered.shape[1]} genes")

        return filtered

    def _mad_filter(self, expr: pd.DataFrame) -> pd.DataFrame:
        """
        Median Absolute Deviation filtering.
        More robust to outliers than variance.
        """

        mad = expr.mad(axis=0)

        filtered = expr.loc[:, mad > self.min_variance]

        print(f"[INFO] MAD filtering: {expr.shape[1]} → {filtered.shape[1]} genes")

        return filtered

    def _top_n(self, expr: pd.DataFrame) -> pd.DataFrame:
        """
        Keep top N most variable genes.
        """

        variances = expr.var(axis=0)

        top_genes = variances.sort_values(ascending=False).head(self.top_n_genes).index

        filtered = expr.loc[:, top_genes]

        print(f"[INFO] Top-{self.top_n_genes} genes selected")

        return filtered