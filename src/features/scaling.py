from sklearn.preprocessing import StandardScaler
import pandas as pd


class ExpressionScaler:
    """
    Standardize gene expression features.

    Centers each gene to mean=0 and variance=1.
    """

    def __init__(self):
        self.scaler = StandardScaler()

    def fit_transform(self, X):

        X_scaled = self.scaler.fit_transform(X)

        X_scaled = pd.DataFrame(
            X_scaled,
            index=X.index,
            columns=X.columns
        )

        return X_scaled