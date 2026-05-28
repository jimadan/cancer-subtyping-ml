import pandas as pd
from sklearn.decomposition import PCA


class PCAReducer:

    def __init__(self, n_components=30):

        self.n_components = n_components
        self.model = PCA(n_components=n_components)

    def fit_transform(self, X):

        X_pca = self.model.fit_transform(X)

        columns = [
            f"PC{i+1}"
            for i in range(self.n_components)
        ]

        X_pca = pd.DataFrame(
            X_pca,
            index=X.index,
            columns=columns
        )

        return X_pca