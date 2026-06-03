import pandas as pd
from sklearn.decomposition import PCA


class PCAReducer:

    def __init__(self, n_components=30):

        self.n_components = n_components
        self.model = PCA(n_components=n_components, random_state=42)

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
    
    def fit_pca_dynamic(self, X, variance=0.95):

        if (self.n_components != variance):
            self.n_components = variance
            self.model = PCA(n_components=variance, random_state=42)
        
        X_pca = self.model.fit_transform(X)

        return X_pca, self.model.n_components_