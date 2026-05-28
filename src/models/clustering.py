import pandas as pd
import numpy as np

from sklearn.cluster import (
    KMeans,
    AgglomerativeClustering
)

from sklearn.mixture import GaussianMixture

from sklearn.metrics import (
    silhouette_score,
    davies_bouldin_score
)


class ClusteringModel:
    """
    Unsupervised clustering for tumor subtyping using gene expression data.

    Supports:
    - KMeans
    - Hierarchical clustering
    - Gaussian Mixture Models (GMM)

    Includes:
    - silhouette evaluation
    - Davies-Bouldin score
    - cluster assignment output
    """

    def __init__(
        self,
        method="kmeans",
        n_clusters=3,
        random_state=42
    ):

        self.method = method
        self.n_clusters = n_clusters
        self.random_state = random_state

        self.model = None
        self.labels_ = None

    def fit(self, X: pd.DataFrame):

        # -------------------------
        # KMeans
        # -------------------------
        if self.method == "kmeans":

            self.model = KMeans(
                n_clusters=self.n_clusters,
                random_state=self.random_state,
                n_init=20
            )

            self.labels_ = self.model.fit_predict(X)

        # -------------------------
        # Hierarchical clustering
        # -------------------------
        elif self.method == "hierarchical":

            self.model = AgglomerativeClustering(
                n_clusters=self.n_clusters
            )

            self.labels_ = self.model.fit_predict(X)

        # -------------------------
        # Gaussian Mixture Model
        # -------------------------
        elif self.method == "gmm":

            self.model = GaussianMixture(
                n_components=self.n_clusters,
                covariance_type="full",
                random_state=self.random_state,
                n_init=10
            )

            self.model.fit(X)

            self.labels_ = self.model.predict(X)

        else:
            raise ValueError(
                "Method must be "
                "'kmeans', 'hierarchical', or 'gmm'"
            )

        return self

    def get_labels(self):
        return self.labels_

    def predict_proba(self, X: pd.DataFrame):
        """
        Returns soft cluster probabilities for GMM.
        Useful for biological uncertainty estimation.
        """

        if self.method != "gmm":
            raise ValueError(
                "predict_proba only available for GMM"
            )

        return self.model.predict_proba(X)

    def evaluate(self, X: pd.DataFrame):
        """
        Core clustering evaluation metrics.
        """

        if self.labels_ is None:
            return None

        if len(set(self.labels_)) < 2:
            return None

        silhouette = silhouette_score(X, self.labels_)

        db_score = davies_bouldin_score(X, self.labels_)

        results = {
            "method": self.method,
            "n_clusters": self.n_clusters,
            "silhouette_score": round(silhouette, 4),
            "davies_bouldin_score": round(db_score, 4)
        }

        # -------------------------
        # Add GMM-specific metrics
        # -------------------------
        if self.method == "gmm":

            results["bic"] = round(self.model.bic(X), 2)
            results["aic"] = round(self.model.aic(X), 2)

        return results