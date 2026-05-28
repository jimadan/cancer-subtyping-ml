import umap
from sklearn.metrics import silhouette_score
from src.models.clustering import ClusteringModel


def run_clustering_pipeline(X_pca):
    """
    Clustering pipeline:
    - UMAP for visualization only
    - PCA space used for clustering
    - k selection using silhouette score
    """

    # -------------------------
    # UMAP (VISUALIZATION ONLY)
    # -------------------------
    print("\n[STEP] UMAP...")

    reducer = umap.UMAP(
        n_neighbors=15,
        min_dist=0.1,
        n_components=2,   # IMPORTANT: 2D for plots
        random_state=42
    )

    X_umap = reducer.fit_transform(X_pca)

    # -------------------------
    # CLUSTERING + K SELECTION
    # -------------------------
    print("\n[STEP] Clustering (k search 2–5)...")

    best_k = None
    best_score = -1
    best_model = None
    best_labels = None

    results = []

    for k in range(2, 6):

        print(f"[CLUSTERING] Testing k={k}")

        model = ClusteringModel(method="gmm", n_clusters=k)
        model.fit(X_pca)

        labels = model.get_labels()

        score = silhouette_score(X_pca, labels)

        results.append((k, score))

        if score > best_score:
            best_k = k
            best_score = score
            best_model = model
            best_labels = labels

    print("[INFO] Best k:", best_k)
    print("[INFO] Silhouette:", best_score)

    # -------------------------
    # FINAL MODEL
    # -------------------------
    labels = best_labels

    eval_results = best_model.evaluate(X_pca)

    print("[INFO] Clustering results:", eval_results)

    return X_umap, labels, best_k, best_model