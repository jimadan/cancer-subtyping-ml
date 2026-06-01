from src.features.feature_filtering import FeatureFilter
from src.features.scaling import ExpressionScaler
from src.features.dimensionality_reduction import PCAReducer
from src.features.cox_feature_selection import cox_feature_selection


def run_feature_pipeline(expr, clinical):

    print("\n[STEP] Feature filtering...")

    filter_model = FeatureFilter(method="variance", min_variance=0.0, top_n_genes=5000) #0.5, 200
    expr_filtered = filter_model.fit_transform(expr)

    print("\n[STEP] Survival feature selection...")

    expr_survival, ranked_genes = cox_feature_selection(
        expr_filtered,
        clinical,
        top_n=500
    )

    print("[INFO] Cox-selected expression shape:", expr_survival.shape)

    print("\n[STEP] Scaling...")

    scaler = ExpressionScaler()
    X_scaled = scaler.fit_transform(expr_survival)

    print("\n[STEP] PCA...")

    pca = PCAReducer(n_components=20)
    X_pca = pca.fit_transform(X_scaled)

    return expr_filtered, expr_survival, X_pca
