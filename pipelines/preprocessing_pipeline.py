from pathlib import Path

import pandas as pd

from src.data.loader.data_loader import load_expression, load_clinical
from src.data.expression_preprocess import ExpressionPreprocessor
from src.data.clinical_preprocess import ClinicalPreprocessor
from src.data.alignment import align_datasets


def run_preprocessing(cache=True):

    out_dir = Path("data/processed")
    expr_path = out_dir / "expression_processed.csv"
    clinical_path = out_dir / "clinical_processed.csv"

    # -------------------------
    # CHECK CACHE FIRST
    # -------------------------
    if cache and expr_path.exists() and clinical_path.exists():

        print("\n[INFO] Loading cached preprocessed data...")

        expr = pd.read_csv(expr_path, index_col=0)
        clinical = pd.read_csv(clinical_path, index_col=0)

        print("[INFO] Cached expression shape:", expr.shape)
        print("[INFO] Cached clinical shape:", clinical.shape)

        return expr, clinical

    # -------------------------
    # OTHERWISE RUN PIPELINE
    # -------------------------
    print("\n[STEP] Loading raw data...")

    expr_raw = load_expression("data/raw/HiSeqV2.txt")
    clinical_raw = load_clinical("data/raw/BRCA_survival.txt")

    print("\n[STEP] Preprocessing...")

    expr_processor = ExpressionPreprocessor(log_transform=True, min_variance=0.0)
    clinical_processor = ClinicalPreprocessor()

    expr = expr_processor.fit_transform(expr_raw)
    clinical = clinical_processor.transform(clinical_raw)

    print("\n[STEP] Aligning datasets...")

    expr, clinical = align_datasets(expr, clinical)

    # -------------------------
    # SAVE CACHE
    # -------------------------
    if cache:
        out_dir.mkdir(parents=True, exist_ok=True)

        print("\n[INFO] Saving processed data to cache...")

        expr.to_csv(expr_path)
        clinical.to_csv(clinical_path)

    return expr, clinical
