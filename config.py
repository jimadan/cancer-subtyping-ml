from pathlib import Path

# ---------------------------------------
# Project root (stable reference)
# ---------------------------------------

BASE_DIR = Path(__file__).resolve().parent

# ---------------------------------------
# Data paths
# ---------------------------------------

DATA_DIR = BASE_DIR / "data"
RAW_DIR = DATA_DIR / "raw"

EXPRESSION_FILE = RAW_DIR / "HiSeqV2.txt"
CLINICAL_FILE = RAW_DIR / "BRCA_survival.txt"

# ---------------------------------------
# Output paths
# ---------------------------------------

RESULTS_DIR = BASE_DIR / "results"
FIGURES_DIR = RESULTS_DIR / "figures"
GSEA_DIR = RESULTS_DIR / "gsea"
AGING_DIR = RESULTS_DIR / "aging"