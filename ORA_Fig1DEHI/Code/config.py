"""Shared paths and parameters for the Figure 7D-F / S5C-D gene set enrichment pipeline."""
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent          # "Figure 7D-F_S5C-D_Gene Set Enrichment"
CODE_DIR = Path(__file__).resolve().parent
INTERMEDIATE_DIR = CODE_DIR / "intermediate"
GENE_LIST_DIR = INTERMEDIATE_DIR / "gene_lists"
RANKED_DIR = INTERMEDIATE_DIR / "ranked"
ORA_RAW_DIR = INTERMEDIATE_DIR / "ora_raw"

SOURCE_XLSX = ROOT / "Q1 and Q3-IRF.xlsx"
ALL_GENES_SHEET = "DAR-DEG overlap"
IRF5_BOUND_SHEET = "IRF5 Binding DAR-DEG overlap"

ORA_RESULTS_XLSX = ROOT / "ORA_Results.xlsx"
GSEA_RESULTS_XLSX = ROOT / "GSEA_Results.xlsx"
ORA_CLUSTERED_XLSX = ROOT / "ORA_Clustered.xlsx"
IRF5_COMPARISON_XLSX = ROOT / "ORA_IRF5_Comparison_Summary.xlsx"

# Quadrant I = concordant RNA-up/ATAC-up genes ("Q1", up in iHPC).
# Quadrant III = concordant RNA-down/ATAC-down genes ("Q3", up in TLR7 monocytes).
QUADRANTS = {
    "Q1": {"code": "I", "label": "iHPC", "file_tag": "Q1_iHPC"},
    "Q3": {"code": "III", "label": "TLR7 Mo", "file_tag": "Q3_TLR7Mo"},
}
ANALYSES = ["all", "irf5"]  # "all" genes vs. IRF5-bound subset
ANALYSIS_LABEL = {"all": "All Genes", "irf5": "IRF5 Bound"}

GENE_SET_LIBRARIES = [
    "MSigDB_Hallmark_2020",
    "GO_Biological_Process_2021",
    "KEGG_2019_Mouse",
]

# Mouse gene symbols are upper-cased to match human-symbol Enrichr libraries
# (standard approximation for mouse/human ortholog symbols that differ only in case).
ENRICHR_ORGANISM = "human"

SIG_ADJP = 0.05
# Enrichr's default background gene universe size, used to back-calculate
# Fold Enrichment = (k/n) / (K/N) from the reported Overlap ("k/K") string.
BACKGROUND_SIZE = 20000

GSEA_PERMUTATIONS = 1000
GSEA_MIN_SIZE = 5
GSEA_MAX_SIZE = 1000
GSEA_SEED = 0

# Overlap-coefficient (|A∩B| / min(|A|,|B|), on ORA overlap-gene sets)
# threshold for merging two significant ORA terms into the same redundancy
# cluster. A term joins a cluster only if it clears this threshold against
# EVERY existing member (complete linkage). 0.44 best reproduces the cluster
# count/size distribution seen in the original ORA_Clustered.xlsx, but the
# original clustering procedure was not recorded, so exact membership can
# differ -- treat clusters as a redundancy-reduction aid, not ground truth.
CLUSTER_SIMILARITY_THRESHOLD = 0.44

for d in (INTERMEDIATE_DIR, GENE_LIST_DIR, RANKED_DIR, ORA_RAW_DIR):
    d.mkdir(parents=True, exist_ok=True)


def sheet_name(quadrant_key, analysis):
    q = QUADRANTS[quadrant_key]
    return f"{quadrant_key} {q['label']} - {ANALYSIS_LABEL[analysis]}"


def key_name(quadrant_key, analysis):
    return f"{quadrant_key}_{analysis}"
