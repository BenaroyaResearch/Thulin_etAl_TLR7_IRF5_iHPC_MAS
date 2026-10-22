"""
This code was generated using Claude Code, Anthropic (Opus 5.5) and carefully by a bioinformatician.

Fig. 2C: kinetics of R848-induced, IRF5-dependent genes in WT and IRF5ko
bone marrow (BM) monocytes in vitro.

Overview
    The heatmap shows the subset of genes that were (i) R848-induced in an
    IRF5-dependent manner in BM monocytes (Fig. 2B, Fig. S3C-D) and (ii) upregulated with
    increased chromatin accessibility in iHPCs compared to Ly6Chi monocytes
    (Fig. 1C, quadrant 1, green dots). These 92 genes are the overlap of the
    two gene lists in the Venn diagram of Fig. 2B.

Inputs (Data/)
    - Data_VennDiagram_Fig2B.xlsx
        sheet "Venn Overlap_Fig 2B": one column of gene symbols. The first row
        of the sheet is a free-text description; the column header
        ("Venn Overlap") is on the second row.
    - P689-1_InVitro_R848_WT-IRF5ko_Mo_normalized log2 expression.csv
        normalized log2 expression, genes x samples. First column "gene";
        sample columns are named <genotype>_<treatment>_<time>_<library>,
        e.g. WT_R848_6h_lib112676.

Method
    Samples: WT and IRF5ko monocytes, untreated at 0h and R848-treated for
    6h, 24h and 48h. Each gene (row) is z-scored across these samples and rows
    are hierarchically clustered (seaborn defaults: average linkage, Euclidean
    distance). Columns are ordered by genotype, then time (WT 0h/6h/24h/48h,
    then IRF5ko 0h/6h/24h/48h), are not clustered, and carry two color
    annotation bars (genotype, time).

Output (Generated Files/)
    - Fig2C_VennOverlap_heatmap.pdf

Usage
    python Code/01_generate_venn_overlap_heatmap.py
    (paths are resolved relative to this script, so it can be run from any
    working directory as long as the Code/, Data/ and Generated Files/
    folders stay side by side)
"""
import matplotlib
matplotlib.rcParams["pdf.fonttype"] = 42   # embed text as editable TrueType, not Type 3 bitmap
matplotlib.rcParams["ps.fonttype"] = 42

import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ---- paths (relative to this script: <project>/Code/<this file>) ----
PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(PROJECT_DIR, "Data")
OUT_DIR = os.path.join(PROJECT_DIR, "Generated Files")

VENN_PATH = os.path.join(DATA_DIR, "Data_VennDiagram_Fig2B.xlsx")
VENN_SHEET = "Venn Overlap_Fig 2B"
VENN_COLUMN = "Venn Overlap"
EXPR_PATH = os.path.join(DATA_DIR, "P689-1_InVitro_R848_WT-IRF5ko_Mo_normalized log2 expression.csv")
OUT_PATH = os.path.join(OUT_DIR, "Fig2C_VennOverlap_heatmap.pdf")

# ---- plot colors ----
COLOR_WT = "#1b9e77"      # teal
COLOR_IRF5KO = "#d95f02"  # orange
TIME_COLORS = {"0h": "#f7fbff", "6h": "#deebf7", "24h": "#9ecae1", "48h": "#3182bd"}  # light -> dark blue

GENOTYPE_ORDER = {"WT": 0, "IRF5ko": 1}
TIME_ORDER = {"0h": 0, "6h": 1, "24h": 2, "48h": 3}

# ---- 1. gene list from the Venn overlap sheet ----
# skiprows=1 skips the description row so the "Venn Overlap" header is read as the column name
venn = pd.read_excel(VENN_PATH, sheet_name=VENN_SHEET, skiprows=1)
genes = venn[VENN_COLUMN].dropna().astype(str).str.strip().tolist()
print(f"{len(genes)} genes read from '{VENN_SHEET}'")

# ---- 2. expression matrix: 0h untreated + R848-treated 6h/24h/48h samples ----
expr = pd.read_csv(EXPR_PATH).set_index("gene")

def sample_info(col):
    """Split a sample name like 'WT_R848_6h_lib112676' into (genotype, treatment, time)."""
    genotype, treatment, time = col.split("_")[:3]
    return genotype, treatment, time

sample_cols = [c for c in expr.columns if sample_info(c)[2] in TIME_ORDER]
sample_cols.sort(key=lambda c: (GENOTYPE_ORDER[sample_info(c)[0]], TIME_ORDER[sample_info(c)[2]]))

print(f"{len(sample_cols)} samples selected:")
for c in sample_cols:
    print(" ", c)

present = [g for g in genes if g in expr.index]
missing = [g for g in genes if g not in expr.index]
if missing:
    print(f"{len(missing)} gene(s) not found in expression matrix, dropped: {missing}")

sub = expr.loc[present, sample_cols]
zero_var = sub.index[sub.std(axis=1) == 0]
if len(zero_var):
    print(f"{len(zero_var)} gene(s) with zero variance across selected samples, dropped: {list(zero_var)}")
    sub = sub.drop(index=zero_var)

# ---- 3. row z-score (per gene, across the selected samples) ----
z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1), axis=0)

# ---- 4. column annotation bars: genotype + time ----
genotypes = [sample_info(c)[0] for c in sample_cols]
times = [sample_info(c)[2] for c in sample_cols]
col_colors = pd.DataFrame({
    "Genotype": [COLOR_WT if g == "WT" else COLOR_IRF5KO for g in genotypes],
    "Time": [TIME_COLORS[t] for t in times],
}, index=sample_cols)

# ---- 5. clustered heatmap ----
n_genes = z.shape[0]
fig_h = max(6, min(40, n_genes * 0.16 + 2))   # scale height with gene count

g = sns.clustermap(
    z,
    row_cluster=True,
    col_cluster=False,
    cmap="vlag",
    center=0,
    robust=True,          # color limits set from the 2nd-98th percentile of z-scores
    col_colors=col_colors,
    xticklabels=False,
    yticklabels=True,
    figsize=(9, fig_h),
    dendrogram_ratio=(0.12, 0.02),
    cbar_pos=(0.02, 0.83, 0.03, 0.15),
)
g.ax_heatmap.set_ylabel("")
g.ax_heatmap.set_xlabel("")
g.ax_heatmap.set_yticklabels(list(z.index[g.dendrogram_row.reordered_ind]), fontsize=6)
g.ax_cbar.set_title("Row z-score", fontsize=7)

# ---- 6. legends for the annotation bars ----
genotype_handles = [mpatches.Patch(color=COLOR_WT, label="WT"),
                    mpatches.Patch(color=COLOR_IRF5KO, label="IRF5ko")]
leg1 = g.fig.legend(handles=genotype_handles, title="Genotype",
                    bbox_to_anchor=(1.0, 0.95), loc="upper left",
                    bbox_transform=g.fig.transFigure,
                    frameon=False, fontsize=7, title_fontsize=7)
g.fig.add_artist(leg1)

time_handles = [mpatches.Patch(color=TIME_COLORS[t], label=t) for t in TIME_ORDER]
leg2 = g.fig.legend(handles=time_handles, title="Time",
                    bbox_to_anchor=(1.0, 0.80), loc="upper left",
                    bbox_transform=g.fig.transFigure,
                    frameon=False, fontsize=7, title_fontsize=7)
g.fig.add_artist(leg2)

g.fig.suptitle(f"Venn Overlap genes (n={n_genes}), 0h/6h/24h/48h, WT vs IRF5ko monocytes", y=1.02, fontsize=11)

# ---- 7. save ----
os.makedirs(OUT_DIR, exist_ok=True)
g.fig.savefig(OUT_PATH, bbox_inches="tight")
plt.close(g.fig)
print(f"wrote {OUT_PATH}")
