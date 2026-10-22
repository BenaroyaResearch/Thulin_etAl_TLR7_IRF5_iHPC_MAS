# Figure 7D-F / S5C-D — Gene Set Enrichment pipeline

This folder's original scripts were not saved anywhere accessible (no `.py`/`.ipynb`
files or Claude Code session history existed for this analysis), so this code is a
**reconstruction**, reverse-engineered from the structure and values already present
in `GSEA_Results.xlsx`, `ORA_Results.xlsx`, `ORA_Clustered.xlsx`, and
`ORA_IRF5_Comparison_Summary.xlsx`. Each step below was validated against those
existing files; see "Validation" for exactly how closely each stage reproduces them.

## Pipeline

Input: `Q1 and Q3-IRF.xlsx` (RNA-seq/ATAC-seq DAR-DEG overlap, already labeled with
quadrant calls — Quadrant I = concordant RNA-up/ATAC-up genes = "Q1, up in iHPC";
Quadrant III = concordant RNA-down/ATAC-down genes = "Q3, up in TLR7 monocytes").

| Step | Script | Output |
|---|---|---|
| 1 | `01_prepare_gene_lists.py` | Gene lists + ranked (`.rnk`) lists for Q1/Q3 × All-genes/IRF5-bound, under `intermediate/` |
| 2 | `02_run_ora.py` | `ORA_Results.xlsx` — Enrichr ORA (MSigDB_Hallmark_2020, GO_Biological_Process_2021, KEGG_2019_Mouse), significant terms only (Adj. P < 0.05) |
| 3 | `03_run_gsea.py` | `GSEA_Results.xlsx` — prerank GSEA on the same gene universes/libraries; ES/NES merged onto the ORA-significant terms |
| 4 | `04_consistency_analysis.py` | Appends a "Consistency" sheet per quadrant (shared/unique significant pathways, Jaccard index) to both workbooks |
| 5 | `05_cluster_ora_pathways.py` | `ORA_Clustered.xlsx` — redundant ORA terms grouped by gene-overlap similarity; auto-drafted "For plot" curated shortlist |
| 6 | `06_irf5_comparison_summary.py` | `ORA_IRF5_Comparison_Summary.xlsx` — All-genes vs. IRF5-bound side-by-side stats for the curated "For plot" list |
| 7 | `07_make_figures.py` | All `Q1_iHPC_*` / `Q3_TLR7Mo_*` PNG/PDF figures |

Run everything with `python run_all.py` (needs internet access — step 2/3 query the
live Enrichr API via `gseapy`). **This overwrites the existing output files** —
copy them elsewhere first if you want to keep the current versions.

Mouse gene symbols are upper-cased before querying Enrichr, the standard
approximation for matching mouse genes to human-symbol Enrichr libraries (mouse and
human orthologs usually differ only in case).

## Manual curation step

The "For plot" sheet in `ORA_Clustered.xlsx` is a hand-picked shortlist of
non-redundant pathways chosen for the final figure — it is not fully automatable.
`05_cluster_ora_pathways.py` auto-drafts it (top cluster representatives by
Adj. P-value) as a starting point; **edit that sheet by hand** to match the
intended published panel, then re-run steps 6 and 7.

## Validation against the existing output files

- **ORA (`02_run_ora.py`)**: exact match. Re-running Enrichr on the reconstructed
  Q1/Q3 All-genes/IRF5-bound gene lists reproduced the existing `ORA_Results.xlsx`
  exactly — same significant-term counts (156 / 241 / 61 / 6) and identical
  Overlap, Fold Enrichment, Odds Ratio, Combined Score values.
- **Consistency (`04_consistency_analysis.py`)**: exact match, including the
  "% of total pathways shared in IRF5 analysis" metric, which is
  `shared / significant-in-All-genes` (not `/ significant-in-IRF5-bound`,
  despite the label).
- **IRF5 comparison summary (`06_irf5_comparison_summary.py`)**: exact match against
  `ORA_IRF5_Comparison_Summary.xlsx` for every pathway/stat checked.
- **Figures (`07_make_figures.py`)**: visually reproduces the lollipop-plot style
  (viridis color = -log10(adj. p), point size = gene count, black dot at x=0 for
  non-significant terms) and the ORA_plot dot-plot variant (Reds colormap, no stems).
- **GSEA (`03_run_gsea.py`)**: methodology-level match only, not bit-exact. The
  source `rnaseq_log2FC` values are rounded to 2 decimals, so a large fraction of
  genes tie in rank (~55% for Q1). GSEA's ES/NES depend on how ties are broken,
  and the original tie-break order wasn't recorded, so exact ES/NES values will
  differ slightly from a prior run. Pathway membership and the carried-over ORA
  stats (Overlap, Fold Enrichment) are unaffected.
- **Clustering (`05_cluster_ora_pathways.py`)**: approximate. The exact original
  clustering algorithm/threshold isn't recorded. Complete-linkage clustering on
  overlap-coefficient similarity (threshold 0.44, see `config.py`) reproduces the
  correct cluster *count* (79 for Q1 iHPC All Genes) and a similar size
  distribution (max cluster size 6, matching the original), but exact cluster
  membership can differ. Treat clusters as a redundancy-reduction aid for
  choosing the "For plot" shortlist, not as ground truth.

## Assumptions worth knowing about

- **Background gene universe** for Fold Enrichment back-calculation
  (`BACKGROUND_SIZE` in `config.py`) is fixed at 20,000 (Enrichr's typical default).
- **Enrichr organism** is `"human"` for all three libraries (KEGG_2019_Mouse is
  hosted within Enrichr's standard/human library set, not a separate mouse portal).
- **GSEA universe**: each quadrant/analysis is preranked using only its own gene
  pool (e.g., the 836 Quadrant-I genes for "Q1 All Genes"), not the whole
  transcriptome — inferred because Overlap/Fold Enrichment values in
  `GSEA_Results.xlsx` are identical to the corresponding `ORA_Results.xlsx` row.
